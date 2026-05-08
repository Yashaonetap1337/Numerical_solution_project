from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Optional
GAMMA = 1.4
NGHOST = 2

def cons_to_prim(U: np.ndarray) -> np.ndarray:
    rho = U[0]
    u = U[1] / rho
    v = U[2] / rho
    E = U[3]
    p = (GAMMA - 1.0) * (E - 0.5 * rho * (u * u + v * v))
    W = np.empty_like(U)
    W[0] = rho
    W[1] = u
    W[2] = v
    W[3] = p
    return W

def prim_to_cons(W: np.ndarray) -> np.ndarray:
    rho, u, v, p = (W[0], W[1], W[2], W[3])
    U = np.empty_like(W)
    U[0] = rho
    U[1] = rho * u
    U[2] = rho * v
    U[3] = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    return U

def sound_speed(rho: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.sqrt(GAMMA * np.maximum(p, 1e-12) / np.maximum(rho, 1e-12))

def flux_x(W: np.ndarray) -> np.ndarray:
    rho, u, v, p = (W[0], W[1], W[2], W[3])
    E = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    F = np.empty_like(W)
    F[0] = rho * u
    F[1] = rho * u * u + p
    F[2] = rho * u * v
    F[3] = (E + p) * u
    return F

def flux_y(W: np.ndarray) -> np.ndarray:
    rho, u, v, p = (W[0], W[1], W[2], W[3])
    E = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    G = np.empty_like(W)
    G[0] = rho * v
    G[1] = rho * u * v
    G[2] = rho * v * v + p
    G[3] = (E + p) * v
    return G

def hllc_flux_x(WL: np.ndarray, WR: np.ndarray) -> np.ndarray:
    rL, uL, vL, pL = (WL[0], WL[1], WL[2], WL[3])
    rR, uR, vR, pR = (WR[0], WR[1], WR[2], WR[3])
    aL = sound_speed(rL, pL)
    aR = sound_speed(rR, pR)
    SL = np.minimum(uL - aL, uR - aR)
    SR = np.maximum(uL + aL, uR + aR)
    num = pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)
    den = rL * (SL - uL) - rR * (SR - uR)
    S_star = num / np.where(np.abs(den) < 1e-30, 1e-30, den)
    EL = pL / (GAMMA - 1.0) + 0.5 * rL * (uL * uL + vL * vL)
    ER = pR / (GAMMA - 1.0) + 0.5 * rR * (uR * uR + vR * vR)
    FL = flux_x(WL)
    FR = flux_x(WR)
    UL = np.empty_like(WL)
    UL[0] = rL
    UL[1] = rL * uL
    UL[2] = rL * vL
    UL[3] = EL
    UR = np.empty_like(WR)
    UR[0] = rR
    UR[1] = rR * uR
    UR[2] = rR * vR
    UR[3] = ER
    factor_L = rL * (SL - uL) / np.where(np.abs(SL - S_star) < 1e-30, 1e-30, SL - S_star)
    factor_R = rR * (SR - uR) / np.where(np.abs(SR - S_star) < 1e-30, 1e-30, SR - S_star)
    UstarL = np.empty_like(WL)
    UstarL[0] = factor_L
    UstarL[1] = factor_L * S_star
    UstarL[2] = factor_L * vL
    UstarL[3] = factor_L * (EL / rL + (S_star - uL) * (S_star + pL / (rL * (SL - uL))))
    UstarR = np.empty_like(WR)
    UstarR[0] = factor_R
    UstarR[1] = factor_R * S_star
    UstarR[2] = factor_R * vR
    UstarR[3] = factor_R * (ER / rR + (S_star - uR) * (S_star + pR / (rR * (SR - uR))))
    F = np.empty_like(WL)
    m1 = SL >= 0.0
    m2 = (SL < 0.0) & (S_star >= 0.0)
    m3 = (S_star < 0.0) & (SR >= 0.0)
    m4 = SR < 0.0
    for k in range(4):
        F[k] = np.where(m1, FL[k], np.where(m2, FL[k] + SL * (UstarL[k] - UL[k]), np.where(m3, FR[k] + SR * (UstarR[k] - UR[k]), FR[k])))
    return F

def hllc_flux_y(WL: np.ndarray, WR: np.ndarray) -> np.ndarray:
    WL_swap = np.empty_like(WL)
    WL_swap[0] = WL[0]
    WL_swap[1] = WL[2]
    WL_swap[2] = WL[1]
    WL_swap[3] = WL[3]
    WR_swap = np.empty_like(WR)
    WR_swap[0] = WR[0]
    WR_swap[1] = WR[2]
    WR_swap[2] = WR[1]
    WR_swap[3] = WR[3]
    F_swap = hllc_flux_x(WL_swap, WR_swap)
    F = np.empty_like(F_swap)
    F[0] = F_swap[0]
    F[1] = F_swap[2]
    F[2] = F_swap[1]
    F[3] = F_swap[3]
    return F

def minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))

def reconstruct_x(W: np.ndarray):
    dW_left = W[:, :, 1:-1] - W[:, :, :-2]
    dW_right = W[:, :, 2:] - W[:, :, 1:-1]
    slope = minmod(dW_left, dW_right)
    WL = W[:, :, 1:-2] + 0.5 * slope[:, :, :-1]
    WR = W[:, :, 2:-1] - 0.5 * slope[:, :, 1:]
    return (WL, WR)

def reconstruct_y(W: np.ndarray):
    dW_down = W[:, 1:-1, :] - W[:, :-2, :]
    dW_up = W[:, 2:, :] - W[:, 1:-1, :]
    slope = minmod(dW_down, dW_up)
    WL = W[:, 1:-2, :] + 0.5 * slope[:, :-1, :]
    WR = W[:, 2:-1, :] - 0.5 * slope[:, 1:, :]
    return (WL, WR)

def apply_bc(U: np.ndarray, bc: dict) -> None:
    g = NGHOST
    if bc['left'] == 'wall':
        for k in range(g):
            U[:, :, g - 1 - k] = U[:, :, g + k]
            U[1, :, g - 1 - k] = -U[1, :, g + k]
    elif bc['left'] == 'outflow':
        for k in range(g):
            U[:, :, g - 1 - k] = U[:, :, g]
    elif bc['left'] == 'inflow_post_shock':
        Wstate = bc['inflow_state']
        Ustate = prim_to_cons(Wstate.reshape(4, 1, 1)).reshape(4)
        for k in range(g):
            for c in range(4):
                U[c, :, g - 1 - k] = Ustate[c]
    elif bc['left'] == 'steady_inflow':
        state = bc['steady_inflow_state']
        Wstate = np.array([state['rho_0'], state['u_0'], state['v_0'], state['p_0']])
        Ustate = prim_to_cons(Wstate.reshape(4, 1, 1)).reshape(4)
        for k in range(g):
            for c in range(4):
                U[c, :, g - 1 - k] = Ustate[c]
    elif bc['left'] == 'steady_profile_inflow':
        state = bc['left_state']
        y_shock = bc['left_y_shock']
        y_coords = bc['y_centers']
        ny_inner = len(y_coords)
        W0 = np.array([state['rho_0'], state['u_0'], state['v_0'], state['p_0']])
        W1 = np.array([state['rho_1'], state['u_1'], state['v_1'], state['p_1']])
        U0 = prim_to_cons(W0.reshape(4, 1, 1)).reshape(4)
        U1 = prim_to_cons(W1.reshape(4, 1, 1)).reshape(4)
        for j in range(ny_inner):
            j_full = j + g
            Ufill = U1 if y_coords[j] < y_shock else U0
            for k in range(g):
                for c in range(4):
                    U[c, j_full, g - 1 - k] = Ufill[c]
    if bc['right'] == 'wall':
        nx = U.shape[2]
        for k in range(g):
            U[:, :, nx - g + k] = U[:, :, nx - g - 1 - k]
            U[1, :, nx - g + k] = -U[1, :, nx - g - 1 - k]
    elif bc['right'] == 'outflow':
        nx = U.shape[2]
        for k in range(g):
            U[:, :, nx - g + k] = U[:, :, nx - g - 1]
    if bc['bottom'] == 'wall':
        for k in range(g):
            U[:, g - 1 - k, :] = U[:, g + k, :]
            U[2, g - 1 - k, :] = -U[2, g + k, :]
    elif bc['bottom'] == 'outflow':
        for k in range(g):
            U[:, g - 1 - k, :] = U[:, g, :]
    if bc['top'] == 'wall':
        ny = U.shape[1]
        for k in range(g):
            U[:, ny - g + k, :] = U[:, ny - g - 1 - k, :]
            U[2, ny - g + k, :] = -U[2, ny - g - 1 - k, :]
    elif bc['top'] == 'outflow':
        ny = U.shape[1]
        for k in range(g):
            U[:, ny - g + k, :] = U[:, ny - g - 1, :]

def compute_dt(U: np.ndarray, dx: float, dy: float, cfl: float) -> float:
    W = cons_to_prim(U)
    rho, u, v, p = (W[0], W[1], W[2], W[3])
    a = sound_speed(rho, p)
    smax_x = np.max(np.abs(u) + a)
    smax_y = np.max(np.abs(v) + a)
    dt = cfl / (smax_x / dx + smax_y / dy + 1e-30)
    return dt

def rhs_muscl(U: np.ndarray, dx: float, dy: float, bc: dict) -> np.ndarray:
    apply_bc(U, bc)
    W = cons_to_prim(U)
    WLx, WRx = reconstruct_x(W)
    Fx = hllc_flux_x(WLx, WRx)
    WLy, WRy = reconstruct_y(W)
    Gy = hllc_flux_y(WLy, WRy)
    rhs = np.zeros_like(U)
    rhs[:, 2:-2, 2:-2] = -(Fx[:, 2:-2, 1:] - Fx[:, 2:-2, :-1]) / dx - (Gy[:, 1:, 2:-2] - Gy[:, :-1, 2:-2]) / dy
    return rhs

def apply_wedge_bc(U: np.ndarray, wedge_mask: np.ndarray, wedge_normal_x: np.ndarray, wedge_normal_y: np.ndarray, rho_gas_default: float=1.0, p_gas_default: float=1.0, gamma: float=GAMMA) -> None:
    NY, NX = wedge_mask.shape
    if not np.any(wedge_mask):
        return
    U[0][wedge_mask] = rho_gas_default
    U[1][wedge_mask] = 0.0
    U[2][wedge_mask] = 0.0
    U[3][wedge_mask] = p_gas_default / (gamma - 1.0)
    di = np.round(wedge_normal_x).astype(int)
    dj = np.round(wedge_normal_y).astype(int)
    j_idx, i_idx = np.where(wedge_mask)
    if len(j_idx) == 0:
        return
    j_neigh = j_idx + dj[j_idx, i_idx]
    i_neigh = i_idx + di[j_idx, i_idx]
    valid = (j_neigh >= 0) & (j_neigh < NY) & (i_neigh >= 0) & (i_neigh < NX)
    valid &= ~wedge_mask[np.clip(j_neigh, 0, NY - 1), np.clip(i_neigh, 0, NX - 1)]
    j_w = j_idx[valid]
    i_w = i_idx[valid]
    j_n = j_neigh[valid]
    i_n = i_neigh[valid]
    if len(j_w) == 0:
        return
    nx_n = wedge_normal_x[j_w, i_w]
    ny_n = wedge_normal_y[j_w, i_w]
    rho = U[0, j_n, i_n]
    u = U[1, j_n, i_n] / rho
    v = U[2, j_n, i_n] / rho
    E = U[3, j_n, i_n]
    v_n = u * nx_n + v * ny_n
    u_new = u - 2.0 * v_n * nx_n
    v_new = v - 2.0 * v_n * ny_n
    U[0, j_w, i_w] = rho
    U[1, j_w, i_w] = rho * u_new
    U[2, j_w, i_w] = rho * v_new
    U[3, j_w, i_w] = E

def step_ssprk2(U: np.ndarray, dt: float, dx: float, dy: float, bc: dict, wedge_mask: Optional[np.ndarray]=None, wedge_normal_x: Optional[np.ndarray]=None, wedge_normal_y: Optional[np.ndarray]=None) -> np.ndarray:
    L1 = rhs_muscl(U, dx, dy, bc)
    U1 = U + dt * L1
    if wedge_mask is not None:
        apply_wedge_bc(U1, wedge_mask, wedge_normal_x, wedge_normal_y)
    L2 = rhs_muscl(U1, dx, dy, bc)
    U_new = 0.5 * (U + U1 + dt * L2)
    if wedge_mask is not None:
        apply_wedge_bc(U_new, wedge_mask, wedge_normal_x, wedge_normal_y)
    return U_new

@dataclass
class EulerSim:
    nx: int
    ny: int
    Lx: float
    Ly: float
    cfl: float = 0.45
    bc: dict = field(default_factory=dict)

    def __post_init__(self):
        self.dx = self.Lx / self.nx
        self.dy = self.Ly / self.ny
        self.NX = self.nx + 2 * NGHOST
        self.NY = self.ny + 2 * NGHOST
        self.x = (np.arange(self.nx) + 0.5) * self.dx
        self.y = (np.arange(self.ny) + 0.5) * self.dy
        self.x_full = (np.arange(self.NX) - NGHOST + 0.5) * self.dx
        self.y_full = (np.arange(self.NY) - NGHOST + 0.5) * self.dy
        self.U = np.zeros((4, self.NY, self.NX))
        self.t = 0.0
        if isinstance(self.bc, dict) and self.bc.get('left') == 'steady_profile_inflow':
            self.bc['y_centers'] = self.y
        self.wedge_mask: Optional[np.ndarray] = None
        self.wedge_normal_x: Optional[np.ndarray] = None
        self.wedge_normal_y: Optional[np.ndarray] = None

    def set_wedge(self, x_corner: float, theta_w_deg: float, location: str='top', wedge_length: Optional[float]=None):
        theta_w = np.deg2rad(theta_w_deg)
        X, Y = np.meshgrid(self.x_full, self.y_full)
        if wedge_length is not None:
            x_back = x_corner + wedge_length
        else:
            x_back = self.Lx + 100
        if location == 'top':
            y_surf = self.Ly - np.tan(theta_w) * (X - x_corner)
            mask = (Y > y_surf) & (X > x_corner) & (X <= x_back)
            n_x = -np.sin(theta_w)
            n_y = -np.cos(theta_w)
            normal_x = np.full_like(X, n_x)
            normal_y = np.full_like(X, n_y)
        elif location == 'bottom':
            y_surf = np.tan(theta_w) * (X - x_corner)
            mask = (Y < y_surf) & (X > x_corner) & (X <= x_back)
            n_x = -np.sin(theta_w)
            n_y = np.cos(theta_w)
            normal_x = np.full_like(X, n_x)
            normal_y = np.full_like(X, n_y)
        else:
            raise ValueError(f"location must be 'top' or 'bottom', got {location!r}")
        self.wedge_mask = mask
        self.wedge_normal_x = normal_x
        self.wedge_normal_y = normal_y
        self.wedge_x_corner = x_corner
        self.wedge_theta_w_deg = theta_w_deg
        self.wedge_location = location
        self.wedge_length = wedge_length

    def set_initial_primitive(self, init_func: Callable):
        X, Y = np.meshgrid(self.x_full, self.y_full)
        W = np.zeros((4, self.NY, self.NX))
        rho, u, v, p = init_func(X, Y)
        W[0] = rho
        W[1] = u
        W[2] = v
        W[3] = p
        self.U = prim_to_cons(W)

    def get_primitive(self):
        g = NGHOST
        W = cons_to_prim(self.U[:, g:-g, g:-g])
        return W

    def advance(self, t_end: float, verbose_every: Optional[int]=None, callback: Optional[Callable]=None):
        n = 0
        while self.t < t_end:
            dt = compute_dt(self.U, self.dx, self.dy, self.cfl)
            if self.t + dt > t_end:
                dt = t_end - self.t
            self.U = step_ssprk2(self.U, dt, self.dx, self.dy, self.bc, wedge_mask=self.wedge_mask, wedge_normal_x=self.wedge_normal_x, wedge_normal_y=self.wedge_normal_y)
            self.t += dt
            n += 1
            if verbose_every is not None and n % verbose_every == 0:
                W = self.get_primitive()
                print(f'  step {n:5d}  t={self.t:.4f}  dt={dt:.2e}  rho∈[{W[0].min():.3f},{W[0].max():.3f}]')
            if callback is not None:
                callback(self, n)
        return n

    def advance_n_steps(self, n_steps: int, callback: Optional[Callable]=None):
        for n in range(n_steps):
            dt = compute_dt(self.U, self.dx, self.dy, self.cfl)
            self.U = step_ssprk2(self.U, dt, self.dx, self.dy, self.bc, wedge_mask=self.wedge_mask, wedge_normal_x=self.wedge_normal_x, wedge_normal_y=self.wedge_normal_y)
            self.t += dt
            if callback is not None:
                callback(self, n)
        return n_steps

    def update_bc(self, new_bc: dict):
        if new_bc.get('left') == 'steady_profile_inflow':
            new_bc = dict(new_bc)
            new_bc['y_centers'] = self.y
        self.bc = new_bc
