"""
Решатель 2D уравнений Эйлера для невязкого сжимаемого газа.

Численная схема:
    - Метод конечных объёмов на декартовой сетке.
    - MUSCL-Hancock второго порядка: линейная реконструкция + минимод-лимитер
      + полушаг по времени (predictor-corrector).
    - HLLC римановый решатель для численных потоков.
    - Управление по числу Куранта.

Консервативные переменные: U = (rho, rho*u, rho*v, E),
где E = rho*e + 0.5*rho*(u^2 + v^2),  p = (gamma-1)*rho*e.

Все массивы имеют форму (4, ny, nx). Индекс (j, i) = (y, x).
Граничные ячейки: используем 2 ghost-слоя с каждой стороны.
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Optional


GAMMA = 1.4  # показатель адиабаты для воздуха
NGHOST = 2   # число ghost-ячеек с каждой стороны (нужно для MUSCL)


# =============================================================================
# Преобразования между консервативными и примитивными переменными
# =============================================================================

def cons_to_prim(U: np.ndarray) -> np.ndarray:
    """U=(rho, rho*u, rho*v, E) -> W=(rho, u, v, p)."""
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
    """W=(rho, u, v, p) -> U=(rho, rho*u, rho*v, E)."""
    rho, u, v, p = W[0], W[1], W[2], W[3]
    U = np.empty_like(W)
    U[0] = rho
    U[1] = rho * u
    U[2] = rho * v
    U[3] = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    return U


def sound_speed(rho: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.sqrt(GAMMA * np.maximum(p, 1e-12) / np.maximum(rho, 1e-12))


# =============================================================================
# Физические потоки в декартовой геометрии
# =============================================================================

def flux_x(W: np.ndarray) -> np.ndarray:
    """Физический поток F(U) в направлении x по примитивам."""
    rho, u, v, p = W[0], W[1], W[2], W[3]
    E = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    F = np.empty_like(W)
    F[0] = rho * u
    F[1] = rho * u * u + p
    F[2] = rho * u * v
    F[3] = (E + p) * u
    return F


def flux_y(W: np.ndarray) -> np.ndarray:
    """Физический поток G(U) в направлении y по примитивам."""
    rho, u, v, p = W[0], W[1], W[2], W[3]
    E = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    G = np.empty_like(W)
    G[0] = rho * v
    G[1] = rho * u * v
    G[2] = rho * v * v + p
    G[3] = (E + p) * v
    return G


# =============================================================================
# HLLC римановый решатель
# =============================================================================

def hllc_flux_x(WL: np.ndarray, WR: np.ndarray) -> np.ndarray:
    """
    HLLC численный поток на грани, нормальной к оси x.
    WL, WR — примитивы слева/справа от грани, форма (4, ...).

    Реализация по Toro (2009), глава 10. Используем оценку Davis для волновых
    скоростей S_L, S_R. Контактная скорость S_star — формула Batten et al.
    """
    rL, uL, vL, pL = WL[0], WL[1], WL[2], WL[3]
    rR, uR, vR, pR = WR[0], WR[1], WR[2], WR[3]

    aL = sound_speed(rL, pL)
    aR = sound_speed(rR, pR)

    # Простая оценка волновых скоростей (Davis)
    SL = np.minimum(uL - aL, uR - aR)
    SR = np.maximum(uL + aL, uR + aR)

    # Контактная скорость
    num = pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)
    den = rL * (SL - uL) - rR * (SR - uR)
    S_star = num / np.where(np.abs(den) < 1e-30, 1e-30, den)

    # Полные энергии
    EL = pL / (GAMMA - 1.0) + 0.5 * rL * (uL * uL + vL * vL)
    ER = pR / (GAMMA - 1.0) + 0.5 * rR * (uR * uR + vR * vR)

    # Левый и правый физические потоки
    FL = flux_x(WL)
    FR = flux_x(WR)

    # Состояния U_L, U_R, U*_L, U*_R
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

    # U* состояния (Toro, eq. 10.39)
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

    # Поток HLLC: выбираем по знаку волновых скоростей
    F = np.empty_like(WL)
    # маски для четырёх режимов
    m1 = SL >= 0.0
    m2 = (SL < 0.0) & (S_star >= 0.0)
    m3 = (S_star < 0.0) & (SR >= 0.0)
    m4 = SR < 0.0

    for k in range(4):
        F[k] = np.where(m1, FL[k],
              np.where(m2, FL[k] + SL * (UstarL[k] - UL[k]),
              np.where(m3, FR[k] + SR * (UstarR[k] - UR[k]),
                            FR[k])))
    return F


def hllc_flux_y(WL: np.ndarray, WR: np.ndarray) -> np.ndarray:
    """HLLC поток в направлении y. Меняем местами роли u и v."""
    # Перестановка компонент: (rho, u, v, p) -> (rho, v, u, p) для x-формулы
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

    # обратно: (rho_flux, mom_y_flux, mom_x_flux, E_flux) -> y-поток
    F = np.empty_like(F_swap)
    F[0] = F_swap[0]
    F[1] = F_swap[2]
    F[2] = F_swap[1]
    F[3] = F_swap[3]
    return F


# =============================================================================
# MUSCL-реконструкция с minmod-лимитером
# =============================================================================

def minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Стандартный minmod лимитер."""
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))


def reconstruct_x(W: np.ndarray):
    """
    Линейная реконструкция в x-направлении.
    Возвращает левое и правое состояния на гранях между ячейками.
    Если W имеет форму (4, ny, nx), возвращаются:
        WL_face: (4, ny, nx-1) — состояние слева от грани (i+1/2)
        WR_face: (4, ny, nx-1) — состояние справа от грани (i+1/2)
    Грани нумеруются от 0 до nx-2.
    """
    # Разности
    dW_left = W[:, :, 1:-1] - W[:, :, :-2]    # W[i] - W[i-1]
    dW_right = W[:, :, 2:] - W[:, :, 1:-1]    # W[i+1] - W[i]
    slope = minmod(dW_left, dW_right)         # наклон в ячейке i (i от 1 до nx-2)

    # Состояния на гранях:
    # Левое состояние на грани (i+1/2) — экстраполяция из ячейки i вправо
    # Правое состояние на грани (i+1/2) — экстраполяция из ячейки i+1 влево
    # Грани между ячейками i и i+1, где i от 1 до nx-3 (т.е. внутренние грани)
    WL = W[:, :, 1:-2] + 0.5 * slope[:, :, :-1]   # из ячейки i
    WR = W[:, :, 2:-1] - 0.5 * slope[:, :, 1:]    # из ячейки i+1
    return WL, WR


def reconstruct_y(W: np.ndarray):
    """Линейная реконструкция в y-направлении."""
    dW_down = W[:, 1:-1, :] - W[:, :-2, :]
    dW_up = W[:, 2:, :] - W[:, 1:-1, :]
    slope = minmod(dW_down, dW_up)

    WL = W[:, 1:-2, :] + 0.5 * slope[:, :-1, :]
    WR = W[:, 2:-1, :] - 0.5 * slope[:, 1:, :]
    return WL, WR


# =============================================================================
# Граничные условия
# =============================================================================

def apply_bc(U: np.ndarray, bc: dict) -> None:
    """
    Применяет граничные условия путём заполнения ghost-слоёв.
    bc = {'left': type, 'right': type, 'bottom': type, 'top': type}
    типы: 'wall' (твёрдая стенка с отражением), 'outflow' (zero-gradient),
          'inflow_post_shock' (фиксированное состояние за УВ — подаётся через bc['inflow_state'])
    """
    g = NGHOST
    # Левая
    if bc['left'] == 'wall':
        # отражение: знак нормальной (x) компоненты импульса меняется
        for k in range(g):
            U[:, :, g - 1 - k] = U[:, :, g + k]
            U[1, :, g - 1 - k] = -U[1, :, g + k]
    elif bc['left'] == 'outflow':
        for k in range(g):
            U[:, :, g - 1 - k] = U[:, :, g]
    elif bc['left'] == 'inflow_post_shock':
        Wstate = bc['inflow_state']  # форма (4,)
        Ustate = prim_to_cons(Wstate.reshape(4, 1, 1)).reshape(4)
        for k in range(g):
            for c in range(4):
                U[c, :, g - 1 - k] = Ustate[c]
    elif bc['left'] == 'steady_inflow':
        # Стационарный сверхзвуковой набегающий поток: задаётся однородный 
        # профиль (rho_0, u_0, 0, p_0). Косая УВ возникает естественным образом
        # из взаимодействия с симметрической осью (стенкой клина).
        state = bc['steady_inflow_state']
        # Заполняем ghost-ячейки невозмущённым набегающим потоком
        Wstate = np.array([state['rho_0'], state['u_0'], state['v_0'], state['p_0']])
        Ustate = prim_to_cons(Wstate.reshape(4, 1, 1)).reshape(4)
        for k in range(g):
            for c in range(4):
                U[c, :, g - 1 - k] = Ustate[c]
    elif bc['left'] == 'steady_profile_inflow':
        # Профильный inflow для STEADY wind-tunnel постановки.
        # Линия УВ на левой границе: y_shock_left.
        # y > y_shock_left: свежий поток (state_0)
        # y < y_shock_left: за УВ (state_1)
        state = bc['left_state']
        y_shock = bc['left_y_shock']
        y_coords = bc['y_centers']  # массив y-координат ячеек (без ghost)
        ny_inner = len(y_coords)
        W0 = np.array([state['rho_0'], state['u_0'], state['v_0'], state['p_0']])
        W1 = np.array([state['rho_1'], state['u_1'], state['v_1'], state['p_1']])
        U0 = prim_to_cons(W0.reshape(4, 1, 1)).reshape(4)
        U1 = prim_to_cons(W1.reshape(4, 1, 1)).reshape(4)
        for j in range(ny_inner):
            j_full = j + g  # индекс в массиве с ghost
            Ufill = U1 if y_coords[j] < y_shock else U0
            for k in range(g):
                for c in range(4):
                    U[c, j_full, g - 1 - k] = Ufill[c]

    # Правая
    if bc['right'] == 'wall':
        nx = U.shape[2]
        for k in range(g):
            U[:, :, nx - g + k] = U[:, :, nx - g - 1 - k]
            U[1, :, nx - g + k] = -U[1, :, nx - g - 1 - k]
    elif bc['right'] == 'outflow':
        nx = U.shape[2]
        for k in range(g):
            U[:, :, nx - g + k] = U[:, :, nx - g - 1]

    # Нижняя
    if bc['bottom'] == 'wall':
        for k in range(g):
            U[:, g - 1 - k, :] = U[:, g + k, :]
            U[2, g - 1 - k, :] = -U[2, g + k, :]
    elif bc['bottom'] == 'outflow':
        for k in range(g):
            U[:, g - 1 - k, :] = U[:, g, :]

    # Верхняя
    if bc['top'] == 'wall':
        ny = U.shape[1]
        for k in range(g):
            U[:, ny - g + k, :] = U[:, ny - g - 1 - k, :]
            U[2, ny - g + k, :] = -U[2, ny - g - 1 - k, :]
    elif bc['top'] == 'outflow':
        ny = U.shape[1]
        for k in range(g):
            U[:, ny - g + k, :] = U[:, ny - g - 1, :]


# =============================================================================
# Один шаг по времени MUSCL-Hancock
# =============================================================================

def compute_dt(U: np.ndarray, dx: float, dy: float, cfl: float) -> float:
    """Шаг по времени по условию CFL."""
    W = cons_to_prim(U)
    rho, u, v, p = W[0], W[1], W[2], W[3]
    a = sound_speed(rho, p)
    smax_x = np.max(np.abs(u) + a)
    smax_y = np.max(np.abs(v) + a)
    dt = cfl / (smax_x / dx + smax_y / dy + 1e-30)
    return dt


def rhs_muscl(U: np.ndarray, dx: float, dy: float, bc: dict) -> np.ndarray:
    """
    Правая часть -dU/dt (т.е. величина, на которую U *уменьшается* за единицу времени)
    с реконструкцией MUSCL и HLLC. Без полушага по времени.
    
    Возвращает массив той же формы, что U; для ghost-ячеек значение нулевое
    (они всё равно перезаполнятся apply_bc).
    """
    apply_bc(U, bc)
    W = cons_to_prim(U)

    # x-направление
    WLx, WRx = reconstruct_x(W)  # формы (4, ny, nx-3)
    Fx = hllc_flux_x(WLx, WRx)   # (4, ny, nx-3) — потоки на гранях между i и i+1, для i = 1..nx-3

    # y-направление
    WLy, WRy = reconstruct_y(W)  # (4, ny-3, nx)
    Gy = hllc_flux_y(WLy, WRy)   # (4, ny-3, nx) — потоки на гранях между j и j+1, для j = 1..ny-3

    # Дивергенция: dU/dt = -(F_{i+1/2} - F_{i-1/2})/dx - (G_{j+1/2} - G_{j-1/2})/dy
    # Внутренние ячейки, для которых определены оба потока: 
    #   x: i от 2 до nx-3 (включительно) -> Fx[i-1] - Fx[i-2] для соседних граней
    #   y: j от 2 до ny-3
    rhs = np.zeros_like(U)
    rhs[:, 2:-2, 2:-2] = -(Fx[:, 2:-2, 1:] - Fx[:, 2:-2, :-1]) / dx \
                        -(Gy[:, 1:, 2:-2] - Gy[:, :-1, 2:-2]) / dy
    return rhs


def apply_wedge_bc(U: np.ndarray, wedge_mask: np.ndarray,
                    wedge_normal_x: np.ndarray, wedge_normal_y: np.ndarray,
                    rho_gas_default: float = 1.0, p_gas_default: float = 1.0,
                    gamma: float = GAMMA) -> None:
    """
    In-place применение BC для immersed wedge (твёрдый клин внутри расчётной области).
    
    Векторизованная реализация: для каждой клин-ячейки находим газового соседа 
    в направлении нормали (выходящей из клина), копируем его состояние с 
    инверсией нормальной компоненты скорости.
    """
    NY, NX = wedge_mask.shape
    if not np.any(wedge_mask):
        return
    
    # Сначала инициализируем все клин-ячейки стоячим газом (для изолированных)
    U[0][wedge_mask] = rho_gas_default
    U[1][wedge_mask] = 0.0
    U[2][wedge_mask] = 0.0
    U[3][wedge_mask] = p_gas_default / (gamma - 1.0)
    
    # Целочисленный сдвиг к ближайшему газовому соседу
    di = np.round(wedge_normal_x).astype(int)  # +1 (вправо), 0, или -1 (влево)
    dj = np.round(wedge_normal_y).astype(int)  # +1 (вверх), 0, или -1 (вниз)
    
    # Координаты клин-ячеек
    j_idx, i_idx = np.where(wedge_mask)
    if len(j_idx) == 0:
        return
    
    # Координаты соседей
    j_neigh = j_idx + dj[j_idx, i_idx]
    i_neigh = i_idx + di[j_idx, i_idx]
    
    # Фильтр валидности
    valid = ((j_neigh >= 0) & (j_neigh < NY)
             & (i_neigh >= 0) & (i_neigh < NX))
    valid &= ~wedge_mask[np.clip(j_neigh, 0, NY-1), np.clip(i_neigh, 0, NX-1)]
    
    j_w = j_idx[valid]
    i_w = i_idx[valid]
    j_n = j_neigh[valid]
    i_n = i_neigh[valid]
    
    if len(j_w) == 0:
        return
    
    # Копируем с отражением нормальной скорости
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


def step_ssprk2(U: np.ndarray, dt: float, dx: float, dy: float, bc: dict,
                 wedge_mask: Optional[np.ndarray] = None,
                 wedge_normal_x: Optional[np.ndarray] = None,
                 wedge_normal_y: Optional[np.ndarray] = None) -> np.ndarray:
    """SSP-RK2 (Heun) шаг по времени. Дает 2-й порядок и сохраняет монотонность с MUSCL.
    
    При наличии wedge_mask — после каждого подшага применяется immersed wedge BC.
    """
    L1 = rhs_muscl(U, dx, dy, bc)
    U1 = U + dt * L1
    if wedge_mask is not None:
        apply_wedge_bc(U1, wedge_mask, wedge_normal_x, wedge_normal_y)
    L2 = rhs_muscl(U1, dx, dy, bc)
    U_new = 0.5 * (U + U1 + dt * L2)
    if wedge_mask is not None:
        apply_wedge_bc(U_new, wedge_mask, wedge_normal_x, wedge_normal_y)
    return U_new


# =============================================================================
# Класс-обёртка
# =============================================================================

@dataclass
class EulerSim:
    """2D-симуляция уравнений Эйлера."""
    nx: int
    ny: int
    Lx: float
    Ly: float
    cfl: float = 0.45
    bc: dict = field(default_factory=dict)

    def __post_init__(self):
        # Полная сетка с ghost-ячейками
        self.dx = self.Lx / self.nx
        self.dy = self.Ly / self.ny
        self.NX = self.nx + 2 * NGHOST
        self.NY = self.ny + 2 * NGHOST
        # координаты центров физических ячеек
        self.x = (np.arange(self.nx) + 0.5) * self.dx
        self.y = (np.arange(self.ny) + 0.5) * self.dy
        # все ячейки включая ghost (для удобства задания НУ)
        self.x_full = (np.arange(self.NX) - NGHOST + 0.5) * self.dx
        self.y_full = (np.arange(self.NY) - NGHOST + 0.5) * self.dy

        self.U = np.zeros((4, self.NY, self.NX))
        self.t = 0.0
        
        # Если в bc используется steady_profile_inflow — передаём y-координаты
        if isinstance(self.bc, dict) and self.bc.get('left') == 'steady_profile_inflow':
            self.bc['y_centers'] = self.y
        
        # Wedge mask: по умолчанию None (нет встроенного клина)
        self.wedge_mask: Optional[np.ndarray] = None
        self.wedge_normal_x: Optional[np.ndarray] = None
        self.wedge_normal_y: Optional[np.ndarray] = None
    
    def set_wedge(self, x_corner: float, theta_w_deg: float,
                   location: str = 'top',
                   wedge_length: Optional[float] = None):
        """
        Установить immersed wedge.
        
        Аргументы:
            x_corner       : x-координата вершины клина
            theta_w_deg    : угол клина к набегающему потоку
            location       : 'top' или 'bottom' (где клин)
            wedge_length   : длина наклонной грани клина по оси x.
                             Если None — клин тянется до правой границы.
                             Иначе клин ограничен сверху x_corner + wedge_length.
        """
        theta_w = np.deg2rad(theta_w_deg)
        # Создаём 2D-сетку всех ячеек включая ghost
        X, Y = np.meshgrid(self.x_full, self.y_full)
        
        if wedge_length is not None:
            x_back = x_corner + wedge_length
        else:
            x_back = self.Lx + 100  # фактически бесконечность
        
        if location == 'top':
            # Поверхность клина: y_surface(x) = Ly - tan(θ_w) * (x - x_corner)
            # для x_corner ≤ x ≤ x_back.
            # При x > x_back — задняя грань клина (вертикальная) до Ly.
            y_surf = self.Ly - np.tan(theta_w) * (X - x_corner)
            mask = ((Y > y_surf) & (X > x_corner) & (X <= x_back))
            # Если wedge_length ограничен — добавляем заднюю вертикальную грань:
            # клин тянется по y от y_back = Ly - tan*(wedge_length) до Ly при x_back-1дх ≤ x ≤ x_back
            # На самом деле проще: клин = область с y > y_surf и x_corner < x < x_back
            # Сзади клин обрывается — это "задняя грань", тоже стенка.
            # Нормаль к поверхности клина (наружу в газ): (-sin θ, -cos θ)
            n_x = -np.sin(theta_w)
            n_y = -np.cos(theta_w)
            normal_x = np.full_like(X, n_x)
            normal_y = np.full_like(X, n_y)
        elif location == 'bottom':
            y_surf = np.tan(theta_w) * (X - x_corner)
            mask = ((Y < y_surf) & (X > x_corner) & (X <= x_back))
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
        """init_func(x, y) -> (rho, u, v, p) — возвращает 4-кортеж массивов или скаляров."""
        X, Y = np.meshgrid(self.x_full, self.y_full)
        W = np.zeros((4, self.NY, self.NX))
        rho, u, v, p = init_func(X, Y)
        W[0] = rho
        W[1] = u
        W[2] = v
        W[3] = p
        self.U = prim_to_cons(W)

    def get_primitive(self):
        """Возвращает примитивы только в физических ячейках (без ghost)."""
        g = NGHOST
        W = cons_to_prim(self.U[:, g:-g, g:-g])
        return W  # (4, ny, nx)

    def advance(self, t_end: float, verbose_every: Optional[int] = None,
                callback: Optional[Callable] = None):
        """Интегрирует до t_end."""
        n = 0
        while self.t < t_end:
            dt = compute_dt(self.U, self.dx, self.dy, self.cfl)
            if self.t + dt > t_end:
                dt = t_end - self.t
            self.U = step_ssprk2(self.U, dt, self.dx, self.dy, self.bc,
                                  wedge_mask=self.wedge_mask,
                                  wedge_normal_x=self.wedge_normal_x,
                                  wedge_normal_y=self.wedge_normal_y)
            self.t += dt
            n += 1
            if verbose_every is not None and n % verbose_every == 0:
                W = self.get_primitive()
                print(f"  step {n:5d}  t={self.t:.4f}  dt={dt:.2e}  "
                      f"rho∈[{W[0].min():.3f},{W[0].max():.3f}]")
            if callback is not None:
                callback(self, n)
        return n
    
    def advance_n_steps(self, n_steps: int,
                         callback: Optional[Callable] = None):
        """Интегрирует ровно n_steps шагов."""
        for n in range(n_steps):
            dt = compute_dt(self.U, self.dx, self.dy, self.cfl)
            self.U = step_ssprk2(self.U, dt, self.dx, self.dy, self.bc,
                                  wedge_mask=self.wedge_mask,
                                  wedge_normal_x=self.wedge_normal_x,
                                  wedge_normal_y=self.wedge_normal_y)
            self.t += dt
            if callback is not None:
                callback(self, n)
        return n_steps
    
    def update_bc(self, new_bc: dict):
        """Заменить граничные условия на лету (для quasi-static change клина).
        
        Сохраняет y_centers если требуется steady_profile_inflow.
        """
        if new_bc.get('left') == 'steady_profile_inflow':
            new_bc = dict(new_bc)
            new_bc['y_centers'] = self.y
        self.bc = new_bc