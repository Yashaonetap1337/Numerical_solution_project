from __future__ import annotations
import numpy as np
from .solver import GAMMA, EulerSim

def make_uniform_inflow_bc(M: float, rho_0: float=1.0, p_0: float=1.0, gamma: float=GAMMA) -> dict:
    a_0 = np.sqrt(gamma * p_0 / rho_0)
    state = {'rho_0': rho_0, 'p_0': p_0, 'u_0': M * a_0, 'v_0': 0.0, 'a_0': a_0, 'rho_1': rho_0, 'p_1': p_0, 'u_1': M * a_0, 'v_1': 0.0}
    return {'left': 'steady_profile_inflow', 'left_state': state, 'left_y_shock': -1.0, 'left_beta_1': 0.0, 'right': 'outflow', 'bottom': 'wall', 'top': 'wall'}

def make_uniform_initial_condition(M: float, rho_0: float=1.0, p_0: float=1.0, gamma: float=GAMMA):
    a_0 = np.sqrt(gamma * p_0 / rho_0)
    u_0 = M * a_0

    def init(X, Y):
        rho = np.full_like(X, rho_0)
        u = np.full_like(X, u_0)
        v = np.zeros_like(X)
        p = np.full_like(X, p_0)
        return (rho, u, v, p)
    return init

def setup_wedge_simulation(M: float, theta_w_deg: float, nx: int=200, ny: int=100, Lx: float=3.0, Ly: float=1.0, x_corner: float=0.5, wedge_length: float=None, cfl: float=0.4) -> EulerSim:
    bc = make_uniform_inflow_bc(M)
    sim = EulerSim(nx=nx, ny=ny, Lx=Lx, Ly=Ly, cfl=cfl, bc=bc)
    init_func = make_uniform_initial_condition(M)
    sim.set_initial_primitive(init_func)
    if wedge_length is None:
        wedge_length = (Lx - x_corner) * 0.55
    sim.set_wedge(x_corner=x_corner, theta_w_deg=theta_w_deg, location='top', wedge_length=wedge_length)
    return sim

def update_wedge_angle(sim: EulerSim, new_theta_w_deg: float):
    sim.set_wedge(x_corner=sim.wedge_x_corner, theta_w_deg=new_theta_w_deg, location=sim.wedge_location, wedge_length=sim.wedge_length)
