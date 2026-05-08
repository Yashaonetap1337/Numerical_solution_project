"""
Постановка задачи Иванов-Chpoun: ОДНОРОДНЫЙ сверхзвуковой поток + ФИЗИЧЕСКИЙ КЛИН.

В отличие от setup.py (где поле "за УВ" задаётся через профильный inflow), здесь:
    - Левая граница: однородный поток M=M_0, направление +x.
    - Нижняя стенка y=0: плоскость симметрии.
    - В верхней половине области: физический клин (immersed wedge),
      нижняя поверхность которого образует косую УВ при обтекании.
    - Угол клина θ_w задаётся через геометрию маски.

Геометрия клина (top-mounted wedge):
    Поверхность клина — прямая от точки (x_corner, Ly) под углом θ_w вниз-вправо.
    Уравнение: y = Ly - tan(θ_w) * (x - x_corner)  для x ≥ x_corner.
    Ячейки ВЫШЕ этой прямой и x ≥ x_corner — внутри клина (твёрдое тело).
    Ячейки НИЖЕ — газ.

Логика гистерезиса:
    Вращая клин (меняя θ_w), мы изменяем геометрию обтекания. Если поле уже было 
    в MR-конфигурации (большое θ_w), маховский стержень "помнит" свою структуру 
    при медленном уменьшении θ_w — это и есть классический Chpoun-Ivanov hysteresis.
"""
from __future__ import annotations
import numpy as np
from .solver import GAMMA, EulerSim


def make_uniform_inflow_bc(M: float, rho_0: float = 1.0, p_0: float = 1.0,
                            gamma: float = GAMMA) -> dict:
    """
    Граничные условия для wedge-in-flow задачи:
    - left:   однородный сверхзвуковой поток (M_0, +x)
    - right:  outflow
    - bottom: wall (плоскость симметрии)
    - top:    wall (внутри top-wedge BC применяется immersed BC; над клином 
              — обычная стенка, неважно что именно, всё равно перекрыто клином)
    """
    a_0 = np.sqrt(gamma * p_0 / rho_0)
    state = {
        'rho_0': rho_0, 'p_0': p_0,
        'u_0': M * a_0, 'v_0': 0.0, 'a_0': a_0,
        # Заглушки для полей за УВ — не используются в этой постановке,
        # но требуются apply_bc для ключа left='steady_profile_inflow'
        'rho_1': rho_0, 'p_1': p_0, 'u_1': M * a_0, 'v_1': 0.0,
    }
    return {
        'left': 'steady_profile_inflow',
        'left_state': state,
        'left_y_shock': -1.0,  # ниже всех y → весь профиль = свежий поток
        'left_beta_1': 0.0,
        'right': 'outflow',
        'bottom': 'wall',
        'top': 'wall',
    }


def make_uniform_initial_condition(M: float,
                                     rho_0: float = 1.0, p_0: float = 1.0,
                                     gamma: float = GAMMA):
    """Начальное условие: однородный поток везде."""
    a_0 = np.sqrt(gamma * p_0 / rho_0)
    u_0 = M * a_0
    
    def init(X, Y):
        rho = np.full_like(X, rho_0)
        u = np.full_like(X, u_0)
        v = np.zeros_like(X)
        p = np.full_like(X, p_0)
        return rho, u, v, p
    
    return init


def setup_wedge_simulation(M: float, theta_w_deg: float,
                            nx: int = 200, ny: int = 100,
                            Lx: float = 3.0, Ly: float = 1.0,
                            x_corner: float = 0.5,
                            wedge_length: float = None,
                            cfl: float = 0.4) -> EulerSim:
    """
    Полный setup wedge-in-flow задачи.
    
    Аргументы:
        wedge_length : длина наклонной грани клина по x. 
                       По умолчанию = (Lx - x_corner) / 2 (клин конечной длины).
                       Это даёт пространство ЗА клином для развития потока.
    """
    bc = make_uniform_inflow_bc(M)
    sim = EulerSim(nx=nx, ny=ny, Lx=Lx, Ly=Ly, cfl=cfl, bc=bc)
    
    init_func = make_uniform_initial_condition(M)
    sim.set_initial_primitive(init_func)
    
    if wedge_length is None:
        wedge_length = (Lx - x_corner) * 0.55
    
    sim.set_wedge(x_corner=x_corner, theta_w_deg=theta_w_deg, location='top',
                   wedge_length=wedge_length)
    
    return sim


def update_wedge_angle(sim: EulerSim, new_theta_w_deg: float):
    """Изменить угол клина на лету (для quasi-static rotation)."""
    sim.set_wedge(x_corner=sim.wedge_x_corner,
                   theta_w_deg=new_theta_w_deg,
                   location=sim.wedge_location,
                   wedge_length=sim.wedge_length)