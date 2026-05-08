"""
Постановка задачи Иванов-Chpoun: ОДНОРОДНЫЙ сверхзвуковой поток + ФИЗИЧЕСКИЙ КЛИН.

См. историю исправлений в комментариях к update_wedge_angle.
"""
from __future__ import annotations
import numpy as np
from .solver import GAMMA, EulerSim


def make_uniform_inflow_bc(M: float, rho_0: float = 1.0, p_0: float = 1.0,
                            gamma: float = GAMMA) -> dict:
    a_0 = np.sqrt(gamma * p_0 / rho_0)
    state = {
        'rho_0': rho_0, 'p_0': p_0,
        'u_0': M * a_0, 'v_0': 0.0, 'a_0': a_0,
        'rho_1': rho_0, 'p_1': p_0, 'u_1': M * a_0, 'v_1': 0.0,
    }
    return {
        'left': 'steady_profile_inflow',
        'left_state': state,
        'left_y_shock': -1.0,
        'left_beta_1': 0.0,
        'right': 'outflow',
        'bottom': 'wall',
        'top': 'wall',
    }


def make_uniform_initial_condition(M: float,
                                     rho_0: float = 1.0, p_0: float = 1.0,
                                     gamma: float = GAMMA):
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
    bc = make_uniform_inflow_bc(M)
    sim = EulerSim(nx=nx, ny=ny, Lx=Lx, Ly=Ly, cfl=cfl, bc=bc)

    init_func = make_uniform_initial_condition(M)
    sim.set_initial_primitive(init_func)

    if wedge_length is None:
        wedge_length = (Lx - x_corner) * 0.55

    sim.set_wedge(x_corner=x_corner, theta_w_deg=theta_w_deg, location='top',
                   wedge_length=wedge_length)

    return sim


def update_wedge_angle(sim: EulerSim, new_theta_w_deg: float,
                        max_reinit_iters: int = 50,
                        verbose: bool = False):
    """
    Изменить угол клина на лету с ИТЕРАТИВНОЙ реинициализацией ячеек,
    переходящих из клина обратно в газ.

    Алгоритм flood-fill: на каждой итерации заполняются released-ячейки,
    у которых уже есть хотя бы один валидный сосед (8-связно). После
    заполнения они сами становятся валидными для следующей итерации.
    Так волна валидных значений распространяется внутрь стрипа любой
    толщины.
    """
    if sim.wedge_mask is None:
        sim.set_wedge(x_corner=sim.wedge_x_corner,
                      theta_w_deg=new_theta_w_deg,
                      location=sim.wedge_location,
                      wedge_length=sim.wedge_length)
        return

    old_mask = sim.wedge_mask.copy()

    sim.set_wedge(x_corner=sim.wedge_x_corner,
                   theta_w_deg=new_theta_w_deg,
                   location=sim.wedge_location,
                   wedge_length=sim.wedge_length)

    new_mask = sim.wedge_mask
    released = old_mask & ~new_mask

    if not released.any():
        return

    U = sim.U
    NY, NX = old_mask.shape

    # Стартово валидны те ячейки, которые были и остались fluid
    valid = (~old_mask) & (~new_mask)
    pending = released.copy()

    n_initial_pending = int(pending.sum())
    iters_used = 0

    for it in range(max_reinit_iters):
        if not pending.any():
            break
        iters_used = it + 1

        filled_this_iter = np.zeros_like(pending)
        j_idx, i_idx = np.where(pending)

        for j, i in zip(j_idx, i_idx):
            accum = np.zeros(4)
            count = 0
            for dj in (-1, 0, 1):
                for di in (-1, 0, 1):
                    if dj == 0 and di == 0:
                        continue
                    jn, in_ = j + dj, i + di
                    if 0 <= jn < NY and 0 <= in_ < NX:
                        if valid[jn, in_]:
                            accum += U[:, jn, in_]
                            count += 1
            if count > 0:
                U[:, j, i] = accum / count
                filled_this_iter[j, i] = True

        if not filled_this_iter.any():
            break

        valid |= filled_this_iter
        pending &= ~filled_this_iter

    n_orphan = int(pending.sum())

    if verbose:
        print(f"  [update_wedge_angle] released={n_initial_pending}, "
              f"reinit за {iters_used} итераций, orphan={n_orphan}")

    if n_orphan > 0:
        j_o, i_o = np.where(pending)
        fluid_cells = (~new_mask)
        if fluid_cells.any():
            fallback = np.array([
                np.median(U[k][fluid_cells]) for k in range(4)
            ])
            for j, i in zip(j_o, i_o):
                U[:, j, i] = fallback
        if n_orphan > 5:
            print(f"  [update_wedge_angle] WARNING: {n_orphan} orphan cells "
                  f"после {iters_used} итераций — использован median fallback.")