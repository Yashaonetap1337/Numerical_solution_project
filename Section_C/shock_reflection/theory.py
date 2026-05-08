from __future__ import annotations
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from dataclasses import dataclass
from typing import Optional
GAMMA = 1.4
M_0_CRITICAL = 2.202

def oblique_deflection(M: float, beta: float, gamma: float=GAMMA) -> float:
    Mn = M * np.sin(beta)
    if Mn <= 1.0:
        return 0.0
    num = 2.0 * (Mn ** 2 - 1.0) / np.tan(beta)
    den = M ** 2 * (gamma + np.cos(2.0 * beta)) + 2.0
    return np.arctan(num / den)

def oblique_M2(M: float, beta: float, gamma: float=GAMMA) -> float:
    Mn = M * np.sin(beta)
    if Mn <= 1.0:
        return M
    Mn2_sq = (1.0 + (gamma - 1.0) / 2.0 * Mn ** 2) / (gamma * Mn ** 2 - (gamma - 1.0) / 2.0)
    Mn2 = np.sqrt(Mn2_sq)
    theta = oblique_deflection(M, beta, gamma)
    return Mn2 / np.sin(beta - theta)

def oblique_pressure_ratio(M: float, beta: float, gamma: float=GAMMA) -> float:
    Mn = M * np.sin(beta)
    if Mn <= 1.0:
        return 1.0
    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (Mn ** 2 - 1.0)

def normal_shock_pressure_ratio(M: float, gamma: float=GAMMA) -> float:
    if M <= 1.0:
        return 1.0
    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (M ** 2 - 1.0)

def theta_max_for_M(M: float, gamma: float=GAMMA) -> float:
    if M <= 1.0:
        return 0.0
    res = minimize_scalar(lambda b: -oblique_deflection(M, b, gamma), bounds=(np.arcsin(1.0 / M) + 1e-07, np.pi / 2 - 1e-07), method='bounded', options={'xatol': 1e-10})
    return -res.fun

def beta_at_max_deflection(M: float, gamma: float=GAMMA) -> float:
    if M <= 1.0:
        return 0.0
    res = minimize_scalar(lambda b: -oblique_deflection(M, b, gamma), bounds=(np.arcsin(1.0 / M) + 1e-07, np.pi / 2 - 1e-07), method='bounded', options={'xatol': 1e-10})
    return res.x

def beta_for_deflection(M: float, theta: float, branch: str='weak', gamma: float=GAMMA) -> Optional[float]:
    th_max = theta_max_for_M(M, gamma)
    if theta > th_max + 1e-09:
        return None
    if abs(theta - th_max) < 1e-09:
        return beta_at_max_deflection(M, gamma)
    beta_at_max = beta_at_max_deflection(M, gamma)
    if branch == 'weak':
        try:
            return brentq(lambda b: oblique_deflection(M, b, gamma) - theta, np.arcsin(1.0 / M) + 1e-07, beta_at_max, xtol=1e-10)
        except ValueError:
            return None
    elif branch == 'strong':
        try:
            return brentq(lambda b: oblique_deflection(M, b, gamma) - theta, beta_at_max, np.pi / 2 - 1e-07, xtol=1e-10)
        except ValueError:
            return None
    else:
        raise ValueError(f"branch must be 'weak' or 'strong', got {branch!r}")

def detachment_wedge_angle(M: float, gamma: float=GAMMA) -> Optional[float]:
    if M <= 1.0:
        return None

    def residual(theta_w_deg: float) -> float:
        theta_w = np.deg2rad(theta_w_deg)
        beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
        if beta_1 is None:
            return 10000000000.0
        M_2 = oblique_M2(M, beta_1, gamma)
        if M_2 <= 1.0:
            return 10000000000.0
        return theta_w - theta_max_for_M(M_2, gamma)
    th_max_M = np.rad2deg(theta_max_for_M(M, gamma))
    if th_max_M < 1.0:
        return None
    try:
        return brentq(residual, 0.5, th_max_M - 0.01, xtol=1e-07)
    except ValueError:
        return None

def sonic_wedge_angle(M: float, gamma: float=GAMMA) -> Optional[float]:
    if M <= 1.0:
        return None

    def residual(theta_w_deg: float) -> float:
        theta_w = np.deg2rad(theta_w_deg)
        beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
        if beta_1 is None:
            return -10000000000.0
        M_2 = oblique_M2(M, beta_1, gamma)
        if M_2 <= 1.0:
            return -10000000000.0
        beta_R = beta_for_deflection(M_2, theta_w, 'weak', gamma)
        if beta_R is None:
            return -10000000000.0
        M_3 = oblique_M2(M_2, beta_R, gamma)
        return M_3 - 1.0
    th_max_M = np.rad2deg(theta_max_for_M(M, gamma))
    if th_max_M < 1.0:
        return None
    try:
        return brentq(residual, 0.5, th_max_M - 0.01, xtol=1e-07)
    except ValueError:
        return None

def von_neumann_wedge_angle(M: float, gamma: float=GAMMA) -> Optional[float]:
    if M < M_0_CRITICAL:
        return None
    p_normal_M = normal_shock_pressure_ratio(M, gamma)

    def residual(theta_w_deg: float) -> float:
        theta_w = np.deg2rad(theta_w_deg)
        beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
        if beta_1 is None:
            return None
        p_1_p_0 = oblique_pressure_ratio(M, beta_1, gamma)
        M_2 = oblique_M2(M, beta_1, gamma)
        if M_2 <= 1.0:
            return None
        beta_R = beta_for_deflection(M_2, theta_w, 'weak', gamma)
        if beta_R is None:
            return None
        p_2_p_1 = oblique_pressure_ratio(M_2, beta_R, gamma)
        p_2_p_0 = p_2_p_1 * p_1_p_0
        return p_2_p_0 - p_normal_M
    th_max_M = np.rad2deg(theta_max_for_M(M, gamma))
    phi_arr = np.linspace(0.5, th_max_M - 0.001, 500)
    valid_phi, valid_res = ([], [])
    for p in phi_arr:
        r = residual(p)
        if r is not None and abs(r) < 1000000000.0:
            valid_phi.append(p)
            valid_res.append(r)
    if len(valid_phi) < 2:
        return None
    valid_phi = np.array(valid_phi)
    valid_res = np.array(valid_res)
    sign_changes = np.where(np.diff(np.sign(valid_res)) != 0)[0]
    if len(sign_changes) > 0:
        i = sign_changes[0]
        try:
            return brentq(lambda x: residual(x) if residual(x) is not None else 10000000000.0, valid_phi[i], valid_phi[i + 1], xtol=1e-07)
        except (ValueError, RuntimeError):
            pass
    abs_res = np.abs(valid_res)
    i_min = np.argmin(abs_res)
    if abs_res[i_min] < 5.0:
        return float(valid_phi[i_min])
    return None

def beta_1_at_detachment(M: float, gamma: float=GAMMA) -> Optional[float]:
    tw = detachment_wedge_angle(M, gamma)
    if tw is None:
        return None
    beta_1 = beta_for_deflection(M, np.deg2rad(tw), 'weak', gamma)
    return None if beta_1 is None else float(np.rad2deg(beta_1))

def beta_1_at_von_neumann(M: float, gamma: float=GAMMA) -> Optional[float]:
    tw = von_neumann_wedge_angle(M, gamma)
    if tw is None:
        return None
    beta_1 = beta_for_deflection(M, np.deg2rad(tw), 'weak', gamma)
    return None if beta_1 is None else float(np.rad2deg(beta_1))

def beta_1_at_sonic(M: float, gamma: float=GAMMA) -> Optional[float]:
    tw = sonic_wedge_angle(M, gamma)
    if tw is None:
        return None
    beta_1 = beta_for_deflection(M, np.deg2rad(tw), 'weak', gamma)
    return None if beta_1 is None else float(np.rad2deg(beta_1))

@dataclass
class SteadyShockState:
    M_0: float
    theta_w: float
    beta_1: float
    M_1: float
    p_1_p_0: float

def steady_state_at(M: float, theta_w_deg: float, gamma: float=GAMMA) -> Optional[SteadyShockState]:
    theta_w = np.deg2rad(theta_w_deg)
    beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
    if beta_1 is None:
        return None
    M_1 = oblique_M2(M, beta_1, gamma)
    p_1_p_0 = oblique_pressure_ratio(M, beta_1, gamma)
    return SteadyShockState(M_0=M, theta_w=theta_w_deg, beta_1=float(np.rad2deg(beta_1)), M_1=M_1, p_1_p_0=p_1_p_0)

def phase_diagram_curves(M_range=(2.0, 7.0), npts=80, gamma=GAMMA):
    M_arr = np.linspace(M_range[0], M_range[1], npts)
    theta_d, theta_vN = ([], [])
    for M in M_arr:
        td = detachment_wedge_angle(M, gamma)
        tvN = von_neumann_wedge_angle(M, gamma)
        theta_d.append(td if td is not None else np.nan)
        theta_vN.append(tvN if tvN is not None else np.nan)
    return (M_arr, np.array(theta_d), np.array(theta_vN))
if __name__ == '__main__':
    print(f"{'M':>5} {'θ_w^D':>10} {'β_1^D':>10} {'θ_w^s':>10} {'θ_w^N':>10} {'β_1^N':>10}")
    print('-' * 70)
    for M in [2.0, 2.202, 2.5, 3.0, 4.0, 4.96, 5.0, 6.0, 7.0]:
        td = detachment_wedge_angle(M)
        ts = sonic_wedge_angle(M)
        tvN = von_neumann_wedge_angle(M)
        b1d = beta_1_at_detachment(M)
        b1vN = beta_1_at_von_neumann(M)
        td_s = f'{td:.3f}°' if td else '  ---'
        ts_s = f'{ts:.3f}°' if ts else '  ---'
        tvN_s = f'{tvN:.3f}°' if tvN else '  ---'
        b1d_s = f'{b1d:.3f}°' if b1d else '  ---'
        b1vN_s = f'{b1vN:.3f}°' if b1vN else '  ---'
        print(f'{M:>5.2f} {td_s:>10} {b1d_s:>10} {ts_s:>10} {tvN_s:>10} {b1vN_s:>10}')
