"""
Аналитические критерии перехода RR↔MR для STEADY wind-tunnel геометрии.

Геометрия (Ben-Dor 2007 §2.4, стр. 75-83):
- Однородный сверхзвуковой поток с числом Маха M_0 = M.
- Клин с углом θ_w к набегающему потоку.
- Поток отклоняется на угол θ_1 = θ_w через падающую косую УВ.
- Угол падающей УВ к потоку: β_1.

Все углы хранятся и возвращаются в ГРАДУСАХ (исключая внутренние вычисления).

КРИТЕРИИ ПЕРЕХОДА (steady, Ben-Dor §1.5 + §2.4):

1. Detachment (θ_w = θ_w^D):
   На R-полар максимальное отклонение θ_2_max(M_2) равно θ_1.
   Условие: θ_w = θ_max(M_2(M, β_1(θ_w))).
   Это верхняя граница: при θ_w > θ_w^D возможно только MR.

2. Sonic (θ_w = θ_w^s):
   Поток за отражённой УВ (M_3 относительно точки отражения) равен 1.
   Численно очень близок к detachment (разница ~0.05° для умеренных M).

3. Mechanical-equilibrium / von Neumann (θ_w = θ_w^N):
   R-полар пересекает p-ось ИМЕННО в точке нормальной УВ I-полар.
   Эквивалентно: давление за отражённой УВ (weak branch) при отклонении 
   θ_R = θ_w равно давлению за нормальной УВ при числе Маха M.
   Существует только при M ≥ M_0c ≈ 2.202.

ВЕРИФИКАЦИЯ (см. tests/test_theory.py):
- M=5.0:   θ_w^D = 27.77°,  θ_w^N = 20.86°
- M=4.96:  β_1^D = 39.33°,  β_1^N = 30.88°  (Chpoun et al. 1995)
- M<2.202: von_neumann_wedge_angle = None (vN не существует)
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from dataclasses import dataclass
from typing import Optional

GAMMA = 1.4
M_0_CRITICAL = 2.202  # ниже этого числа Маха vN-критерий не существует (γ=1.4)


# =============================================================================
# Базовые соотношения косых УВ
# =============================================================================

def oblique_deflection(M: float, beta: float, gamma: float = GAMMA) -> float:
    """Угол отклонения потока через косую УВ. β,θ — в радианах. Возвращает θ."""
    Mn = M * np.sin(beta)
    if Mn <= 1.0:
        return 0.0
    num = 2.0 * (Mn**2 - 1.0) / np.tan(beta)
    den = M**2 * (gamma + np.cos(2.0 * beta)) + 2.0
    return np.arctan(num / den)


def oblique_M2(M: float, beta: float, gamma: float = GAMMA) -> float:
    """Число Маха за косой УВ."""
    Mn = M * np.sin(beta)
    if Mn <= 1.0:
        return M
    Mn2_sq = (1.0 + (gamma - 1.0)/2.0 * Mn**2) / (gamma * Mn**2 - (gamma - 1.0)/2.0)
    Mn2 = np.sqrt(Mn2_sq)
    theta = oblique_deflection(M, beta, gamma)
    return Mn2 / np.sin(beta - theta)


def oblique_pressure_ratio(M: float, beta: float, gamma: float = GAMMA) -> float:
    """Отношение давлений за/перед косой УВ."""
    Mn = M * np.sin(beta)
    if Mn <= 1.0:
        return 1.0
    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (Mn**2 - 1.0)


def normal_shock_pressure_ratio(M: float, gamma: float = GAMMA) -> float:
    """Отношение давлений за/перед нормальной УВ."""
    if M <= 1.0:
        return 1.0
    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (M**2 - 1.0)


def theta_max_for_M(M: float, gamma: float = GAMMA) -> float:
    """Максимальное отклонение потока через косую УВ при заданном M (рад.)."""
    if M <= 1.0:
        return 0.0
    res = minimize_scalar(lambda b: -oblique_deflection(M, b, gamma),
                          bounds=(np.arcsin(1.0/M) + 1e-7, np.pi/2 - 1e-7),
                          method='bounded',
                          options={'xatol': 1e-10})
    return -res.fun


def beta_at_max_deflection(M: float, gamma: float = GAMMA) -> float:
    """Угол косой УВ при максимальном отклонении (рад.)."""
    if M <= 1.0:
        return 0.0
    res = minimize_scalar(lambda b: -oblique_deflection(M, b, gamma),
                          bounds=(np.arcsin(1.0/M) + 1e-7, np.pi/2 - 1e-7),
                          method='bounded',
                          options={'xatol': 1e-10})
    return res.x


def beta_for_deflection(M: float, theta: float, branch: str = 'weak',
                         gamma: float = GAMMA) -> Optional[float]:
    """Найти угол β косой УВ при отклонении θ. branch ∈ {'weak', 'strong'}."""
    th_max = theta_max_for_M(M, gamma)
    if theta > th_max + 1e-9:
        return None
    if abs(theta - th_max) < 1e-9:
        return beta_at_max_deflection(M, gamma)
    beta_at_max = beta_at_max_deflection(M, gamma)
    if branch == 'weak':
        try:
            return brentq(lambda b: oblique_deflection(M, b, gamma) - theta,
                          np.arcsin(1.0/M) + 1e-7, beta_at_max, xtol=1e-10)
        except ValueError:
            return None
    elif branch == 'strong':
        try:
            return brentq(lambda b: oblique_deflection(M, b, gamma) - theta,
                          beta_at_max, np.pi/2 - 1e-7, xtol=1e-10)
        except ValueError:
            return None
    else:
        raise ValueError(f"branch must be 'weak' or 'strong', got {branch!r}")


# =============================================================================
# Критерии перехода RR↔MR
# =============================================================================

def detachment_wedge_angle(M: float, gamma: float = GAMMA) -> Optional[float]:
    """
    Detachment criterion (steady, Ben-Dor §1.5.1, §2.4):
    Угол клина θ_w^D, при котором достигается максимум отклонения через 
    отражённую УВ: θ_2_max(M_2) = θ_1 = θ_w.

    Возвращает θ_w^D в ГРАДУСАХ. None, если решения нет.

    Метод: brentq на условии θ_w − θ_max(M_2(M, β_1(θ_w))) = 0.
    """
    if M <= 1.0:
        return None
    
    def residual(theta_w_deg: float) -> float:
        theta_w = np.deg2rad(theta_w_deg)
        beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
        if beta_1 is None:
            return 1e10
        M_2 = oblique_M2(M, beta_1, gamma)
        if M_2 <= 1.0:
            return 1e10
        return theta_w - theta_max_for_M(M_2, gamma)

    # Интервал поиска: θ_w = 0..θ_max(M)
    th_max_M = np.rad2deg(theta_max_for_M(M, gamma))
    if th_max_M < 1.0:
        return None
    try:
        return brentq(residual, 0.5, th_max_M - 0.01, xtol=1e-7)
    except ValueError:
        return None


def sonic_wedge_angle(M: float, gamma: float = GAMMA) -> Optional[float]:
    """
    Sonic criterion (Ben-Dor §1.5.3):
    Угол клина θ_w^s, при котором поток за отражённой УВ (M_3 относительно 
    точки отражения R) становится звуковым: M_3 = 1.

    Возвращает θ_w^s в ГРАДУСАХ. None, если решения нет.
    """
    if M <= 1.0:
        return None

    def residual(theta_w_deg: float) -> float:
        theta_w = np.deg2rad(theta_w_deg)
        beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
        if beta_1 is None:
            return -1e10
        M_2 = oblique_M2(M, beta_1, gamma)
        if M_2 <= 1.0:
            return -1e10
        # Отражённая УВ (weak branch): возвращает поток на угол θ_w
        beta_R = beta_for_deflection(M_2, theta_w, 'weak', gamma)
        if beta_R is None:
            return -1e10
        M_3 = oblique_M2(M_2, beta_R, gamma)
        return M_3 - 1.0

    th_max_M = np.rad2deg(theta_max_for_M(M, gamma))
    if th_max_M < 1.0:
        return None
    try:
        return brentq(residual, 0.5, th_max_M - 0.01, xtol=1e-7)
    except ValueError:
        return None


def von_neumann_wedge_angle(M: float, gamma: float = GAMMA) -> Optional[float]:
    """
    Mechanical-equilibrium / von Neumann criterion (Ben-Dor §1.5.2, ур. 1.32):
    R-полар пересекает p-ось ИМЕННО в точке нормальной УВ I-полар.
    Численно: давление за отражённой УВ (weak branch) при отклонении θ_R = θ_w 
    равно давлению за нормальной УВ при числе Маха M_0 = M.

    Возвращает θ_w^N в ГРАДУСАХ.
    None, если M < M_0_CRITICAL = 2.202 (vN не существует).
    
    Около критической точки M ≈ M_0c кривая vN сходится к detachment-кривой
    (это точка K на фазовой диаграмме). В этой области residual может НЕ 
    менять знак (только касается нуля), поэтому используется fallback на 
    минимум |residual|.
    """
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

    # Грубый скан
    th_max_M = np.rad2deg(theta_max_for_M(M, gamma))
    phi_arr = np.linspace(0.5, th_max_M - 0.001, 500)
    valid_phi, valid_res = [], []
    for p in phi_arr:
        r = residual(p)
        if r is not None and abs(r) < 1e9:
            valid_phi.append(p)
            valid_res.append(r)
    if len(valid_phi) < 2:
        return None
    valid_phi = np.array(valid_phi)
    valid_res = np.array(valid_res)
    
    # Поиск смены знака
    sign_changes = np.where(np.diff(np.sign(valid_res)) != 0)[0]
    if len(sign_changes) > 0:
        i = sign_changes[0]
        try:
            return brentq(lambda x: residual(x) if residual(x) is not None else 1e10,
                          valid_phi[i], valid_phi[i+1], xtol=1e-7)
        except (ValueError, RuntimeError):
            pass
    
    # Fallback: residual не пересекает 0, но может касаться (близко к точке K).
    # При M ≈ M_0c кривая vN сходится с detachment-кривой, и residual только 
    # слегка касается 0. Используем точку min|residual| как приближение vN.
    abs_res = np.abs(valid_res)
    i_min = np.argmin(abs_res)
    if abs_res[i_min] < 5.0:
        return float(valid_phi[i_min])
    
    return None


# =============================================================================
# Удобные API: углы падения β_1 на критических кривых
# =============================================================================

def beta_1_at_detachment(M: float, gamma: float = GAMMA) -> Optional[float]:
    """Угол падающей УВ β_1 на границе detachment (в градусах)."""
    tw = detachment_wedge_angle(M, gamma)
    if tw is None:
        return None
    beta_1 = beta_for_deflection(M, np.deg2rad(tw), 'weak', gamma)
    return None if beta_1 is None else float(np.rad2deg(beta_1))


def beta_1_at_von_neumann(M: float, gamma: float = GAMMA) -> Optional[float]:
    """Угол падающей УВ β_1 на границе von Neumann (в градусах)."""
    tw = von_neumann_wedge_angle(M, gamma)
    if tw is None:
        return None
    beta_1 = beta_for_deflection(M, np.deg2rad(tw), 'weak', gamma)
    return None if beta_1 is None else float(np.rad2deg(beta_1))


def beta_1_at_sonic(M: float, gamma: float = GAMMA) -> Optional[float]:
    """Угол падающей УВ β_1 на границе sonic (в градусах)."""
    tw = sonic_wedge_angle(M, gamma)
    if tw is None:
        return None
    beta_1 = beta_for_deflection(M, np.deg2rad(tw), 'weak', gamma)
    return None if beta_1 is None else float(np.rad2deg(beta_1))


# =============================================================================
# Вспомогательное состояние (для дашборда / визуализации)
# =============================================================================

@dataclass
class SteadyShockState:
    """Состояние steady-отражения при заданных (M, θ_w)."""
    M_0: float          # число Маха набегающего потока
    theta_w: float      # угол клина (град)
    beta_1: float       # угол падающей УВ (град)
    M_1: float          # число Маха за падающей УВ
    p_1_p_0: float      # отношение давлений за падающей УВ


def steady_state_at(M: float, theta_w_deg: float,
                     gamma: float = GAMMA) -> Optional[SteadyShockState]:
    """Вычислить состояние при steady-отражении (M, θ_w)."""
    theta_w = np.deg2rad(theta_w_deg)
    beta_1 = beta_for_deflection(M, theta_w, 'weak', gamma)
    if beta_1 is None:
        return None
    M_1 = oblique_M2(M, beta_1, gamma)
    p_1_p_0 = oblique_pressure_ratio(M, beta_1, gamma)
    return SteadyShockState(
        M_0=M,
        theta_w=theta_w_deg,
        beta_1=float(np.rad2deg(beta_1)),
        M_1=M_1,
        p_1_p_0=p_1_p_0,
    )


# =============================================================================
# Phase diagram curves
# =============================================================================

def phase_diagram_curves(M_range=(2.0, 7.0), npts=80, gamma=GAMMA):
    """Возвращает (M_arr, theta_d_arr, theta_vN_arr) — теоретические кривые 
    detachment и von Neumann для построения фазовой диаграммы."""
    M_arr = np.linspace(M_range[0], M_range[1], npts)
    theta_d, theta_vN = [], []
    for M in M_arr:
        td = detachment_wedge_angle(M, gamma)
        tvN = von_neumann_wedge_angle(M, gamma)
        theta_d.append(td if td is not None else np.nan)
        theta_vN.append(tvN if tvN is not None else np.nan)
    return M_arr, np.array(theta_d), np.array(theta_vN)


# =============================================================================
# CLI
# =============================================================================

if __name__ == '__main__':
    print(f"{'M':>5} {'θ_w^D':>10} {'β_1^D':>10} {'θ_w^s':>10} {'θ_w^N':>10} {'β_1^N':>10}")
    print("-" * 70)
    for M in [2.0, 2.202, 2.5, 3.0, 4.0, 4.96, 5.0, 6.0, 7.0]:
        td = detachment_wedge_angle(M)
        ts = sonic_wedge_angle(M)
        tvN = von_neumann_wedge_angle(M)
        b1d = beta_1_at_detachment(M)
        b1vN = beta_1_at_von_neumann(M)
        td_s = f"{td:.3f}°" if td else "  ---"
        ts_s = f"{ts:.3f}°" if ts else "  ---"
        tvN_s = f"{tvN:.3f}°" if tvN else "  ---"
        b1d_s = f"{b1d:.3f}°" if b1d else "  ---"
        b1vN_s = f"{b1vN:.3f}°" if b1vN else "  ---"
        print(f"{M:>5.2f} {td_s:>10} {b1d_s:>10} {ts_s:>10} {tvN_s:>10} {b1vN_s:>10}")