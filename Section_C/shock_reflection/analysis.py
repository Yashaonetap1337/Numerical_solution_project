"""
Детектор RR/MR и тройной точки для wedge-in-flow геометрии.

ФИЗИЧЕСКАЯ ОСНОВА:
В MR-конфигурации в тройной точке T сходятся падающая УВ I, отражённая R,
маховский стержень S и контактный разрыв (slip line, vortex sheet) V.

Контактный разрыв V:
    - плотность ρ скачет (поток разделён на «за RR-системой» и «за стержнем»);
    - давление p непрерывно (slip-line изобарическая);
    - касательная скорость может скачком меняться.

В RR контактного разрыва нет.

ИДЕЯ ДЕТЕКТОРА:
1. Считаются |∇ρ| и |∇p|, нормируются на свои максимумы.
2. Индикатор slip_score = |∇ρ|_n · (1 − |∇p|_n) ∈ [0,1] максимален именно
   на контактных разрывах (на УВ |∇p| велик, и (1-|∇p|_n) подавляет score).
3. Из области поиска исключаются:
    - расширенный клин (избегаем ложного сигнала от энтропийного слоя
      у поверхности клина — численный артефакт wedge BC);
    - верхняя часть домена y > y_search_max_fraction · Ly (отсекает
      spurious-сигналы из области wake'а за концом клина и т.п.);
    - 3 ячейки на левой и правой границах.
4. Тройная точка = argmax slip_score в оставшейся области.
   Если max slip < slip_threshold → RR (нет slip-line).

КАЛИБРОВКА (на 4 ground-truth кейсах):
    M=3, θ=25°: TRUE (5.675, 0.025), детектор (5.58, 0.07), ошибка ~0.1
    M=3, θ=28°: TRUE (4.575, 0.375), детектор (4.88, 0.17), ошибка ~0.3
    M=3, θ=30°: TRUE (4.125, 0.425), детектор (4.28, 0.33), ошибка ~0.15
    M=5, θ=35°: TRUE (4.075, 0.675), детектор (4.17, 0.68), ошибка ~0.1

Все ошибки в пределах ~6 ячеек на сетке dy=0.05 — сравнимо с толщиной
размытого УВ после численной диффузии.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional
from scipy.ndimage import binary_dilation


@dataclass
class ReflectionDiagnosis:
    type: str               # 'RR' или 'MR'
    x_TP: Optional[float]
    y_TP: Optional[float]
    h_stem: float           # = y_TP, длина маховского стержня
    h_relative: float       # h_stem / Ly
    confidence: float = 1.0


def detect_reflection(rho: np.ndarray, x: np.ndarray, y: np.ndarray,
                      p: Optional[np.ndarray] = None,
                      wedge_mask: Optional[np.ndarray] = None,
                      slip_threshold: float = 0.20,
                      wedge_dilation_iters: int = 4,
                      y_search_max_fraction: float = 0.18,
                      ) -> ReflectionDiagnosis:
    """
    Детектор RR/MR через признак slip-line. Тройная точка = argmax
    slip_score в области поиска (исключая клин + верх + границы).

    Параметры:
        rho, p              : поля плотности и давления (без ghost).
        x, y                : 1D координаты центров ячеек.
        wedge_mask          : маска клина (True где клин). None → клин
                              не исключается из поиска.
        slip_threshold      : минимальный slip_score для принятия точки
                              как тройной (иначе RR). 0.20 — стандарт.
        wedge_dilation_iters: насколько ячеек расширить wedge_mask при
                              исключении из поиска. 4 ячейки = 0.20 ед.
                              на стандартной сетке — закрывает энтропийный
                              слой у поверхности клина.
        y_search_max_fraction: верхняя граница поиска по y, в долях Ly.
                              0.18 (= y < 0.9 при Ly=5) отсекает
                              spurious-сигналы в средней части домена.
    """
    if p is None:
        p = rho ** 1.4   # для unit-тестов

    ny, nx = rho.shape
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    Ly = y[-1] + dy

    drho_dx = np.gradient(rho, dx, axis=1)
    drho_dy = np.gradient(rho, dy, axis=0)
    G_rho = np.sqrt(drho_dx**2 + drho_dy**2)

    dp_dx = np.gradient(p, dx, axis=1)
    dp_dy = np.gradient(p, dy, axis=0)
    G_p = np.sqrt(dp_dx**2 + dp_dy**2)

    g_rho_norm = G_rho / (G_rho.max() + 1e-12)
    g_p_norm   = G_p   / (G_p.max()   + 1e-12)

    slip_score = g_rho_norm * (1.0 - g_p_norm)

    exclude = np.zeros_like(slip_score, dtype=bool)
    if wedge_mask is not None:
        exclude |= binary_dilation(wedge_mask, iterations=wedge_dilation_iters)

    j_max_search = int(y_search_max_fraction * Ly / dy)
    if j_max_search < ny:
        exclude[j_max_search:, :] = True

    exclude[:, :3]  = True
    exclude[:, -3:] = True

    slip_masked = slip_score.copy()
    slip_masked[exclude] = 0.0

    j_TP, i_TP = np.unravel_index(np.argmax(slip_masked), slip_masked.shape)
    max_slip = float(slip_masked[j_TP, i_TP])

    if max_slip < slip_threshold:
        return ReflectionDiagnosis(
            type='RR', x_TP=None, y_TP=None,
            h_stem=0.0, h_relative=0.0,
            confidence=min(1.0, 0.5 + (slip_threshold - max_slip) / slip_threshold),
        )

    x_TP_val = float(x[i_TP])
    y_TP_val = float(y[j_TP])

    return ReflectionDiagnosis(
        type='MR',
        x_TP=x_TP_val, y_TP=y_TP_val,
        h_stem=y_TP_val,
        h_relative=y_TP_val / Ly,
        confidence=min(1.0, max_slip / 0.40),
    )


def schlieren_field(rho: np.ndarray, dx: float, dy: float,
                    contrast: float = 15.0) -> np.ndarray:
    """Численный шлирен — exp(-contrast · |∇ρ| / max|∇ρ|)."""
    drhodx = np.gradient(rho, dx, axis=1)
    drhody = np.gradient(rho, dy, axis=0)
    g = np.sqrt(drhodx**2 + drhody**2)
    g /= (g.max() + 1e-12)
    return np.exp(-contrast * g)