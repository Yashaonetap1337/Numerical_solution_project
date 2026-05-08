"""
Детектор тройной точки для wedge-in-flow геометрии (физический клин сверху).

Геометрия: клин в верхней части, поток налетает слева, нижняя стенка — 
плоскость симметрии. От вершины клина идёт косая УВ.

RR vs MR в этой геометрии:
    - RR: косая УВ от вершины клина встречает нижнюю стенку в ОДНОЙ ТОЧКЕ P 
          и отражается под симметричным углом.
    - MR: возле нижней стенки УВ "обрезается" вертикальным маховским стержнем.
          Тройная точка T находится НА НЕКОТОРОЙ ВЫСОТЕ над нижней стенкой.

ДЕТЕКТОР (через профиль давления у нижней стенки):
    1. Извлекаем профиль p(x) у y=0 (среднее по 3-5 нижним строкам).
    2. Находим максимум давления p_max (там где УВ ударяется в стенку).
    3. Берём вертикальную колонку x=x_max и идём вверх по y.
    4. Высота стержня h_TP = первое y, где давление падает ниже порога.
    5. Если h_TP/Ly > min_threshold → MR, иначе RR.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional


@dataclass
class ReflectionDiagnosis:
    type: str
    x_TP: Optional[float]
    y_TP: Optional[float]
    h_stem: float
    h_relative: float
    confidence: float = 1.0


def detect_reflection(rho: np.ndarray, x: np.ndarray, y: np.ndarray,
                       p: Optional[np.ndarray] = None,
                       wedge_mask: Optional[np.ndarray] = None,
                       threshold_factor: float = 0.5,
                       min_h_rel: float = 0.10) -> ReflectionDiagnosis:
    """
    Детектор для wedge-in-flow геометрии.
    
    Аргументы:
        rho, p          : поля плотности и давления (без ghost-слоёв)
        x, y            : координаты центров ячеек
        wedge_mask      : маска клина (True где клин)
        threshold_factor: доля от p_max для определения "стержень есть"
        min_h_rel       : минимальная h/Ly для классификации MR
    """
    ny, nx = rho.shape
    
    # Если поле p не передано — приближаем (грубо)
    if p is None:
        p = rho**1.4  # эта аппроксимация пригодится только для unit-тестов
    
    # Профиль давления у нижней стенки (нижние 5 строк)
    p_bot = p[:5, :].mean(axis=0)
    
    p_max = p_bot.max()
    p_floor = p_bot[:nx//6].mean()  # давление в свежем потоке слева
    
    if p_max - p_floor < 0.3:
        return ReflectionDiagnosis(
            type='unknown', x_TP=None, y_TP=None,
            h_stem=0.0, h_relative=0.0, confidence=0.3,
        )
    
    # Колонка с максимальным давлением у нижней стенки — там где УВ ударяет
    i_peak = int(np.argmax(p_bot))
    x_TP_candidate = float(x[i_peak])
    
    # Идём ВВЕРХ по этой колонке
    p_col = p[:, i_peak]
    threshold = p_floor + threshold_factor * (p_max - p_floor)
    
    h_cells = 0
    for j in range(ny):
        # Игнорируем ячейки внутри клина
        if wedge_mask is not None and wedge_mask[j, i_peak]:
            break
        if p_col[j] > threshold:
            h_cells = j + 1
        else:
            break
    
    Ly = y[-1] + (y[1] - y[0])  # полная высота области
    h_rel = h_cells * (y[1] - y[0]) / Ly
    y_TP_val = y[h_cells - 1] if h_cells > 0 else 0.0
    
    if h_rel < min_h_rel:
        return ReflectionDiagnosis(
            type='RR', x_TP=None, y_TP=None,
            h_stem=0.0, h_relative=0.0, confidence=0.85,
        )
    
    return ReflectionDiagnosis(
        type='MR',
        x_TP=x_TP_candidate,
        y_TP=float(y_TP_val),
        h_stem=float(y_TP_val),
        h_relative=float(h_rel),
        confidence=1.0,
    )


def schlieren_field(rho: np.ndarray, dx: float, dy: float,
                     contrast: float = 15.0) -> np.ndarray:
    """Численный шлирен — exp(-contrast · |∇ρ| / max|∇ρ|)."""
    drhodx = np.gradient(rho, dx, axis=1)
    drhody = np.gradient(rho, dy, axis=0)
    g = np.sqrt(drhodx**2 + drhody**2)
    g /= (g.max() + 1e-12)
    return np.exp(-contrast * g)