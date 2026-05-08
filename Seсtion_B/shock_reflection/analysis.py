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

def detect_reflection(rho: np.ndarray, x: np.ndarray, y: np.ndarray, p: Optional[np.ndarray]=None, wedge_mask: Optional[np.ndarray]=None, threshold_factor: float=0.5, min_h_rel: float=0.1) -> ReflectionDiagnosis:
    ny, nx = rho.shape
    if p is None:
        p = rho ** 1.4
    p_bot = p[:5, :].mean(axis=0)
    p_max = p_bot.max()
    p_floor = p_bot[:nx // 6].mean()
    if p_max - p_floor < 0.3:
        return ReflectionDiagnosis(type='unknown', x_TP=None, y_TP=None, h_stem=0.0, h_relative=0.0, confidence=0.3)
    i_peak = int(np.argmax(p_bot))
    x_TP_candidate = float(x[i_peak])
    p_col = p[:, i_peak]
    threshold = p_floor + threshold_factor * (p_max - p_floor)
    h_cells = 0
    for j in range(ny):
        if wedge_mask is not None and wedge_mask[j, i_peak]:
            break
        if p_col[j] > threshold:
            h_cells = j + 1
        else:
            break
    Ly = y[-1] + (y[1] - y[0])
    h_rel = h_cells * (y[1] - y[0]) / Ly
    y_TP_val = y[h_cells - 1] if h_cells > 0 else 0.0
    if h_rel < min_h_rel:
        return ReflectionDiagnosis(type='RR', x_TP=None, y_TP=None, h_stem=0.0, h_relative=0.0, confidence=0.85)
    return ReflectionDiagnosis(type='MR', x_TP=x_TP_candidate, y_TP=float(y_TP_val), h_stem=float(y_TP_val), h_relative=float(h_rel), confidence=1.0)

def schlieren_field(rho: np.ndarray, dx: float, dy: float, contrast: float=15.0) -> np.ndarray:
    drhodx = np.gradient(rho, dx, axis=1)
    drhody = np.gradient(rho, dy, axis=0)
    g = np.sqrt(drhodx ** 2 + drhody ** 2)
    g /= g.max() + 1e-12
    return np.exp(-contrast * g)
