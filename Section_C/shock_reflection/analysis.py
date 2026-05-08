from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional
from scipy.ndimage import binary_dilation

@dataclass
class ReflectionDiagnosis:
    type: str
    x_TP: Optional[float]
    y_TP: Optional[float]
    h_stem: float
    h_relative: float
    confidence: float = 1.0

def detect_reflection(rho: np.ndarray, x: np.ndarray, y: np.ndarray, p: Optional[np.ndarray]=None, wedge_mask: Optional[np.ndarray]=None, slip_threshold: float=0.2, wedge_dilation_iters: int=4, y_search_max_fraction: float=0.18) -> ReflectionDiagnosis:
    if p is None:
        p = rho ** 1.4
    ny, nx = rho.shape
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    Ly = y[-1] + dy
    drho_dx = np.gradient(rho, dx, axis=1)
    drho_dy = np.gradient(rho, dy, axis=0)
    G_rho = np.sqrt(drho_dx ** 2 + drho_dy ** 2)
    dp_dx = np.gradient(p, dx, axis=1)
    dp_dy = np.gradient(p, dy, axis=0)
    G_p = np.sqrt(dp_dx ** 2 + dp_dy ** 2)
    g_rho_norm = G_rho / (G_rho.max() + 1e-12)
    g_p_norm = G_p / (G_p.max() + 1e-12)
    slip_score = g_rho_norm * (1.0 - g_p_norm)
    exclude = np.zeros_like(slip_score, dtype=bool)
    if wedge_mask is not None:
        exclude |= binary_dilation(wedge_mask, iterations=wedge_dilation_iters)
    j_max_search = int(y_search_max_fraction * Ly / dy)
    if j_max_search < ny:
        exclude[j_max_search:, :] = True
    exclude[:, :3] = True
    exclude[:, -3:] = True
    slip_masked = slip_score.copy()
    slip_masked[exclude] = 0.0
    j_TP, i_TP = np.unravel_index(np.argmax(slip_masked), slip_masked.shape)
    max_slip = float(slip_masked[j_TP, i_TP])
    if max_slip < slip_threshold:
        return ReflectionDiagnosis(type='RR', x_TP=None, y_TP=None, h_stem=0.0, h_relative=0.0, confidence=min(1.0, 0.5 + (slip_threshold - max_slip) / slip_threshold))
    x_TP_val = float(x[i_TP])
    y_TP_val = float(y[j_TP])
    return ReflectionDiagnosis(type='MR', x_TP=x_TP_val, y_TP=y_TP_val, h_stem=y_TP_val, h_relative=y_TP_val / Ly, confidence=min(1.0, max_slip / 0.4))

def schlieren_field(rho: np.ndarray, dx: float, dy: float, contrast: float=15.0) -> np.ndarray:
    drhodx = np.gradient(rho, dx, axis=1)
    drhody = np.gradient(rho, dy, axis=0)
    g = np.sqrt(drhodx ** 2 + drhody ** 2)
    g /= g.max() + 1e-12
    return np.exp(-contrast * g)
