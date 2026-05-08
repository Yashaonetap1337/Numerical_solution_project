"""
Пакет shock_reflection для исследования отражения косой УВ от клина.

УРОВЕНЬ A: M=3, θ ∈ {15°, 25°, 35°}.
Постановка: ФИЗИЧЕСКИЙ КЛИН в потоке (immersed boundary, по Ben-Dor Fig. 2.30).
"""
from .solver import EulerSim, GAMMA, NGHOST, apply_wedge_bc

# Постановка с физическим клином (Ivanov-Chpoun, Ben-Dor §2.4)
from .setup_wedge import (
    setup_wedge_simulation, update_wedge_angle,
    make_uniform_inflow_bc, make_uniform_initial_condition,
)

from .theory import (
    detachment_wedge_angle, von_neumann_wedge_angle, sonic_wedge_angle,
    beta_1_at_detachment, beta_1_at_von_neumann,
    oblique_deflection, oblique_M2, theta_max_for_M, beta_for_deflection,
    M_0_CRITICAL,
)
from .analysis import detect_reflection, schlieren_field, ReflectionDiagnosis

__all__ = [
    'EulerSim', 'GAMMA', 'NGHOST', 'apply_wedge_bc',
    'setup_wedge_simulation', 'update_wedge_angle',
    'make_uniform_inflow_bc', 'make_uniform_initial_condition',
    'detachment_wedge_angle', 'von_neumann_wedge_angle', 'sonic_wedge_angle',
    'beta_1_at_detachment', 'beta_1_at_von_neumann',
    'oblique_deflection', 'oblique_M2', 'theta_max_for_M', 'beta_for_deflection',
    'M_0_CRITICAL',
    'detect_reflection', 'schlieren_field', 'ReflectionDiagnosis',
]