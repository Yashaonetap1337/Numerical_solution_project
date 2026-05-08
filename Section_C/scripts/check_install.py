#!/usr/bin/env python
"""
Быстрый sanity check (~30 секунд):
1. Проверяет что теоретические значения корректны.
2. Запускает один маленький расчёт с физическим клином чтобы убедиться,
   что solver работает.

Запуск:
    python scripts/check_install.py
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from shock_reflection import (
    detachment_wedge_angle, von_neumann_wedge_angle, beta_1_at_detachment,
    setup_wedge_simulation, detect_reflection,
)


def main():
    print("="*60)
    print("CHECK 1: Теоретические значения")
    print("="*60)
    
    # Проверка известных значений (Ben-Dor стр. 83 + Chpoun et al. 1995)
    checks = [
        ('M=5: theta_w^D ~ 27.77 deg',
         detachment_wedge_angle(5.0), 27.77, 0.05),
        ('M=5: theta_w^N ~ 20.86 deg',
         von_neumann_wedge_angle(5.0), 20.86, 0.05),
        ('M=4.96: beta_1^D = 39.33 deg (Chpoun)',
         beta_1_at_detachment(4.96), 39.33, 0.05),
        ('M=2: theta_w^N = None (M<2.202)',
         von_neumann_wedge_angle(2.0), None, None),
    ]
    
    all_ok = True
    for name, got, expected, tol in checks:
        if expected is None:
            ok = got is None
        else:
            ok = got is not None and abs(got - expected) < tol
        status = '[OK]' if ok else '[FAIL]'
        if expected is None:
            print(f"  {status} {name}: got={got}")
        else:
            print(f"  {status} {name}: got={got:.3f} deg, expected={expected}+-{tol} deg")
        if not ok:
            all_ok = False
    
    if not all_ok:
        print("\n[!] Some theoretical values failed!")
        sys.exit(1)
    
    print("\nAll theoretical values OK.")
    
    print()
    print("="*60)
    print("CHECK 2: Minimum calculation (M=3, theta_w=25 deg, physical wedge)")
    print("="*60)
    
    # Маленький расчёт с физическим клином для проверки solver
    sim = setup_wedge_simulation(
        M=3.0, theta_w_deg=25.0,
        nx=120, ny=60, Lx=3.0, Ly=1.0,
        x_corner=0.5,
        cfl=0.4,
    )
    
    print(f"  Wedge cells: {sim.wedge_mask.sum()}")
    print("  Running 300 solver steps...")
    
    t0 = time.time()
    sim.advance_n_steps(300)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")
    
    W = sim.get_primitive()
    rho = W[0]
    p = W[3]
    rho_max = rho.max()
    rho_min = rho.min()
    print(f"  rho range: [{rho_min:.3f}, {rho_max:.3f}]")
    
    if rho_max < 2.0:
        print("\n[!] Density did not grow - something wrong with calculation!")
        sys.exit(1)
    
    # Проверим тип отражения
    wedge_mask_inner = sim.wedge_mask[2:-2, 2:-2]
    diag = detect_reflection(rho, sim.x, sim.y, p=p, wedge_mask=wedge_mask_inner)
    print(f"  Reflection type: {diag.type}, h_TP/Ly={diag.h_relative:.3f}")
    
    print()
    print("="*60)
    print("ALL CHECKS PASSED [OK]")
    print("="*60)
    print()
    print("You can now run the main experiment:")
    print("  python scripts/run_levelA_physical_wedge.py")


if __name__ == '__main__':
    main()