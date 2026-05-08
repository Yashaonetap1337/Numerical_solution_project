#!/usr/bin/env python
"""
УРОВЕНЬ C (loop): прохождение по конкретному пути углов с substep-walking.

Между milestone-углами симуляция движется не одним прыжком, а серией мелких
подшагов Δθ_substep с короткой релаксацией после каждого. Это:
    1. Делает released-стрип тоньше одной ячейки (orphan'ов нет).
    2. Сохраняет «инерцию» состояния между milestone'ами лучше:
       система плавно отслеживает квази-равновесие.
Картинки сохраняются ТОЛЬКО в milestone-точках (ваш ANGLE_PATH).
"""
import sys
import os
import time
import csv
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from run_levelA import visualize_result

from shock_reflection import (
    setup_wedge_simulation, update_wedge_angle,
    detect_reflection,
    detachment_wedge_angle, von_neumann_wedge_angle,
)


# ============================================================================
# ПАРАМЕТРЫ
# ============================================================================

M = 5.0

# Milestone-углы (ваш путь). Картинки сохраняются ТОЛЬКО здесь.
ANGLE_PATH = [
    18.00,   # старт (RR-зона: 18 < theta_w^N = 20.87)
    22.15,   # forward, в зоне двойного решения
    24.00,   # forward
    27.50,   # forward, под самым theta_w^D
    30.00,   # пик, глубоко в MR
    27.50,   # backward
    24.00,   # backward
    22.15,   # backward
    18.00,   # возврат (ожидается RR)
]

# Substep-walking: между milestone'ами делаем мелкие шаги
DTHETA_SUBSTEP        = 0.5      # шаг по углу при substep-walking
N_RELAX_PER_SUBSTEP   = 200      # релаксация после каждого substep'а
N_RELAX_AT_MILESTONE  = 1000     # доп. релаксация в milestone'е перед снимком

# Численные параметры
NX, NY    = 200, 100
LX, LY    = 10.0, 5.0
X_CORNER  = 0.5
WEDGE_LENGTH = 4.0
CFL = 0.4

N_INITIAL_SETTLE = 3000

OUTPUT_DIR = './results_levelC_loop'
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================================
# ВСПОМОГАТЕЛЬНОЕ
# ============================================================================

def step_direction(milestone_idx, n_total):
    """forward для шагов до пика, backward после."""
    n_peak = int(np.argmax(ANGLE_PATH))
    if milestone_idx <= n_peak:
        return 'forward'
    else:
        return 'backward'


def walk_to_angle(sim, theta_target, theta_current,
                   dtheta_sub, n_relax_sub, verbose=False):
    """
    Квази-статический поворот клина от theta_current до theta_target через
    серию мелких substep'ов размера dtheta_sub. После каждого substep'а
    запускается n_relax_sub шагов решателя.
    """
    delta = theta_target - theta_current
    if abs(delta) < 1e-9:
        return  # ничего делать
    n_sub = max(1, int(np.ceil(abs(delta) / dtheta_sub)))

    if verbose:
        print(f"    walk: {theta_current:.2f}° -> {theta_target:.2f}°  "
              f"({n_sub} substeps по {abs(delta)/n_sub:.3f}°)")

    for k in range(1, n_sub + 1):
        frac = k / n_sub
        theta_intermediate = theta_current + frac * delta
        update_wedge_angle(sim, theta_intermediate)
        sim.advance_n_steps(n_relax_sub)


def take_snapshot(sim, theta, milestone_idx, direction,
                   output_dir, log_rows):
    """Снимок поля + запись в лог. Возвращает (png_path, diag)."""
    W = sim.get_primitive()
    rho, u, v, p = W
    wedge_mask_inner = sim.wedge_mask[2:-2, 2:-2]
    diag = detect_reflection(rho, sim.x, sim.y,
                             p=p, wedge_mask=wedge_mask_inner)

    result = {
        'M': M, 'theta_w': float(theta),
        'rho': rho, 'u': u, 'v': v, 'p': p,
        'x': sim.x, 'y': sim.y,
        'wedge_mask': wedge_mask_inner,
        'wedge_length': sim.wedge_length,
        'x_corner': sim.wedge_x_corner,
        'Lx': LX, 'Ly': LY, 'nx': NX, 'ny': NY,
    }

    png_path = os.path.join(
        output_dir,
        f'step{milestone_idx:02d}_{direction}_theta{theta:g}.png'
    )
    visualize_result(result, png_path)

    log_rows.append({
        'step': milestone_idx,
        'theta_w': f"{theta:.4f}",
        'direction': direction,
        'type': diag.type,
        'x_TP': '' if diag.x_TP is None else f"{diag.x_TP:.4f}",
        'y_TP': '' if diag.y_TP is None else f"{diag.y_TP:.4f}",
        'h_stem': f"{diag.h_stem:.4f}",
        'h_relative': f"{diag.h_relative:.4f}",
        'confidence': f"{diag.confidence:.3f}",
    })

    return diag


# ============================================================================
# ОСНОВНОЙ ЦИКЛ
# ============================================================================

def run_loop():
    print(f"\n{'='*72}")
    print(f"  УРОВЕНЬ C (loop) с substep-walking, M = {M}")
    print(f"  Milestones: " + " -> ".join(f"{a:g}" for a in ANGLE_PATH))
    print(f"  DTHETA_SUBSTEP = {DTHETA_SUBSTEP}, "
          f"N_RELAX_PER_SUBSTEP = {N_RELAX_PER_SUBSTEP}")
    print(f"  N_RELAX_AT_MILESTONE = {N_RELAX_AT_MILESTONE}, "
          f"N_INITIAL_SETTLE = {N_INITIAL_SETTLE}")
    print(f"  Сетка: {NX}x{NY}, область {LX}x{LY}")
    print(f"{'='*72}")

    theta_d  = detachment_wedge_angle(M)
    theta_vN = von_neumann_wedge_angle(M)
    print(f"\n  Теория: theta_w^N = {theta_vN:.3f}°,  "
          f"theta_w^D = {theta_d:.3f}°")
    print(f"  Зона двойного решения: "
          f"[{theta_vN:.2f}°, {theta_d:.2f}°]\n")

    # ----- Cold-start -----
    sim = setup_wedge_simulation(
        M=M, theta_w_deg=ANGLE_PATH[0],
        nx=NX, ny=NY, Lx=LX, Ly=LY,
        x_corner=X_CORNER, wedge_length=WEDGE_LENGTH,
        cfl=CFL,
    )
    print(f"  Initial settle при theta = {ANGLE_PATH[0]}° "
          f"({N_INITIAL_SETTLE} шагов)...")
    t0 = time.time()
    sim.advance_n_steps(N_INITIAL_SETTLE)
    print(f"    готово за {time.time() - t0:.0f} s, t_sim = {sim.t:.3f}\n")

    log_rows = []

    # Snapshot для milestone 1 (стартовая точка)
    direction = step_direction(0, len(ANGLE_PATH))
    diag = take_snapshot(sim, ANGLE_PATH[0], 1, direction,
                          OUTPUT_DIR, log_rows)
    x_str = f"{diag.x_TP:.3f}" if diag.x_TP is not None else "  --  "
    y_str = f"{diag.y_TP:.3f}" if diag.y_TP is not None else "  --  "
    print(f"  step  1  [{direction:>8}]  theta = {ANGLE_PATH[0]:6.2f}°  "
          f"->  {diag.type:>4}  (x_TP = {x_str}, y_TP = {y_str}, "
          f"h/Ly = {diag.h_relative:.3f})")

    # ----- Прогон по milestone'ам -----
    t_loop_start = time.time()

    for ms_idx in range(1, len(ANGLE_PATH)):
        theta_target  = ANGLE_PATH[ms_idx]
        theta_current = ANGLE_PATH[ms_idx - 1]

        t_walk_start = time.time()

        # Substep-walking от текущего угла до целевого
        walk_to_angle(sim, theta_target, theta_current,
                      DTHETA_SUBSTEP, N_RELAX_PER_SUBSTEP)

        # Доп. релаксация в milestone'е для полного равновесия
        sim.advance_n_steps(N_RELAX_AT_MILESTONE)

        t_walk = time.time() - t_walk_start

        # Snapshot
        direction = step_direction(ms_idx, len(ANGLE_PATH))
        diag = take_snapshot(sim, theta_target, ms_idx + 1, direction,
                              OUTPUT_DIR, log_rows)

        x_str = f"{diag.x_TP:.3f}" if diag.x_TP is not None else "  --  "
        y_str = f"{diag.y_TP:.3f}" if diag.y_TP is not None else "  --  "
        print(f"  step {ms_idx+1:>2}  [{direction:>8}]  "
              f"theta = {theta_target:6.2f}°  ->  {diag.type:>4}  "
              f"(x_TP = {x_str}, y_TP = {y_str}, "
              f"h/Ly = {diag.h_relative:.3f},  "
              f"walk+relax = {t_walk:.0f} s)")

    print(f"\n  Полный цикл: {time.time() - t_loop_start:.0f} s")

    # ----- CSV -----
    csv_path = os.path.join(OUTPUT_DIR, 'triple_point_log.csv')
    fieldnames = ['step', 'theta_w', 'direction', 'type',
                  'x_TP', 'y_TP', 'h_stem', 'h_relative', 'confidence']
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(log_rows)
    print(f"\n  Лог тройной точки: {csv_path}")

    return log_rows, theta_vN, theta_d


# ============================================================================
# СВОДКА ПО ГИСТЕРЕЗИСУ
# ============================================================================

def print_hysteresis_summary(log_rows, theta_vN, theta_d):
    print("\n" + "="*72)
    print("СРАВНЕНИЕ ПОВТОРНО ПОСЕЩАЕМЫХ УГЛОВ (тест на гистерезис)")
    print("="*72)
    print(f"  Зона двойного решения: "
          f"theta_w^N = {theta_vN:.2f}° < theta < theta_w^D = {theta_d:.2f}°\n")

    by_angle = {}
    for row in log_rows:
        theta = float(row['theta_w'])
        by_angle.setdefault(theta, []).append(row)

    duplicates = {th: rows for th, rows in by_angle.items() if len(rows) > 1}

    if not duplicates:
        print("  (нет повторно посещаемых углов)")
        return

    print(f"  {'theta_w':>8}  {'forward':>10}  {'backward':>10}  "
          f"{'есть гистерезис?':>18}")
    print("  " + "-"*52)
    for theta in sorted(duplicates.keys()):
        rows = duplicates[theta]
        fwd = next((r for r in rows if r['direction'] == 'forward'), None)
        bwd = next((r for r in rows if r['direction'] == 'backward'), None)
        fwd_type = fwd['type'] if fwd else '?'
        bwd_type = bwd['type'] if bwd else '?'
        in_dual_zone = (theta_vN < theta < theta_d)
        if fwd_type != bwd_type and in_dual_zone:
            verdict = 'ДА (RR vs MR)'
        elif fwd_type == bwd_type:
            verdict = 'нет (одинаковы)'
        else:
            verdict = '?'
        print(f"  {theta:>7.2f}°  {fwd_type:>10}  {bwd_type:>10}  "
              f"{verdict:>18}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    log_rows, theta_vN, theta_d = run_loop()
    print_hysteresis_summary(log_rows, theta_vN, theta_d)
    print(f"\n  Все картинки в: {OUTPUT_DIR}/")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Уровень C: прохождение по конкретному пути углов '
                    'с substep-walking'
    )
    parser.add_argument('--quick', action='store_true',
                        help='Быстрый режим (для отладки)')
    args = parser.parse_args()

    if args.quick:
        N_INITIAL_SETTLE     = 1000
        N_RELAX_PER_SUBSTEP  = 100
        N_RELAX_AT_MILESTONE = 300
        DTHETA_SUBSTEP       = 1.0
        print(f"[QUICK MODE]: settle={N_INITIAL_SETTLE}, "
              f"sub={N_RELAX_PER_SUBSTEP}, ms={N_RELAX_AT_MILESTONE}, "
              f"dtheta_sub={DTHETA_SUBSTEP}")

    main()