"""
Верификация теста CJ детонации (Мейдер, Test 2)
Сравнение численного решения с аналитикой ZND.

Запуск:
    python verify_cj.py results/numerical_solution.csv

Или для отслеживания фронта:
    python verify_cj.py --front results/front_position.csv
"""

import sys
import math
import csv

# ============================================================
# Аналитические параметры (идеальный газ, γ=3)
# ============================================================
GAMMA  = 3.0
RHO0   = 1.84    # г/см³
D_CJ   = 0.88    # см/мкс

P_CJ   = D_CJ**2 * RHO0 / (GAMMA + 1)
RHO_CJ = RHO0 * (GAMMA + 1) / GAMMA
U_CJ   = D_CJ / (GAMMA + 1)
P_VN   = 2 * RHO0 * D_CJ**2 / (GAMMA + 1)   # пик фон Неймана = 2*P_CJ

print("=" * 55)
print("  АНАЛИТИЧЕСКИЕ ЗНАЧЕНИЯ (ZND, γ=3)")
print("=" * 55)
print(f"  D_CJ      = {D_CJ:.4f}  см/мкс")
print(f"  P_CJ      = {P_CJ:.6f} Мбар  (плато за фронтом)")
print(f"  P_vN      = {P_VN:.6f} Мбар  (пик Фон Неймана)")
print(f"  ρ_CJ      = {RHO_CJ:.4f}  г/см³")
print(f"  u_CJ      = {U_CJ:.4f}  см/мкс")
print()

if len(sys.argv) < 2:
    print("Передайте путь к файлу результатов.")
    sys.exit(0)

if "--front" in sys.argv:
    # =========================================================
    # Режим 1: анализ положения фронта x_front(t)
    # =========================================================
    fname = sys.argv[sys.argv.index("--front") + 1]
    times, fronts, pressures = [], [], []
    with open(fname) as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["t"]))
            fronts.append(float(row["x_front"]))
            pressures.append(float(row["P_max"]))

    print("=" * 55)
    print("  АНАЛИЗ СКОРОСТИ ДЕТОНАЦИОННОГО ФРОНТА")
    print("=" * 55)

    # Линейная регрессия x_front(t)
    n = len(times)
    if n > 2:
        t0 = times[0]; x0 = fronts[0]
        t1 = times[-1]; x1 = fronts[-1]
        D_num = (x1 - x0) / (t1 - t0)
        err   = abs(D_num - D_CJ) / D_CJ * 100
        print(f"  D_computed = {D_num:.4f} см/мкс")
        print(f"  D_CJ       = {D_CJ:.4f} см/мкс")
        print(f"  Отклонение = {err:.2f}%  ", end="")
        print("✓ ПРОЙДЕНО" if err < 2.0 else "✗ ПРОВАЛЕНО (норма < 2%)")
        print()

    # Давление на фронте
    P_max_mean = sum(pressures[n//2:]) / len(pressures[n//2:])
    err_p = abs(P_max_mean - P_CJ) / P_CJ * 100
    print(f"  P_max (ср.) = {P_max_mean:.6f} Мбар")
    print(f"  P_CJ        = {P_CJ:.6f} Мбар")
    print(f"  Отклонение  = {err_p:.2f}%  ", end="")
    print("(для тонкого фронта P_max ≈ P_vN = 2·P_CJ)")

else:
    # =========================================================
    # Режим 2: анализ финального профиля P(x), ρ(x), W(x)
    # =========================================================
    fname = sys.argv[1]
    rows = []
    with open(fname) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})

    # Берём срез по y≈0.5 (центральный j)
    y_mid = sorted(set(r["y"] for r in rows))[len(set(r["y"] for r in rows))//2]
    profile = sorted([r for r in rows if abs(r["y"] - y_mid) < 1e-6],
                     key=lambda r: r["x"])

    xs   = [r["x"]   for r in profile]
    rhos = [r["rho"] for r in profile]
    us   = [r["u"]   for r in profile]
    ps   = [r["p"]   for r in profile]

    P_max = max(ps)
    x_front = xs[ps.index(P_max)]

    print("=" * 55)
    print("  АНАЛИЗ ФИНАЛЬНОГО ПРОФИЛЯ")
    print("=" * 55)
    print(f"  Положение фронта: x = {x_front:.4f} см")
    print(f"  P_max            = {P_max:.6f} Мбар")
    print(f"  P_CJ             = {P_CJ:.6f} Мбар")
    print(f"  P_vN             = {P_VN:.6f} Мбар")
    print()

    # Проверяем плато за фронтом (x < x_front - 0.05)
    behind = [r for r in profile if r["x"] < x_front - 0.05]
    if behind:
        P_plateau = sum(r["p"] for r in behind) / len(behind)
        rho_plateau = sum(r["rho"] for r in behind) / len(behind)
        u_plateau   = sum(r["u"]   for r in behind) / len(behind)
        err_P = abs(P_plateau - P_CJ) / P_CJ * 100
        err_u = abs(u_plateau - U_CJ) / U_CJ * 100
        err_r = abs(rho_plateau - RHO_CJ) / RHO_CJ * 100
        print("  Плато за фронтом (продукты взрыва):")
        print(f"    P    = {P_plateau:.6f}  (аналит. {P_CJ:.6f}, откл. {err_P:.1f}%)")
        print(f"    ρ    = {rho_plateau:.4f}     (аналит. {RHO_CJ:.4f},  откл. {err_r:.1f}%)")
        print(f"    u    = {u_plateau:.4f}     (аналит. {U_CJ:.4f},  откл. {err_u:.1f}%)")
        print()
        ok = err_P < 5 and err_u < 5 and err_r < 5
        print("  Критерий плато:", "✓ ПРОЙДЕНО" if ok else "✗ ПРОВАЛЕНО")

    # Профиль — вывод таблицей вблизи фронта
    print()
    print("  Профиль вблизи фронта (±10 ячеек):")
    print(f"  {'x':>8} {'P':>12} {'rho':>10} {'u':>10}")
    print("  " + "-"*44)
    for r in profile:
        if abs(r["x"] - x_front) < 0.06:
            marker = " ←фронт" if abs(r["x"] - x_front) < 0.006 else ""
            print(f"  {r['x']:8.4f} {r['p']:12.6f} {r['rho']:10.4f} "
                  f"{r['u']:10.4f}{marker}")