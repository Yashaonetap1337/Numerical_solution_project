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
from shock_reflection import setup_wedge_simulation, detect_reflection, detachment_wedge_angle, von_neumann_wedge_angle
CASES = [{'M': 3.0, 'theta_w': 25.0, 'expected': 'MR', 'note': 'дуальная зона M=3'}, {'M': 3.0, 'theta_w': 28.0, 'expected': 'MR', 'note': 'дуальная зона M=3'}, {'M': 3.0, 'theta_w': 30.0, 'expected': 'MR', 'note': 'дуальная зона M=3, ближе к D'}, {'M': 5.0, 'theta_w': 30.0, 'expected': 'MR', 'note': 'выше theta_w^D=27.77 для M=5'}, {'M': 5.0, 'theta_w': 35.0, 'expected': 'MR', 'note': 'глубоко в MR, M=5'}]
NX, NY = (200, 100)
LX, LY = (10.0, 5.0)
X_CORNER = 0.5
WEDGE_LENGTH = 4.0
CFL = 0.4
N_STEPS = 2500
OUTPUT_DIR = './results_level_c_detector'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def visualize_with_triple_point(result, diag, output_path):
    rho = result['rho']
    x = result['x']
    y = result['y']
    mask = result['wedge_mask']
    Lx, Ly = (result['Lx'], result['Ly'])
    theta_w = result['theta_w']
    M_val = result['M']
    x_corner = result['x_corner']
    wedge_length = result['wedge_length']
    rho_plot = np.where(mask, np.nan, rho)
    cmap = plt.cm.inferno.copy()
    cmap.set_bad(color='white')
    fig, ax = plt.subplots(figsize=(11, 4.5))
    im = ax.imshow(rho_plot, origin='lower', extent=[0.0, Lx, 0.0, Ly], cmap=cmap, aspect='equal', vmin=np.nanmin(rho_plot), vmax=np.nanmax(rho_plot))
    theta_rad = np.deg2rad(theta_w)
    x_back = x_corner + wedge_length
    y_back = Ly - np.tan(theta_rad) * wedge_length
    ax.plot([x_corner, x_back], [Ly, y_back], 'k-', lw=1.5)
    ax.plot([x_back, x_back], [y_back, Ly], 'k-', lw=1.5)
    if diag.type == 'MR' and diag.x_TP is not None:
        ax.plot(diag.x_TP, diag.y_TP, 'o', markersize=12, markerfacecolor='cyan', markeredgecolor='white', markeredgewidth=2, zorder=10, label=f'TP (детектор): ({diag.x_TP:.3f}, {diag.y_TP:.3f})')
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_xlabel('$x$, безразм.', fontsize=12)
    ax.set_ylabel('$y$, безразм.', fontsize=12)
    title = f'$M = {M_val:g}$,  $\\theta_w = {theta_w:g}^\\circ$  '
    title += f'(детектор: {diag.type}, conf={diag.confidence:.2f})'
    ax.set_title(title, fontsize=12)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label('$\\rho / \\rho_0$, безразм.', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'    Картинка сохранена: {output_path}')

def save_fields_csv(rho, u, v, p, x, y, wedge_mask, csv_path):
    ny, nx = rho.shape
    X, Y = np.meshgrid(x, y)
    data = np.column_stack([np.tile(np.arange(nx), ny), np.repeat(np.arange(ny), nx), X.ravel(), Y.ravel(), rho.ravel(), u.ravel(), v.ravel(), p.ravel(), wedge_mask.ravel().astype(int)])
    fmt = ['%d', '%d', '%.4f', '%.4f', '%.6e', '%.6e', '%.6e', '%.6e', '%d']
    np.savetxt(csv_path, data, delimiter=',', fmt=fmt, header='i,j,x,y,rho,u,v,p,in_wedge', comments='')

def run_case(M, theta_w, expected, note, n_steps, output_dir):
    print(f"\n{'=' * 60}")
    print(f'  M = {M}, theta_w = {theta_w}°  ({expected}, {note})')
    print(f"{'=' * 60}")
    sim = setup_wedge_simulation(M=M, theta_w_deg=theta_w, nx=NX, ny=NY, Lx=LX, Ly=LY, x_corner=X_CORNER, wedge_length=WEDGE_LENGTH, cfl=CFL)
    print(f'  Запускаю {n_steps} шагов...')
    t0 = time.time()
    sim.advance_n_steps(n_steps)
    elapsed = time.time() - t0
    print(f'    готово за {elapsed:.0f} s')
    W = sim.get_primitive()
    rho, u, v, p = W
    wedge_mask_inner = sim.wedge_mask[2:-2, 2:-2]
    diag = detect_reflection(rho, sim.x, sim.y, p=p, wedge_mask=wedge_mask_inner)
    result = {'M': M, 'theta_w': float(theta_w), 'rho': rho, 'u': u, 'v': v, 'p': p, 'x': sim.x, 'y': sim.y, 'wedge_mask': wedge_mask_inner, 'wedge_length': sim.wedge_length, 'x_corner': sim.wedge_x_corner, 'Lx': LX, 'Ly': LY, 'nx': NX, 'ny': NY}
    img_path = os.path.join(output_dir, f'M{int(M)}_theta{theta_w:g}_image.png')
    visualize_with_triple_point(result, diag, img_path)
    csv_path = os.path.join(output_dir, f'M{int(M)}_theta{theta_w:g}_fields.csv')
    print(f'  Сохраняю поля в CSV...')
    save_fields_csv(rho, u, v, p, sim.x, sim.y, wedge_mask_inner, csv_path)
    print(f'  Детектор: type={diag.type}, x_TP={diag.x_TP}, y_TP={diag.y_TP}, h/Ly={diag.h_relative:.3f}, conf={diag.confidence:.2f}')
    return {'M': M, 'theta_w': theta_w, 'expected': expected, 'note': note, 'detector_type': diag.type, 'detector_x_TP': diag.x_TP if diag.x_TP is not None else '', 'detector_y_TP': diag.y_TP if diag.y_TP is not None else '', 'detector_h_relative': diag.h_relative, 'detector_confidence': diag.confidence, 'wall_time_sec': elapsed}

def main():
    print('=' * 70)
    print('КАЛИБРОВОЧНЫЕ ПРОГОНКИ для детектора тройной точки')
    print('=' * 70)
    print('\nКейсы:')
    for c in CASES:
        td = detachment_wedge_angle(c['M'])
        tn = von_neumann_wedge_angle(c['M'])
        print(f"  M={c['M']:.0f}, theta={c['theta_w']:5.1f}°  (theta_w^N={tn:.2f}°, theta_w^D={td:.2f}°)  [{c['note']}]")
    print(f'\nСетка: {NX}x{NY},  область {LX}x{LY},  N_steps = {N_STEPS}\n')
    results = []
    t_total = time.time()
    for case in CASES:
        r = run_case(case['M'], case['theta_w'], case['expected'], case['note'], N_STEPS, OUTPUT_DIR)
        results.append(r)
    print(f'\nВсё за {time.time() - t_total:.0f} s')
    summary_path = os.path.join(OUTPUT_DIR, 'summary.csv')
    with open(summary_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['case_id', 'M', 'theta_w', 'expected_type', 'note', 'detector_type', 'detector_x_TP', 'detector_y_TP', 'detector_h_relative', 'detector_confidence', 'true_x_TP', 'true_y_TP', 'comments', 'wall_time_sec'])
        for i, r in enumerate(results):
            writer.writerow([i + 1, r['M'], r['theta_w'], r['expected'], r['note'], r['detector_type'], r['detector_x_TP'], r['detector_y_TP'], f"{r['detector_h_relative']:.4f}", f"{r['detector_confidence']:.3f}", '', '', '', f"{r['wall_time_sec']:.1f}"])
    print(f'\nСводка: {summary_path}')
    print(f'Все файлы в: {OUTPUT_DIR}/')
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--quick', action='store_true', help='Быстрый режим: 1500 шагов')
    args = parser.parse_args()
    if args.quick:
        N_STEPS = 1500
        print(f'[QUICK MODE]: N_STEPS = {N_STEPS}')
    main()
