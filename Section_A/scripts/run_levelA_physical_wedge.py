import sys
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shock_reflection import setup_wedge_simulation, detect_reflection, detachment_wedge_angle, von_neumann_wedge_angle
M = 3.0
THETAS = [15.0, 25.0, 35.0]
NX, NY = (200, 100)
LX, LY = (10.0, 5.0)
X_CORNER = 0.5
WEDGE_LENGTH = None
CFL = 0.4
N_STEPS = 3600
OUTPUT_DIR = './results_levelA_wedge'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def run_single(M, theta_w, nx, ny, Lx, Ly, x_corner, wedge_length, n_steps, cfl, verbose=True):
    if verbose:
        print(f"\n{'=' * 65}")
        print(f'  M={M}, theta_w={theta_w}°, сетка {nx}×{ny}, x_corner={x_corner}, wedge_length={wedge_length}')
        print(f"{'=' * 65}")
    sim = setup_wedge_simulation(M=M, theta_w_deg=theta_w, nx=nx, ny=ny, Lx=Lx, Ly=Ly, x_corner=x_corner, wedge_length=wedge_length, cfl=cfl)
    if verbose:
        print(f'  Ячеек клина: {sim.wedge_mask.sum()}')
        print(f'  Запускаю {n_steps} шагов решателя...')
    t_start = time.time()
    n_done = 0
    while n_done < n_steps:
        n_block = min(200, n_steps - n_done)
        sim.advance_n_steps(n_block)
        n_done += n_block
        if verbose and n_done % 400 == 0:
            W = sim.get_primitive()
            print(f'    {n_done}/{n_steps}: rho∈[{W[0].min():.2f},{W[0].max():.2f}], t_sim={sim.t:.4f}, wall={time.time() - t_start:.0f}s')
    elapsed = time.time() - t_start
    W = sim.get_primitive()
    wedge_mask_inner = sim.wedge_mask[2:-2, 2:-2]
    return {'M': M, 'theta_w': theta_w, 'nx': nx, 'ny': ny, 'Lx': Lx, 'Ly': Ly, 'x_corner': x_corner, 'wedge_length': sim.wedge_length, 'rho': W[0], 'u': W[1], 'v': W[2], 'p': W[3], 'x': sim.x, 'y': sim.y, 'wedge_mask': wedge_mask_inner, 't_final': sim.t, 'n_steps': n_steps, 'wall_time': elapsed}

def visualize_result(result, output_path):
    rho = result['rho']
    p = result['p']
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
    extent = [0.0, Lx, 0.0, Ly]
    fig, ax = plt.subplots(figsize=(11, 4.5))
    im = ax.imshow(rho_plot, origin='lower', extent=extent, cmap=cmap, aspect='equal', vmin=np.nanmin(rho_plot), vmax=np.nanmax(rho_plot))
    theta_rad = np.deg2rad(theta_w)
    x_back = x_corner + wedge_length
    y_back = Ly - np.tan(theta_rad) * wedge_length
    ax.plot([x_corner, x_back], [Ly, y_back], 'k-', lw=1.5)
    ax.plot([x_back, x_back], [y_back, Ly], 'k-', lw=1.5)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_xlabel('$x$', fontsize=12)
    ax.set_ylabel('$y$', fontsize=12)
    ax.set_title(f'$M = {M_val:g}$,  $\\theta_w = {theta_w:g}^\\circ$', fontsize=13)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label('$\\rho / \\rho_0$', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'    Картинка сохранена: {output_path}')
    diag = detect_reflection(rho, x, y, p=p, wedge_mask=mask)
    return diag

def main():
    print('=' * 70)
    print('УРОВЕНЬ A: M=3, три угла клина — ФИЗИЧЕСКИЙ КЛИН (immersed boundary)')
    print('=' * 70)
    theta_d = detachment_wedge_angle(M)
    theta_vN = von_neumann_wedge_angle(M)
    print(f'\nТеоретические критические углы для M={M}:')
    print(f'  theta_w^N (von Neumann) = {theta_vN:.3f}°  — ниже только RR')
    print(f'  theta_w^D (detachment)  = {theta_d:.3f}°  — выше только MR')
    print(f'  Зона двойного решения: {theta_vN:.2f}° < theta_w < {theta_d:.2f}° (ширина {theta_d - theta_vN:.2f}°)')
    print()
    print(f'Тестируем углы: {THETAS}')
    print(f'Ожидание (по теории, cold-start):')
    for tw in THETAS:
        if tw < theta_vN:
            expected = 'RR (только)'
        elif tw < theta_d:
            expected = 'MR (cold-start в зоне двойного решения)'
        else:
            expected = 'MR (отрыв ударной волны от клина)'
        print(f'  theta_w={tw}°: {expected}')
    results = []
    for theta_w in THETAS:
        result = run_single(M, theta_w, NX, NY, LX, LY, X_CORNER, WEDGE_LENGTH, N_STEPS, CFL)
        out_path = os.path.join(OUTPUT_DIR, f'M{M:.0f}_theta{theta_w:.0f}.png')
        diag = visualize_result(result, out_path)
        result['diagnosis'] = diag
        results.append(result)
    print()
    print('=' * 70)
    print('СВОДКА РЕЗУЛЬТАТОВ')
    print('=' * 70)
    print(f"{'theta_w':>8} {'тип':>7} {'h_TP/Ly':>10} {'ожидание':>32}  совпадение")
    print('-' * 70)
    for r in results:
        d = r['diagnosis']
        tw = r['theta_w']
        if tw < theta_vN:
            expected = 'RR'
        elif tw < theta_d:
            expected = 'MR'
        else:
            expected = 'MR'
        match = '+' if d.type in expected else '?'
        print(f'{tw:>7.1f}° {d.type:>7} {d.h_relative:>10.3f} {expected:>32}  {match}')
    print()
    print(f'Все картинки в: {OUTPUT_DIR}/')
if __name__ == '__main__':
    main()
