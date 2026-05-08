import sys
import os
import time
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_levelA import run_single, visualize_result, NX, NY, LX, LY, X_CORNER, WEDGE_LENGTH, N_STEPS, CFL
from shock_reflection import detachment_wedge_angle, von_neumann_wedge_angle
ANGLE_GRIDS = {3.0: [16.0, 18.0, 20.0, 21.0, 22.0, 24.0, 26.0, 28.0], 5.0: [16.0, 18.0, 20.0, 21.0, 22.0, 24.0, 26.0, 28.0]}
OUTPUT_DIR = './results_levelB'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def run_at_angle(M, theta_w, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n{'=' * 70}")
    print(f'  ОДИН ЭКСПЕРИМЕНТ:  M = {M}, theta_w = {theta_w}°')
    print(f"{'=' * 70}")
    theta_d = detachment_wedge_angle(M)
    theta_vN = von_neumann_wedge_angle(M)
    if theta_vN is not None:
        print(f'  theory: theta_w_N = {theta_vN:.3f} deg')
    else:
        print(f'  theory: theta_w_N — не существует (M < 2.202)')
    if theta_d is not None:
        print(f'  theory: theta_w_D = {theta_d:.3f} deg')
    if theta_vN is None:
        expected = '?'
    elif theta_w < theta_vN:
        expected = 'RR'
    else:
        expected = 'MR'
    print(f'  expected (cold-start): {expected}')
    result = run_single(M, theta_w, NX, NY, LX, LY, X_CORNER, WEDGE_LENGTH, N_STEPS, CFL, verbose=True)
    out_path = os.path.join(output_dir, f'M{int(M)}_theta{theta_w:g}.png')
    diag = visualize_result(result, out_path)
    result['diagnosis'] = diag
    match = '+' if diag.type == expected else '?'
    print(f'\n  detector: type = {diag.type}, h_TP/Ly = {diag.h_relative:.3f}')
    print(f'  match expected ({expected}): {match}')
    return result

def run_sweep(M, theta_list, output_subdir):
    output_dir = os.path.join(OUTPUT_DIR, output_subdir)
    os.makedirs(output_dir, exist_ok=True)
    theta_d = detachment_wedge_angle(M)
    theta_vN = von_neumann_wedge_angle(M)
    print(f"\n{'#' * 70}")
    print(f'#  ПРОГОНКА для M = {M}')
    if theta_vN is not None:
        print(f'#  theory: theta_w_N = {theta_vN:.3f} deg')
    if theta_d is not None:
        print(f'#  theory: theta_w_D = {theta_d:.3f} deg')
    print(f'#  grid sweep: {theta_list}')
    if theta_vN is not None:
        print(f'#  ожидание (cold-start): theta_cr ≈ theta_w_N = {theta_vN:.2f} deg')
    print(f"{'#' * 70}")
    results = []
    t_total = time.time()
    for theta_w in theta_list:
        result = run_single(M, theta_w, NX, NY, LX, LY, X_CORNER, WEDGE_LENGTH, N_STEPS, CFL, verbose=False)
        out_path = os.path.join(output_dir, f'M{int(M)}_theta{theta_w:g}.png')
        diag = visualize_result(result, out_path)
        result['diagnosis'] = diag
        results.append(result)
        wall = result['wall_time']
        print(f'  theta = {theta_w:5.1f} deg  ->  {diag.type:>4}  (h_TP/Ly = {diag.h_relative:.3f}, wall = {wall:.0f} s)')
    print(f'\n  Всего на M={M}: {time.time() - t_total:.0f} s')
    return results

def find_theta_cr(results):
    rr_angles = sorted((r['theta_w'] for r in results if r['diagnosis'].type == 'RR'))
    mr_angles = sorted((r['theta_w'] for r in results if r['diagnosis'].type == 'MR'))
    last_rr = rr_angles[-1] if rr_angles else None
    first_mr = mr_angles[0] if mr_angles else None
    if not rr_angles and (not mr_angles):
        return (None, None, None, None, 'no_data')
    if not rr_angles:
        return (None, first_mr, None, None, 'all_mr')
    if not mr_angles:
        return (last_rr, None, None, None, 'all_rr')
    if last_rr >= first_mr:
        return (last_rr, first_mr, None, None, 'inconsistent')
    theta_cr = (last_rr + first_mr) / 2.0
    uncertainty = (first_mr - last_rr) / 2.0
    return (last_rr, first_mr, theta_cr, uncertainty, 'ok')

def plot_phase_diagram(experimental_points, output_path, M_range=(2.0, 7.0), npts=120):
    M_arr = np.linspace(M_range[0], M_range[1], npts)
    theta_d_arr = np.array([detachment_wedge_angle(M) or np.nan for M in M_arr])
    theta_vN_arr = np.array([von_neumann_wedge_angle(M) or np.nan for M in M_arr])
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(M_arr, theta_d_arr, 'r-', lw=2, label='detachment: $\\theta_w^D(M)$')
    ax.plot(M_arr, theta_vN_arr, 'b-', lw=2, label='von Neumann: $\\theta_w^N(M)$')
    valid = ~(np.isnan(theta_d_arr) | np.isnan(theta_vN_arr))
    if np.any(valid):
        ax.fill_between(M_arr[valid], theta_vN_arr[valid], theta_d_arr[valid], alpha=0.15, color='gray', label='зона двойного решения (RR + MR)')
    if experimental_points:
        M_exp = [p[0] for p in experimental_points]
        theta_exp = [p[1] for p in experimental_points]
        err_exp = [p[2] for p in experimental_points]
        ax.errorbar(M_exp, theta_exp, yerr=err_exp, fmt='o', markersize=10, color='black', ecolor='black', capsize=5, capthick=2, label='эксперимент: $\\theta_{cr}$ (cold-start)')
    ax.set_xlim(M_range)
    ax.set_ylim(15.0, 50.0)
    ax.set_xlabel('$M$, безразм.', fontsize=12)
    ax.set_ylabel('$\\theta_w$, град.', fontsize=12)
    ax.set_title('Фазовая диаграмма перехода RR$\\leftrightarrow$MR', fontsize=13)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\nФазовая диаграмма: {output_path}')

def main():
    print('=' * 70)
    print('УРОВЕНЬ B: серия углов для M = 3 и M = 5, поиск theta_cr')
    print('=' * 70)
    experimental_points = []
    summary_lines = []
    for M, theta_list in ANGLE_GRIDS.items():
        results = run_sweep(M, theta_list, output_subdir=f'M{int(M)}')
        last_rr, first_mr, theta_cr, unc, status = find_theta_cr(results)
        theta_vN = von_neumann_wedge_angle(M)
        print(f'\n  ИТОГ для M = {M}:')
        print(f'    status:     {status}')
        if status == 'ok':
            print(f'    last RR:    {last_rr:.1f} deg')
            print(f'    first MR:   {first_mr:.1f} deg')
            print(f'    theta_cr:   {theta_cr:.1f} +/- {unc:.1f} deg')
            if theta_vN is not None:
                print(f'    theta_w_N:  {theta_vN:.3f} deg (ожидание для cold-start)')
                print(f'    отклонение: {abs(theta_cr - theta_vN):.2f} deg')
            experimental_points.append((M, theta_cr, unc))
            summary_lines.append(f'  M = {M}:  theta_cr = {theta_cr:.1f} +/- {unc:.1f} deg (theory theta_w_N = {theta_vN:.2f} deg)')
        else:
            print(f'    Не удалось определить переход (last_rr={last_rr}, first_mr={first_mr}). Расширьте grid.')
            summary_lines.append(f'  M = {M}:  status = {status}')
    if experimental_points:
        out_path = os.path.join(OUTPUT_DIR, 'phase_diagram.png')
        plot_phase_diagram(experimental_points, out_path)
    print()
    print('=' * 70)
    print('СВОДКА')
    print('=' * 70)
    for line in summary_lines:
        print(line)
    print()
    print(f'Все картинки в: {OUTPUT_DIR}/')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Уровень B: серия углов для M=3 и M=5')
    parser.add_argument('--M', type=float, default=None, help='Число Маха для одного эксперимента (вместе с --theta запускает run_at_angle)')
    parser.add_argument('--theta', type=float, default=None, help='Угол клина (град.) для одного эксперимента')
    args = parser.parse_args()
    if args.M is not None and args.theta is not None:
        run_at_angle(args.M, args.theta, output_dir=os.path.join(OUTPUT_DIR, 'check'))
    elif args.M is not None or args.theta is not None:
        parser.error('--M и --theta нужно указывать вместе')
    else:
        main()
