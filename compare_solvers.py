import pandas as pd
import matplotlib.pyplot as plt
import os


data_dir = 'results_comparison'  
config_file = 'config.txt'       


solvers = {
    'rusanov.csv': 'Rusanov',
    'hll.csv':     'HLL',
    'hllc.csv':    'HLLC',
    'roe.csv':     'Roe',
    'osher.csv':   'Osher'
}


test_case_num = "unknown"

if os.path.exists(config_file):
    try:
        with open(config_file, 'r') as f:
            for line in f:
                clean_line = line.strip().lower()
                if clean_line.startswith('test_case'):
                    parts = clean_line.split('=')
                    if len(parts) > 1:
                        value_part = parts[1].split('#')[0].strip()
                        test_case_num = value_part
                        break
        print(f"Обнаружен Тест №: {test_case_num}")
    except Exception as e:
        print(f"Не удалось прочитать config.txt: {e}")
else:
    print("Внимание: config.txt не найден, номер теста неизвестен.")



fig, axs = plt.subplots(2, 2, figsize=(15, 10))
fig.suptitle(f'Сравнение приближенных решателей Римана (Тест {test_case_num})', fontsize=16)

vars_map = {
    (0, 0): ('rho', 'Плотность (rho)'),
    (0, 1): ('p',   'Давление (p)'),
    (1, 0): ('u',   'Скорость (u)'),
    (1, 1): ('e',   'Внутренняя энергия (e)')
}


styles = ['--', '-.', ':', '--', '-.']
colors = ['green', 'orange', 'red', 'blue', 'purple'] 


exact_file = os.path.join(data_dir, 'analytical.csv')
if os.path.exists(exact_file):
    try:
        df_exact = pd.read_csv(exact_file)
        if 'x' in df_exact.columns:
            df_exact = df_exact.sort_values('x')
            for (r, c), (var_name, title) in vars_map.items():
                if var_name in df_exact.columns:
                    axs[r, c].plot(df_exact['x'], df_exact[var_name], 'k-', linewidth=1.0, label='Exact', zorder=10)
    except Exception as e:
        print(f"Ошибка чтения аналитики: {e}")
else:
    print(f"Файл аналитики {exact_file} не найден.")

i = 0
for filename, label in solvers.items():
    full_path = os.path.join(data_dir, filename)
    if not os.path.exists(full_path):
        continue
    try:
        df = pd.read_csv(full_path)
        if 'x' in df.columns:
            df = df.sort_values('x')
            
            for (r, c), (var_name, title) in vars_map.items():
                if var_name in df.columns:
                    axs[r, c].plot(df['x'], df[var_name], 
                                   linestyle=styles[i % len(styles)], 
                                   color=colors[i % len(colors)],
                                   linewidth=1.5, 
                                   alpha=0.8,
                                   label=label)
                    
                    axs[r, c].set_title(title)
                    axs[r, c].grid(True, linestyle='--', alpha=0.5)
                    axs[r, c].set_xlabel('x')
        i += 1
    except Exception as e:
        print(f"Ошибка при обработке {filename}: {e}")


axs[0, 0].legend()

plt.tight_layout(rect=[0, 0.03, 1, 0.95])


output_filename = f'comparison_solvers_test_{test_case_num}.png'
plt.savefig(output_filename, dpi=300)
print(f"График сохранен как: {output_filename}")

plt.show()