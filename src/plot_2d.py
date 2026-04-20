import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

def plot_pressure_at_times(csv_files, time_labels=None, output_name='pressure_evolution.png'):
    """
    Построение полей давления для нескольких моментов времени
    
    Parameters:
    -----------
    csv_files : list of str
        Список путей к CSV файлам с результатами (в порядке возрастания времени)
    time_labels : list of str, optional
        Подписи времени для каждого графика (если None, берутся имена файлов)
    output_name : str
        Имя выходного файла с графиком
    """
    
    # Параметры геометрии
    x_corner = 2.0
    y_corner = 3.0
    x_det = 0.15
    
    n_files = len(csv_files)
    if n_files == 0:
        print("Ошибка: не указаны CSV файлы")
        return
    
    # Определяем размеры сетки подграфиков
    if n_files <= 3:
        n_rows, n_cols = 1, n_files
    elif n_files <= 6:
        n_rows, n_cols = 2, 3
    else:
        n_rows, n_cols = (n_files + 2) // 3, 3
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_files == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Скрываем лишние подграфики
    for i in range(n_files, len(axes)):
        axes[i].axis('off')
    
    # Для каждого файла строим график давления
    for idx, csv_file in enumerate(csv_files):
        if not os.path.exists(csv_file):
            print(f"Предупреждение: файл {csv_file} не найден, пропускаем")
            continue
            
        # Загрузка данных
        df = pd.read_csv(csv_file)
        
        # Получаем уникальные координаты для построения сетки
        x_unique = np.sort(df['x'].unique())
        y_unique = np.sort(df['y'].unique())
        X, Y = np.meshgrid(x_unique, y_unique)
        
        # Интерполяция давления на сетку
        from scipy.interpolate import griddata
        p_grid = griddata((df['x'], df['y']), df['p'], (X, Y), method='linear')
        
        # Маска для стенки (rho > 100) и невалидных данных
        if 'rho' in df.columns:
            rho_grid = griddata((df['x'], df['y']), df['rho'], (X, Y), method='linear')
            wall_mask = rho_grid > 100
            invalid_mask = rho_grid < 0.1
        else:
            wall_mask = np.zeros_like(p_grid, dtype=bool)
            invalid_mask = np.isnan(p_grid)
        
        p_plot = np.ma.masked_where(wall_mask | invalid_mask, p_grid)
        
        # Получаем время из данных или используем подпись
        if time_labels and idx < len(time_labels):
            time_str = time_labels[idx]
        elif 't' in df.columns:
            t_val = df['t'].iloc[0]
            time_str = f't = {t_val:.3f}'
        else:
            # Пробуем извлечь время из имени файла
            import re
            match = re.search(r't[=_]?([0-9.]+)', csv_file)
            if match:
                time_str = f't = {match.group(1)}'
            else:
                time_str = os.path.basename(csv_file).replace('.csv', '')
        
        # Рисуем график
        ax = axes[idx]
        im = ax.pcolormesh(X, Y, p_plot, shading='auto', cmap='plasma', 
                          vmin=0, vmax=0.4)  # Диапазон давления для детонации
        
        # Отображаем стенку
        ax.add_patch(Rectangle((0, 0), x_corner, y_corner, 
                               facecolor='gray', edgecolor='black', 
                               alpha=0.4, label='Стенка' if idx == 0 else ""))
        
        # Отмечаем угол и детонатор
        if idx == 0:
            ax.axvline(x=x_corner, color='cyan', linestyle='--', alpha=0.5, linewidth=1, label=f'Угол x={x_corner}')
            ax.axhline(y=y_corner, color='cyan', linestyle='--', alpha=0.5, linewidth=1)
            ax.axvline(x=x_det, color='lime', linestyle=':', alpha=0.7, linewidth=1, label=f'Детонатор x={x_det}')
        else:
            ax.axvline(x=x_corner, color='cyan', linestyle='--', alpha=0.5, linewidth=1)
            ax.axhline(y=y_corner, color='cyan', linestyle='--', alpha=0.5, linewidth=1)
            ax.axvline(x=x_det, color='lime', linestyle=':', alpha=0.7, linewidth=1)
        
        ax.set_xlim(0, 9)
        ax.set_ylim(0, 4)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Давление, {time_str}', fontsize=10)
        
        # Добавляем цветовую шкалу для каждого графика
        plt.colorbar(im, ax=ax, label='p', fraction=0.046, pad=0.04)
    
    # Общая легенда (только для первого графика, если есть место)
    if n_files > 0:
        axes[0].legend(loc='upper right', fontsize=8)
    
    plt.suptitle('Эволюция давления при обтекании угла детонационной волной', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_name, dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"График сохранён как {output_name}")

# Пример использования:
if __name__ == "__main__":
    # Укажите здесь свои CSV файлы в порядке возрастания времени
    csv_files = [
        "snapshots/snapshot_time_0_104338.csv",    # t = 0.5
        "snapshots/snapshot_time_1_301464.csv",    # t = 1.0
        "snapshots/snapshot_time_2_302214.csv",    # t = 2.0
        "snapshots/snapshot_time_3_003281.csv",    # t = 3.0
        "snapshots/snapshot_time_3_903021.csv",    # t = 4.0
    ]
    
    # Можно задать свои подписи времени
    time_labels = ["t = 0.5", "t = 1.0", "t = 2.0", "t = 3.0", "t = 4.0"]
    
    plot_pressure_at_times(csv_files, time_labels, "pressure_evolution.png")