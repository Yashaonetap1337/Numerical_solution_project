import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from PIL import Image
import re

def natural_sort_key(filename):
    """
    Естественная сортировка для правильного порядка файлов
    """
    # Извлекаем числа из имени файла для сортировки
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split(r'(\d+)', filename)]

def create_simple_animation(snapshot_dir='results', output_gif='animation.gif', duration=200):
    """
    Простой скрипт для создания GIF анимации из CSV файлов снимков
    """
    # Получаем все CSV файлы из директории
    csv_files = glob.glob(os.path.join(snapshot_dir, "*.csv"))
    
    if not csv_files:
        print(f"В директории {snapshot_dir} не найдено CSV файлов")
        return False
    
    # Сортируем файлы в правильном порядке
    csv_files.sort(key=lambda x: natural_sort_key(os.path.basename(x)))
    
    # Перемещаем final_state в конец, если он есть
    final_files = [f for f in csv_files if 'final' in os.path.basename(f).lower()]
    other_files = [f for f in csv_files if 'final' not in os.path.basename(f).lower()]
    
    if final_files:
        csv_files = other_files + final_files
    
    print(f"Найдено {len(csv_files)} файлов для анимации")
    
    # Создаем временную папку для кадров
    temp_dir = "temp_frames"
    os.makedirs(temp_dir, exist_ok=True)
    
    frame_files = []
    
    # Создаем кадры
    for i, csv_file in enumerate(csv_files):
        try:
            # Читаем данные
            data = pd.read_csv(csv_file)
            
            # Создаем график
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            
            # Графики для каждой переменной
            axes[0, 0].plot(data['x'], data['rho'], 'b-', linewidth=2)
            axes[0, 0].set_title('Плотность (ρ)')
            axes[0, 0].grid(True)
            
            axes[0, 1].plot(data['x'], data['u'], 'r-', linewidth=2)
            axes[0, 1].set_title('Скорость (u)')
            axes[0, 1].grid(True)
            
            axes[1, 0].plot(data['x'], data['p'], 'g-', linewidth=2)
            axes[1, 0].set_title('Давление (p)')
            axes[1, 0].grid(True)
            
            axes[1, 1].plot(data['x'], data['e'], 'm-', linewidth=2)
            axes[1, 1].set_title('Внутренняя энергия (e)')
            axes[1, 1].grid(True)
            
            # Добавляем информацию о времени и имени файла
            filename = os.path.basename(csv_file)
            time_info = ""
            if 'time' in data.columns:
                time = data['time'].iloc[0]
                time_info = f"Время: {time:.4f}"
            
            fig.suptitle(f'{filename}\n{time_info}', fontsize=14)
            
            plt.tight_layout()
            
            # Сохраняем кадр
            frame_file = os.path.join(temp_dir, f"frame_{i:04d}.png")
            plt.savefig(frame_file, dpi=80, bbox_inches='tight')
            frame_files.append(frame_file)
            plt.close(fig)
            
            
        except Exception as e:
            print(f"Ошибка при обработке {csv_file}: {e}")
            continue
    
    # Создаем GIF
    if frame_files:
        images = [Image.open(frame) for frame in frame_files]
        images[0].save(
            output_gif,
            save_all=True,
            append_images=images[1:],
            duration=duration,
            loop=0
        )
        
        # Удаляем временные файлы
        for frame_file in frame_files:
            os.remove(frame_file)
        os.rmdir(temp_dir)
        
        print(f"GIF анимация создана: {output_gif}")
        return True
    else:
        print("Не удалось создать ни одного кадра")
        return False

if __name__ == "__main__":
    create_simple_animation()