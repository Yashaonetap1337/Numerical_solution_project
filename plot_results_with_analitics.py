import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

def plot_results(numerical_file, analytical_file):

    if not os.path.exists(numerical_file):
        print(f"Ошибка: Файл с численным решением '{numerical_file}' не найден.")
        print("Пожалуйста, сначала запустите C++ симуляцию для его создания.")
        return

    if not os.path.exists(analytical_file):
        print(f"Ошибка: Файл с аналитическим решением '{analytical_file}' не найден.")
        return

    print("Файлы найдены. Загрузка данных...")

    try:
        numerical_data = pd.read_csv(numerical_file)
        analytical_data = pd.read_csv(analytical_file)
        print("Данные успешно загружены.")
    except Exception as e:
        print(f"Произошла ошибка при чтении CSV файлов: {e}")
        return


    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    variables_to_plot = [
        ('rho', 'Плотность (ρ)'),
        ('p', 'Давление (p)'),
        ('u', 'Скорость (u)'),
        ('e', 'Внутренняя энергия (e)')
    ]


    for i, (var, title) in enumerate(variables_to_plot):

        ax = axes.flatten()[i]


        if var in analytical_data.columns:
            # Строим аналитическое решение (сплошная черная линия)
            ax.plot(analytical_data['x'], analytical_data[var], 'k-', label='Аналитическое')
        
        if var in numerical_data.columns:
            # Строим численное решение (синие точки)
            ax.plot(numerical_data['x'], numerical_data[var], 'bo', markersize=3, label='Численное')
        

        ax.set_title(title, fontsize=14)
        ax.set_xlabel('Позиция (x)')
        ax.set_ylabel(title.split(' ')[0]) 
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.legend()
        ax.set_xlim(analytical_data['x'].min(), analytical_data['x'].max())


    fig.suptitle('Сравнение численного и аналитического решений', fontsize=18)
    

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])


    output_filename = 'solution_comparison.png'
    plt.savefig(output_filename)
    print(f"\nГрафик успешно сохранен в файл: {output_filename}")
    

    plt.show()


if __name__ == "__main__":
    # Эта часть позволяет запускать скрипт из командной строки с аргументами,
    # но также работает и без них, используя имена файлов по умолчанию.
    parser = argparse.ArgumentParser(description="Построение графиков для задачи о распаде разрыва.")
    
    parser.add_argument('--num', default='numerical_solution.csv',
                        help='Имя файла с численным решением (по умолчанию: numerical_solution.csv)')
    
    parser.add_argument('--analyt', default='analytical_solution.csv',
                        help='Имя файла с аналитическим решением (по умолчанию: analytical_solution.csv)')

    args = parser.parse_args()

    plot_results(args.num, args.analyt)
