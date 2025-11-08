#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

def main():
    # Чтение CSV с разделителем (указано ',', но если на самом деле ';', замените на sep=';')
    df = pd.read_csv('numerical_solution.csv', sep=',')

    # Проверка структуры (в продакшене можно закомментировать)
    print("Столбцы:", df.columns.tolist())
    print("Первые 5 строк:")
    print(df.head())

    # Создание графиков
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Газодинамические переменные в зависимости от x', fontsize=16)

    # 1. Плотность
    axs[0, 0].plot(df['x'], df['rho'], 'b-', linewidth=1.5)
    axs[0, 0].set_xlabel('x')
    axs[0, 0].set_ylabel('Плотность (dens)')
    axs[0, 0].grid(True)

    # 2. Давление
    axs[0, 1].plot(df['x'], df['p'], 'r-', linewidth=1.5)
    axs[0, 1].set_xlabel('x')
    axs[0, 1].set_ylabel('Давление (pres)')
    axs[0, 1].grid(True)

    # 3. Скорость
    axs[1, 0].plot(df['x'], df['u'], 'g-', linewidth=1.5)
    axs[1, 0].set_xlabel('x')
    axs[1, 0].set_ylabel('Скорость (velx)')
    axs[1, 0].grid(True)

    # 4. Внутренняя энергия
    axs[1, 1].plot(df['x'], df['e'], 'm-', linewidth=1.5)
    axs[1, 1].set_xlabel('x')
    axs[1, 1].set_ylabel('Внутр. энергия (eint)')
    axs[1, 1].grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Сохранение в PNG
    output_file = 'gasdynamics_plot.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"График сохранён как '{output_file}'")

    # Опционально: показать график (закомментируйте, если работаете без GUI)
    # plt.show()

if __name__ == '__main__':
    main()