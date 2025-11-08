##Решатель уравнений для газовой динамики

Здесь будет поясняющая реадмишка которую в последствии будем дополнять менять пинать сосать чу хотим делать крч

#Структура проекта

Хз я че то забыл запушить сюда мб тут описание гптшка сделала

project/
├── Makefile                 # Файл для сборки проекта
├── main.cpp                 # Главная функция программы
├── config.txt              # Пример конфигурационного файла
├── types.h                 # Основные типы данных и структуры
├── solver.h, solver.cpp    # Основной модуль симуляции
├── grid.h, grid.cpp        # Инициализация расчетной сетки
├── boundary_conditions.h, boundary_conditions.cpp  # Граничные условия
├── euler_utils.h, euler_utils.cpp  # Утилиты для преобразований
├── godunov.h, godunov.cpp  # Реализация метода Годунова
├── config_reader.h, config_reader.cpp  # Парсер конфигурации
├── test_cases.h, test_cases.cpp    # Стандартные тестовые случаи
├── analytical.cpp          # Работа с аналитическими решениями
└── toro_*_exact.txt        # Файлы с аналитическими решениями

#Сборка

Клонирование репозитория(хз надо проверить)
git clone https://github.com/Yashaonetap1337/Numerical_solution_project
cd Numerical_solution_project

Сборка проекта
make

Очистка проекта
make clean

Полная пересборка
make rebuild


#Запуск
./solver

Описание функций файлов впадлу потом сделаю наверн

