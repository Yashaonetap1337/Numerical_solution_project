#include "types.h" 
#include "solver.h"
#include "test_cases.h"
#include "analytical.h"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>


std::string get_solver_name(RiemannSolverType type) {
    switch (type) {
    case RiemannSolverType::EXACT:    return "exact";
    case RiemannSolverType::ACOUSTIC: return "acoustic";
    case RiemannSolverType::ROE:      return "roe";
    case RiemannSolverType::RUSANOV:  return "rusanov";
    case RiemannSolverType::OSHER:    return "osher";
    case RiemannSolverType::HLL:      return "hll";
    case RiemannSolverType::HLLC:     return "hllc";
    default:                          return "unknown";
    }
}

int main() {
    try {
        // 1. Читаем конфиг
        Config cfg = readConfig("config.txt");
        setup_test_case(cfg);

        std::cout << "Config loaded. Test Case: " << cfg.phys.test_case << std::endl;

        // 2. Проверяем режим работы
        if (cfg.run_riemann_comparison) {
            // === РЕЖИМ СРАВНЕНИЯ (Запуск 5 методов) ===
            std::cout << ">>> STARTING COMPARISON MODE (5 solvers) <<<" << std::endl;

            std::vector<RiemannSolverType> solvers_to_test = {
                RiemannSolverType::RUSANOV,
                RiemannSolverType::ROE,
                RiemannSolverType::OSHER,
                RiemannSolverType::HLL,
                RiemannSolverType::HLLC
            };

            std::string output_dir = "results_comparison";
            // std::filesystem::create_directories(output_dir);

            // Сохраняем аналитику один раз
            save_analytical_solution(cfg, output_dir + "/analytical.csv");

            for (auto solver_type : solvers_to_test) {
                Config current_cfg = cfg;
                current_cfg.riemann_solver_type = solver_type; // Подмена решателя

                std::string name = get_solver_name(solver_type);
                std::string filename = output_dir + "/" + name + ".csv"; // results_comparison/roe.csv

                std::cout << "  Running: " << name << "..." << std::endl;
                run_simulation(current_cfg, filename);
            }
            std::cout << "Comparison finished. Results in folder: " << output_dir << std::endl;

        }
        else {
            // === ОБЫЧНЫЙ РЕЖИМ (Как было раньше) ===
            std::cout << ">>> STARTING SINGLE SIMULATION <<<" << std::endl;

            // Используем имя из solver.cpp или задаем здесь
            std::string output_file = "numerical_solution.csv";

            // Запускаем с настройками, которые реально прописаны в конфиге
            run_simulation(cfg, output_file);

            // Сохраняем аналитику для сверки
            save_analytical_solution(cfg, "analytical_solution.csv");
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}