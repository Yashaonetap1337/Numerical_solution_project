#include "types.h" 
#include "solver.h"
#include "test_cases.h"
#include "analytical.h"

#include <iostream>
#include <locale>


int main() {
    try {
        // чтение конфигурации
        Config cfg = readConfig("config.txt");
        std::cout << "Config loaded successfully. CFL = " << cfg.grid.CFL << std::endl;

        setup_test_case(cfg);
        std::cout << "=== CONFIG DEBUG ===" << std::endl;
        std::cout << "x_min=" << cfg.grid.x_min << " x_max=" << cfg.grid.x_max << std::endl;
        std::cout << "x_diaphragm=" << cfg.grid.x_diaphragm << std::endl;
        std::cout << "Nx=" << cfg.grid.Nx << std::endl;

        // запуск численной симуляции 
        run_simulation(cfg, "numerical_solution.csv");

        // Теперь аналитические снимки создаются автоматически в run_simulation
        std::cout << "All analytical snapshots were generated during simulation." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}