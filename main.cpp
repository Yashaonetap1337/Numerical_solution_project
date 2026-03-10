#include "types.h" 
#include "solver.h"
#include "test_cases.h"
#include "analytical.h"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem> 
#include <chrono> 
#include <fstream> // Для записи в файл

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
        Config cfg = readConfig("config.txt");
        setup_test_case(cfg);

        std::string output_dir = "results";
        // You may need to create the directory manually or use mkdir

        std::string output_file = output_dir + "/numerical_solution.csv";
        auto start = std::chrono::high_resolution_clock::now();
        run_simulation(cfg, output_file);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Simulation time: " << std::chrono::duration<double>(end - start).count() << " s\n";

        save_analytical_solution(cfg, output_dir + "/analytical_solution.csv");

        std::cout << "Done.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}