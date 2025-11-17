#include "analytical.h"
#include "riemann_solver.h"
#include "euler_utils.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <filesystem>



void generate_analytical_solution(const Config& cfg, double t, const std::string& outputFilename) {
    std::ofstream out(outputFilename);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to create file: " + outputFilename);
    }

    // Используем ту же сетку, что и в численном решении
    const int points = cfg.grid.Nx;
    const double dx = (cfg.grid.x_max - cfg.grid.x_min) / points;

    out << "x,rho,u,p,e\n";

    for (int i = 0; i <= points; ++i) {
        double x = cfg.grid.x_min + i * dx;

        // Вычисляем безразмерную координату xi = (x - x_diaphragm) / t
        double xi = 0.0;
        if (t > 1e-9) {
            xi = (x - cfg.grid.x_diaphragm) / t;
        }

        // Решаем задачу Римана для этой точки
        State analytical_state = solve_general_riemann_problem(
            cfg.phys.left,
            cfg.phys.right,
            xi,
            cfg,
            cfg.approx_type
        );

        // Вычисляем внутреннюю энергию
        double e = analytical_state.p / (analytical_state.rho * (cfg.phys.gamma - 1.0));

        out << x << ","
            << analytical_state.rho << ","
            << analytical_state.u << ","
            << analytical_state.p << ","
            << e << "\n";
    }

    out.close();
}

void save_analytical_solution(const Config& cfg, const std::string& outputFilename) {
    generate_analytical_solution(cfg, cfg.grid.t_final, outputFilename);
}

void generate_analytical_snapshot(const Config& cfg, double t, const std::string& filename) {
    // Создаем директорию, если она не существует
    std::string directory = filename.substr(0, filename.find_last_of('/'));
    if (!directory.empty()) {
        std::filesystem::create_directories(directory);
    }

    generate_analytical_solution(cfg, t, filename);
}