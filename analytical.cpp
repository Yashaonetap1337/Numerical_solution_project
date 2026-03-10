#include "analytical.h"
#include "riemann_solver.h"
#include "euler_utils.h"
#include <fstream>
#include <cmath>

void generate_analytical_solution(const Config& cfg, double t, const std::string& outputFilename) {
    std::ofstream out(outputFilename);
    if (!out.is_open()) throw std::runtime_error("Cannot open " + outputFilename);

    const int Nx = cfg.grid.Nx;
    const int Ny = cfg.grid.Ny;
    double dx = (cfg.grid.x_max - cfg.grid.x_min) / Nx;
    double dy = (cfg.grid.y_max - cfg.grid.y_min) / Ny;
    double x_dia = cfg.grid.x_diaphragm;
    double y_dia = cfg.grid.y_diaphragm;
    double gamma = cfg.phys.gamma;

    out << "x,y,rho,u,v,p,E\n";
    for (int i = 0; i <= Nx; ++i) {
        double x = cfg.grid.x_min + i * dx;
        for (int j = 0; j <= Ny; ++j) {
            double y = cfg.grid.y_min + j * dy;
            // Choose left and right states based on y
            const State& left = (y < y_dia) ? cfg.phys.left_bottom : cfg.phys.left_top;
            const State& right = (y < y_dia) ? cfg.phys.right_bottom : cfg.phys.right_top;
            double xi = (t > 1e-12) ? (x - x_dia) / t : 0.0;
            State W_ana = solve_general_riemann_problem(left, right, xi, cfg, cfg.approx_type);
            // v is already set by upwind in riemann solver, but ensure it's from correct side
            if (y < y_dia)
                W_ana.v = (W_ana.u >= 0.0) ? left.v : right.v; // actually already done, but for clarity
            double E = W_ana.p / (gamma - 1.0) + 0.5 * W_ana.rho * (W_ana.u * W_ana.u + W_ana.v * W_ana.v);
            out << x << "," << y << "," << W_ana.rho << "," << W_ana.u << "," << W_ana.v << "," << W_ana.p << "," << E << "\n";
        }
    }
    out.close();
}

void save_analytical_solution(const Config& cfg, const std::string& outputFilename) {
    generate_analytical_solution(cfg, cfg.grid.t_final, outputFilename);
}

void generate_analytical_snapshot(const Config& cfg, double t, const std::string& filename) {
    generate_analytical_solution(cfg, t, filename);
}