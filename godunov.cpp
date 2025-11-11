#include "godunov.h"
#include "riemann_solver.h"
#include "euler_utils.h" 
#include <vector>


// функция, выполняющая один шаг по времени методом Годунова
void godunov_step(Grid& grid, double dt, const Config& cfg, std::vector<Flux>& fluxes) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;



    for (int i = 0; i <= grid.Nx; ++i) {
        const State W_L = grid.W[i + grid.num_fict - 1];
        const State W_R = grid.W[i + grid.num_fict];

        State state_at_interface;

        state_at_interface = solve_general_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
        /*if (i==247 || i == 248 || i == 249 || i == 250 || i == 251 || i == 252) {
            std::cout << "i: " << i << std::endl;
            std::cout << "rho:" << state_at_interface.rho << std::endl;
            std::cout << "u:" << state_at_interface.u << std::endl;
            std::cout << "----------------------------" << std::endl;
        }*/

        fluxes[i] = physToFlux(state_at_interface, gamma);
        /*if (i == 247 || i == 248 || i == 249 || i == 250 || i == 251 || i == 252) {
            std::cout << "rho_f:" << fluxes[i].rho_f << std::endl;
        }*/


    }

    for (int i = 0; i < grid.Nx; ++i) {
        const int cell_idx = i + grid.num_fict;
        const Flux& F_left = fluxes[i];
        const Flux& F_right = fluxes[i + 1];

        grid.U[cell_idx].rho -= (dt / dx) * (F_right.rho_f - F_left.rho_f);
        grid.U[cell_idx].rhou -= (dt / dx) * (F_right.rhou_f - F_left.rhou_f);
        grid.U[cell_idx].E -= (dt / dx) * (F_right.E_f - F_left.E_f);
    }
}
