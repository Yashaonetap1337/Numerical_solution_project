#include "godunov.h"
#include "choice_of_riemann_solvers.h" 
#include "euler_utils.h" 
#include <vector>


// функция, выполняющая один шаг по времени методом Годунова
void godunov_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;

    for (int i = 0; i <= grid.Nx; ++i) {
        const State W_L = grid.W[i + grid.num_fict - 1];
        const State W_R = grid.W[i + grid.num_fict];

        fluxes[i] = compute_interface_flux(W_L, W_R, cfg);
    }
}