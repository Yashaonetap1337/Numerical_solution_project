#include "acoustic.h"
#include "euler_utils.h"
#include "choice_of_riemann_solvers.h" 
#include <vector>
#include <cmath>
#include <algorithm>


void acoustic_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;

    for (int i = 0; i <= grid.Nx; ++i) {
        const State W_L = grid.W[i + grid.num_fict - 1];
        const State W_R = grid.W[i + grid.num_fict];

        const State state_at_interface = solve_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
        fluxes[i] = physToFlux(state_at_interface, gamma);
    }
}