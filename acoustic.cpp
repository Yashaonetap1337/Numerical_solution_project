#include "acoustic.h"
#include "euler_utils.h"
#include <vector>
#include <cmath>
#include <algorithm>

// акустический решатель Римана
static State solve_acoustic_riemann_problem(const State& W_L, const State& W_R, const double gamma) {
    const double a_L = soundSpeed(W_L, gamma);
    const double a_R = soundSpeed(W_R, gamma);

    // cредние значения
    const double rho_avg = 0.5 * (W_L.rho + W_R.rho);
    const double a_avg = 0.5 * (a_L + a_R);

    // давление и скорость в "звездной" области по акустическим формулам
    const double p_star = 0.5 * (W_L.p + W_R.p) - 0.5 * (W_R.u - W_L.u) * rho_avg * a_avg;
    const double u_star = 0.5 * (W_L.u + W_R.u) - 0.5 * (W_R.p - W_L.p) / (rho_avg * a_avg);

    // отбор решения на границе (x/t = 0)
    State solution_state;
    const double S_L = W_L.u - a_L; // приближенная скорость левой волны
    const double S_R = W_R.u + a_R; // приближенная скорость правой волны

    if (S_L > 0.0) {
        solution_state = W_L;
    }
    else if (S_R < 0.0) {
        solution_state = W_R;
    }
    else {
        solution_state.rho = rho_avg;
        solution_state.u = u_star;
        solution_state.p = std::max(1e-9, p_star);
    }

    return solution_state;
}



void acoustic_step(Grid& grid, double dt, const Config& cfg) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;

    std::vector<Flux> fluxes(grid.Nx + 1);

    for (int i = 0; i <= grid.Nx; ++i) {
        const State W_L = grid.W[i + grid.num_fict - 1];
        const State W_R = grid.W[i + grid.num_fict];

        const State state_at_interface = solve_acoustic_riemann_problem(W_L, W_R, gamma);
        fluxes[i] = physToFlux(state_at_interface, gamma);
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
