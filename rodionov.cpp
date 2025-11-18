#include "rodionov.h"
#include "riemann_solver.h"
#include "euler_utils.h"
#include <vector>
#include <cmath>
#include <algorithm>

static double minmod(double a, double b) {
    if (a * b <= 0) {
        return 0.0;
    }
    if (a * a < a * b) {
        return a;
    }
    return b;
}

void rodionov_step(Grid& grid, double dt, const Config& cfg, std::vector<Flux>& fluxes) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    //РЕКОНСТРУКЦИЯ
    std::vector<double> delta_rho(total_cells), delta_u(total_cells), delta_p(total_cells);

    for (int i = 1; i < total_cells - 1; ++i) {
        double drho_L = grid.W[i].rho - grid.W[i - 1].rho;
        double drho_R = grid.W[i + 1].rho - grid.W[i].rho;

        double du_L = grid.W[i].u - grid.W[i - 1].u;
        double du_R = grid.W[i + 1].u - grid.W[i].u;

        double dp_L = grid.W[i].p - grid.W[i - 1].p;
        double dp_R = grid.W[i + 1].p - grid.W[i].p;

        delta_rho[i] = minmod(drho_L, drho_R);
        delta_u[i] = minmod(du_L, du_R);
        delta_p[i] = minmod(dp_L, dp_R);
    }

    if (total_cells > 2) {
        delta_rho[0] = delta_rho[1];
        delta_rho[total_cells - 1] = delta_rho[total_cells - 2];
        delta_u[0] = delta_u[1];
        delta_u[total_cells - 1] = delta_u[total_cells - 2];
        delta_p[0] = delta_p[1];
        delta_p[total_cells - 1] = delta_p[total_cells - 2];
    }

    //2. ПРЕДИКТОР
    std::vector<Conserved> U_pred = grid.U; // копируем текущий слой

    for (int i = 1; i < total_cells - 1; ++i) {
        // Внутренние экстраполированные состояния:
        //   Q_right от (i-1): W[i-1] + 0.5 * delta_[i-1]
        //   Q_left  от (i+1): W[i+1] - 0.5 * delta_[i+1]
        State W_right = {
            grid.W[i - 1].rho + 0.5 * delta_rho[i - 1],
            grid.W[i - 1].u + 0.5 * delta_u[i - 1],
            grid.W[i - 1].p + 0.5 * delta_p[i - 1]
        };
        State W_left = {
            grid.W[i + 1].rho - 0.5 * delta_rho[i + 1],
            grid.W[i + 1].u - 0.5 * delta_u[i + 1],
            grid.W[i + 1].p - 0.5 * delta_p[i + 1]
        };


        Flux F_right = physToFlux(W_right, gamma);
        Flux F_left = physToFlux(W_left, gamma);

        U_pred[i].rho -= (dt / dx) * (F_left.rho_f - F_right.rho_f);
        U_pred[i].rhou -= (dt / dx) * (F_left.rhou_f - F_right.rhou_f);
        U_pred[i].E -= (dt / dx) * (F_left.E_f - F_right.E_f);
    }

    U_pred[0] = grid.U[0];
    U_pred[total_cells - 1] = grid.U[total_cells - 1];

    // Преобразуем U_pred -> W_pred (для корректора)
    std::vector<State> W_pred(total_cells);
    for (int i = 0; i < total_cells; ++i) {
        W_pred[i] = consToPhys(U_pred[i], gamma);
    }

    // КОРРЕКТОР

    // Промежуточное состояние: W_half = 0.5 * (W^n + W_pred)
    std::vector<State> W_half(total_cells);
    for (int i = 0; i < total_cells; ++i) {
        W_half[i].rho = 0.5 * (grid.W[i].rho + W_pred[i].rho);
        W_half[i].u = 0.5 * (grid.W[i].u + W_pred[i].u);
        W_half[i].p = 0.5 * (grid.W[i].p + W_pred[i].p);
    }

    // Вычисляем потоки на гранях (как в Kolgan, но с W_half и теми же delta_*)
    for (int i = 0; i <= grid.Nx; ++i) {
        int cell_i = i + grid.num_fict - 1;
        int cell_i_plus_1 = i + grid.num_fict;

        State W_L_interface;
        W_L_interface.rho = W_half[cell_i].rho + 0.5 * delta_rho[cell_i];
        W_L_interface.u = W_half[cell_i].u + 0.5 * delta_u[cell_i];
        W_L_interface.p = W_half[cell_i].p + 0.5 * delta_p[cell_i];

        State W_R_interface;
        W_R_interface.rho = W_half[cell_i_plus_1].rho - 0.5 * delta_rho[cell_i_plus_1];
        W_R_interface.u = W_half[cell_i_plus_1].u - 0.5 * delta_u[cell_i_plus_1];
        W_R_interface.p = W_half[cell_i_plus_1].p - 0.5 * delta_p[cell_i_plus_1];

        const State state_at_interface = solve_general_riemann_problem(
            W_L_interface, W_R_interface, 0.0, cfg, cfg.approx_type
        );

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
