#include "weno.h"
#include "choice_of_riemann_solvers.h" 
#include "euler_utils.h"
#include <vector>
#include <cmath>

// --- Вспомогательная функция для WENO-реконструкции ОДНОЙ переменной (ЛЕВАЯ) ---
// Вычисляет значение СЛЕВА от границы i+1/2 (т.е., на правой стороне ячейки i)
static double weno5_reconstruct_left(const std::vector<double>& v, int i) {
    // Шаблоны строятся относительно ячейки i
    const double v_m2 = v[i - 2], v_m1 = v[i - 1], v_0 = v[i], v_p1 = v[i + 1], v_p2 = v[i + 2];

    // --- Шаг 1: Полиномы для 3-х шаблонов ---
    const double p0 = (2.0 / 6.0) * v_m2 - (7.0 / 6.0) * v_m1 + (11.0 / 6.0) * v_0;
    const double p1 = (-1.0 / 6.0) * v_m1 + (5.0 / 6.0) * v_0 + (2.0 / 6.0) * v_p1;
    const double p2 = (2.0 / 6.0) * v_0 + (5.0 / 6.0) * v_p1 - (1.0 / 6.0) * v_p2;

    // --- Шаг 2: Индикаторы гладкости ---
    const double beta0 = (13.0 / 12.0) * std::pow(v_m2 - 2 * v_m1 + v_0, 2) + (1.0 / 4.0) * std::pow(v_m2 - 4 * v_m1 + 3 * v_0, 2);
    const double beta1 = (13.0 / 12.0) * std::pow(v_m1 - 2 * v_0 + v_p1, 2) + (1.0 / 4.0) * std::pow(v_m1 - v_p1, 2);
    const double beta2 = (13.0 / 12.0) * std::pow(v_0 - 2 * v_p1 + v_p2, 2) + (1.0 / 4.0) * std::pow(3 * v_0 - 4 * v_p1 + v_p2, 2);

    // --- Шаг 3: Веса ---
    // **ИСПРАВЛЕНО:** Канонические идеальные веса для 5-го порядка
    const double gamma0 = 0.1, gamma1 = 0.6, gamma2 = 0.3;
    const double epsilon = 1e-6;

    const double alpha0 = gamma0 / std::pow(epsilon + beta0, 2);
    const double alpha1 = gamma1 / std::pow(epsilon + beta1, 2);
    const double alpha2 = gamma2 / std::pow(epsilon + beta2, 2);

    const double alpha_sum = alpha0 + alpha1 + alpha2;

    const double w0 = alpha0 / alpha_sum;
    const double w1 = alpha1 / alpha_sum;
    const double w2 = alpha2 / alpha_sum;

    return w0 * p0 + w1 * p1 + w2 * p2;
}

// --- Вспомогательная функция для WENO-реконструкции ОДНОЙ переменной (ПРАВАЯ) ---
// Вычисляет значение СПРАВА от границы i-1/2 (т.е., на левой стороне ячейки i)
static double weno5_reconstruct_right(const std::vector<double>& v, int i) {
    const double v_m2 = v[i - 2], v_m1 = v[i - 1], v_0 = v[i], v_p1 = v[i + 1], v_p2 = v[i + 2];

    // --- Шаг 1: Полиномы для 3-х шаблонов ---
    const double p0 = (-1.0 / 6.0) * v_m2 + (5.0 / 6.0) * v_m1 + (2.0 / 6.0) * v_0;
    const double p1 = (2.0 / 6.0) * v_m1 + (5.0 / 6.0) * v_0 - (1.0 / 6.0) * v_p1;
    const double p2 = (11.0 / 6.0) * v_0 - (7.0 / 6.0) * v_p1 + (2.0 / 6.0) * v_p2;

    // --- Шаг 2: Индикаторы гладкости (те же самые) ---
    const double beta0 = (13.0 / 12.0) * std::pow(v_m2 - 2 * v_m1 + v_0, 2) + (1.0 / 4.0) * std::pow(v_m2 - 4 * v_m1 + 3 * v_0, 2);
    const double beta1 = (13.0 / 12.0) * std::pow(v_m1 - 2 * v_0 + v_p1, 2) + (1.0 / 4.0) * std::pow(v_m1 - v_p1, 2);
    const double beta2 = (13.0 / 12.0) * std::pow(v_0 - 2 * v_p1 + v_p2, 2) + (1.0 / 4.0) * std::pow(3 * v_0 - 4 * v_p1 + v_p2, 2);

    // --- Шаг 3: Веса ---
    // **ИСПРАВЛЕНО:** Зеркальные канонические веса
    const double gamma0 = 0.3, gamma1 = 0.6, gamma2 = 0.1;
    const double epsilon = 1e-6;

    const double alpha0 = gamma0 / std::pow(epsilon + beta0, 2);
    const double alpha1 = gamma1 / std::pow(epsilon + beta1, 2);
    const double alpha2 = gamma2 / std::pow(epsilon + beta2, 2);

    const double alpha_sum = alpha0 + alpha1 + alpha2;

    const double w0 = alpha0 / alpha_sum;
    const double w1 = alpha1 / alpha_sum;
    const double w2 = alpha2 / alpha_sum;

    return w0 * p0 + w1 * p1 + w2 * p2;
}


void weno_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    std::vector<double> rho(total_cells), u(total_cells), p(total_cells);
    for (int i = 0; i < total_cells; ++i) {
        rho[i] = grid.W[i].rho; u[i] = grid.W[i].u; p[i] = grid.W[i].p;
    }

    for (int i = 0; i <= grid.Nx; ++i) {
        int cell_i = i + grid.num_fict - 1;
        int cell_i_plus_1 = i + grid.num_fict;

        // Проверка на выход за границы шаблона (остается без изменений)
        if (cell_i < 2 || cell_i > total_cells - 4) { // cell_i_plus_1 будет в границах автоматически
            const State W_L = grid.W[cell_i];
            const State W_R = grid.W[cell_i + 1];
            const State state_at_interface = solve_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
            fluxes[i] = physToFlux(state_at_interface, gamma);
            continue;
        }

        //  Реконструкция на границе i+1/2
        State W_L_interface, W_R_interface;

        // W_L_interface: Значение СЛЕВА от границы i+1/2.
        // Это экстраполяция из ячейки `cell_i` на ее правую границу.
        W_L_interface.rho = weno5_reconstruct_left(rho, cell_i);
        W_L_interface.u = weno5_reconstruct_left(u, cell_i);
        W_L_interface.p = weno5_reconstruct_left(p, cell_i);

        // W_R_interface: Значение СПРАВА от границы i+1/2.
        // Это экстраполяция из ячейки `cell_i_plus_1` на ее левую границу.
        W_R_interface.rho = weno5_reconstruct_right(rho, cell_i_plus_1);
        W_R_interface.u = weno5_reconstruct_right(u, cell_i_plus_1);
        W_R_interface.p = weno5_reconstruct_right(p, cell_i_plus_1);

        const State state_at_interface = solve_riemann_problem(W_L_interface, W_R_interface, 0.0, cfg, cfg.approx_type);
        fluxes[i] = physToFlux(state_at_interface, gamma);
    }
}
