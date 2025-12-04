#include "flux_correction.h"
#include "euler_utils.h"
#include "types.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

static double sign(double x) {
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
}

// Искусственная вязкость
static void apply_artificial_viscosity(Grid& grid, const Config& cfg, double coeff) {
    int Nx = grid.Nx;
    int num_fict = grid.num_fict;
    int total_cells = grid.Nx + 2 * grid.num_fict;


    for (int i = 0; i < total_cells; ++i) {
        grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);
    }

    std::vector<Conserved> U_smoothed = grid.U;


    for (int i = num_fict; i < Nx + num_fict; ++i) {
        double p_next = grid.W[i + 1].p;
        double p_curr = grid.W[i].p;
        double p_prev = grid.W[i - 1].p;

        double numerator = std::abs(p_next - 2.0 * p_curr + p_prev);
        double denominator = p_next + 2.0 * p_curr + p_prev + 1e-9;
        double sensor = numerator / denominator;

        // Диффузионный член: D = U_{i+1} - 2U_i + U_{i-1}
        Conserved diffusion_term = grid.U[i + 1] - grid.U[i] * 2.0 + grid.U[i - 1];

        // U_new = U_old + coeff * sensor * Diffusion
        Conserved dissipation = diffusion_term * (coeff * sensor);

        U_smoothed[i] = U_smoothed[i] + dissipation;
    }

    grid.U = U_smoothed;
}

// FCT (Диффузия и Антидиффузия)
static void apply_fct_correction(Grid& grid, double Q_coeff) {
    int Nx = grid.Nx;
    int num_fict = grid.num_fict;
    int total_cells = grid.Nx + 2 * grid.num_fict;

    // Вычисляем диффузионные потоки Phi = Q * (U_{i+1} - U_i)
    std::vector<Conserved> diff_fluxes(total_cells);
    for (int i = 0; i < total_cells - 1; ++i) {
        diff_fluxes[i] = (grid.U[i + 1] - grid.U[i]) * Q_coeff;
    }

    // Применяем диффузию: U^diff
    std::vector<Conserved> U_diff = grid.U;
    for (int i = 1; i < total_cells - 1; ++i) {
        U_diff[i] = grid.U[i] + (diff_fluxes[i] - diff_fluxes[i - 1]);
    }

    // Вычисляем антидиффузионные потоки с лимитером
    std::vector<Conserved> delta_U_diff(total_cells);
    for (int i = 0; i < total_cells - 1; ++i) {
        delta_U_diff[i] = U_diff[i + 1] - U_diff[i];
    }

    std::vector<Conserved> antidiff_fluxes(total_cells);


    for (int i = 1; i < total_cells - 2; ++i) {

        // Лямбда для обработки одной компоненты (rho, rhou или E)
        auto limit_component = [&](double phi_val, double delta_L, double delta_R) {
            double s = sign(phi_val);
            double val = std::abs(phi_val);
            double term1 = s * delta_L;
            double term2 = val;
            double term3 = s * delta_R;

            // minmod-like limiter
            double min_term = std::min({ term1, term2, term3 });
            return s * std::max(0.0, min_term);
        };

        antidiff_fluxes[i].rho = limit_component(diff_fluxes[i].rho, delta_U_diff[i - 1].rho, delta_U_diff[i + 1].rho);
        antidiff_fluxes[i].rhou = limit_component(diff_fluxes[i].rhou, delta_U_diff[i - 1].rhou, delta_U_diff[i + 1].rhou);
        antidiff_fluxes[i].E = limit_component(diff_fluxes[i].E, delta_U_diff[i - 1].E, delta_U_diff[i + 1].E);
    }

    // Вычитаем антидиффузию
    for (int i = num_fict; i < Nx + num_fict; ++i) {
        grid.U[i] = U_diff[i] - (antidiff_fluxes[i] - antidiff_fluxes[i - 1]);
    }
}


void apply_flux_correction(Grid& grid, const Config& cfg) {

    if (cfg.flux_correction == FluxCorrectionType::NONE) {
        return;
    }

    for (int i = 0; i < grid.U.size(); ++i) {
        grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);
    }

    switch (cfg.flux_correction) {
    case FluxCorrectionType::VISCOSITY:
        apply_artificial_viscosity(grid, cfg, cfg.viscosity_coeff);
        break;

    case FluxCorrectionType::FCT:
        apply_fct_correction(grid, cfg.viscosity_coeff);
        break;

    default:
        break;
    }
}