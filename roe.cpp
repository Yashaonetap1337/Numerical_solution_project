#include "roe.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

// Улучшенная функция Entropy Fix (Harten)
static double entropy_fix(double lambda, double delta) {
    if (std::abs(lambda) < delta) {
        return (lambda * lambda + delta * delta) / (2.0 * delta);
    }
    return std::abs(lambda);
}

Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg) {
    const double gamma = cfg.phys.gamma;


    double rho_L = std::max(1e-9, W_L.rho);
    double rho_R = std::max(1e-9, W_R.rho);
    double p_L = std::max(1e-9, W_L.p);
    double p_R = std::max(1e-9, W_R.p);


    Conserved U_L;
    U_L.rho = rho_L;
    U_L.rhou = rho_L * W_L.u;
    U_L.E = p_L / (gamma - 1.0) + 0.5 * rho_L * W_L.u * W_L.u;

    Conserved U_R;
    U_R.rho = rho_R;
    U_R.rhou = rho_R * W_R.u;
    U_R.E = p_R / (gamma - 1.0) + 0.5 * rho_R * W_R.u * W_R.u;

    double H_L = (U_L.E + p_L) / rho_L;
    double H_R = (U_R.E + p_R) / rho_R;

    // Roe-средние
    double sqrt_rho_L = std::sqrt(rho_L);
    double sqrt_rho_R = std::sqrt(rho_R);
    double s = sqrt_rho_L + sqrt_rho_R;

    double u_tilde = (sqrt_rho_L * W_L.u + sqrt_rho_R * W_R.u) / s;
    double H_tilde = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / s;

    double a_tilde_sq = (gamma - 1.0) * (H_tilde - 0.5 * u_tilde * u_tilde);
    if (a_tilde_sq < 1e-9) a_tilde_sq = 1e-9;
    double a_tilde = std::sqrt(a_tilde_sq);

    // Собственные значения
    double lambda[3];
    lambda[0] = u_tilde - a_tilde;
    lambda[1] = u_tilde;
    lambda[2] = u_tilde + a_tilde;


    double delta_base = 0.2 * a_tilde;
    double velocity_jump = W_R.u - W_L.u;


    if (velocity_jump > 0) {
        delta_base += 0.2 * velocity_jump;
    }


    double delta = std::max(1e-5, delta_base);


    double delta_u = W_R.u - W_L.u;
    double delta_p = p_R - p_L;
    double delta_rho = rho_R - rho_L;


    double rho_tilde = sqrt_rho_L * sqrt_rho_R;
    double denom = 2.0 * a_tilde_sq;

    double alpha[3];
    alpha[0] = (delta_p - rho_tilde * a_tilde * delta_u) / denom;
    alpha[1] = delta_rho - delta_p / a_tilde_sq;
    alpha[2] = (delta_p + rho_tilde * a_tilde * delta_u) / denom;


    Conserved r[3];
    r[0] = { 1.0, u_tilde - a_tilde, H_tilde - u_tilde * a_tilde };
    r[1] = { 1.0, u_tilde, 0.5 * u_tilde * u_tilde };
    r[2] = { 1.0, u_tilde + a_tilde, H_tilde + u_tilde * a_tilde };


    Flux F_L_flux, F_R_flux;

    F_L_flux.rho_f = U_L.rhou;
    F_L_flux.rhou_f = U_L.rhou * W_L.u + p_L;
    F_L_flux.E_f = (U_L.E + p_L) * W_L.u;

    F_R_flux.rho_f = U_R.rhou;
    F_R_flux.rhou_f = U_R.rhou * W_R.u + p_R;
    F_R_flux.E_f = (U_R.E + p_R) * W_R.u;

    Conserved dissipation = { 0.0, 0.0, 0.0 };
    for (int k = 0; k < 3; ++k) {
        double abs_lambda = entropy_fix(lambda[k], delta);
        dissipation = dissipation + (r[k] * (abs_lambda * alpha[k]));
    }

    Flux F_roe;
    F_roe.rho_f = 0.5 * (F_L_flux.rho_f + F_R_flux.rho_f) - 0.5 * dissipation.rho;
    F_roe.rhou_f = 0.5 * (F_L_flux.rhou_f + F_R_flux.rhou_f) - 0.5 * dissipation.rhou;
    F_roe.E_f = 0.5 * (F_L_flux.E_f + F_R_flux.E_f) - 0.5 * dissipation.E;

    return F_roe;
}