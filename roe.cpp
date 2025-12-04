#include "roe.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <vector>


static double entropy_fix(double lambda, double delta) {
    if (std::abs(lambda) < delta) {
        return (lambda * lambda + delta * delta) / (2.0 * delta);
    }
    return std::abs(lambda);
}


Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg) {
    const double gamma = cfg.phys.gamma;
    const double delta_fix = 0.2;
    // Консервативные переменные и энтальпия
    Conserved U_L = physToCons(W_L, gamma);
    Conserved U_R = physToCons(W_R, gamma);

    double H_L = (U_L.E + W_L.p) / W_L.rho;
    double H_R = (U_R.E + W_R.p) / W_R.rho;

    // Roe-средние
    double sqrt_rho_L = std::sqrt(W_L.rho);
    double sqrt_rho_R = std::sqrt(W_R.rho);
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

    // Разности
    double delta_u = W_R.u - W_L.u;
    double delta_p = W_R.p - W_L.p;
    double delta_rho = W_R.rho - W_L.rho;

    // Коэффициенты alpha
    double rho_tilde = sqrt_rho_L * sqrt_rho_R;
    double alpha[3];
    double denom = 2.0 * a_tilde_sq;

    alpha[0] = (delta_p - rho_tilde * a_tilde * delta_u) / denom;
    alpha[1] = delta_rho - delta_p / a_tilde_sq;
    alpha[2] = (delta_p + rho_tilde * a_tilde * delta_u) / denom;

    // Собственные векторы
    Conserved r[3];
    r[0] = { 1.0, u_tilde - a_tilde, H_tilde - u_tilde * a_tilde };
    r[1] = { 1.0, u_tilde, 0.5 * u_tilde * u_tilde };
    r[2] = { 1.0, u_tilde + a_tilde, H_tilde + u_tilde * a_tilde };

    // Сборка потока
    Flux F_L = physToFlux(W_L, gamma);
    Flux F_R = physToFlux(W_R, gamma);

    Conserved dissipation = { 0.0, 0.0, 0.0 };
    for (int k = 0; k < 3; ++k) {
        double abs_lambda = entropy_fix(lambda[k], delta_fix * a_tilde);
        dissipation = dissipation + (r[k] * (abs_lambda * alpha[k]));
    }

    Flux F_roe;
    F_roe.rho_f = 0.5 * (F_L.rho_f + F_R.rho_f) - 0.5 * dissipation.rho;
    F_roe.rhou_f = 0.5 * (F_L.rhou_f + F_R.rhou_f) - 0.5 * dissipation.rhou;
    F_roe.E_f = 0.5 * (F_L.E_f + F_R.E_f) - 0.5 * dissipation.E;

    return F_roe;
}
