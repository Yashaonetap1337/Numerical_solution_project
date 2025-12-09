#include "types.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <iostream>


Flux solve_rusanov_flux(const State& W_L, const State& W_R, const Config& cfg) {
    const double gamma = cfg.phys.gamma;

    Flux F_L = physToFlux(W_L, gamma);
    Flux F_R = physToFlux(W_R, gamma);

    Flux F_central;
    F_central.rho_f = 0.5 * (F_L.rho_f + F_R.rho_f);
    F_central.rhou_f = 0.5 * (F_L.rhou_f + F_R.rhou_f);
    F_central.E_f = 0.5 * (F_L.E_f + F_R.E_f);

    Conserved U_L = physToCons(W_L, gamma);
    Conserved U_R = physToCons(W_R, gamma);

    double c_L = soundSpeed(W_L, gamma);
    double c_R = soundSpeed(W_R, gamma);

    double lambda_L = std::abs(W_L.u) + c_L;
    double lambda_R = std::abs(W_R.u) + c_R;

    double spectral_radius = std::max(lambda_L, lambda_R);
 
    // Обычно omega >= 1. Берем его из конфига (visc_coeff).
    double omega = (cfg.viscosity_coeff > 0.0) ? cfg.viscosity_coeff : 1.0;

    Conserved Phi;
    double alpha = omega * spectral_radius;

    Phi.rho = alpha * (U_R.rho - U_L.rho);
    Phi.rhou = alpha * (U_R.rhou - U_L.rhou);
    Phi.E = alpha * (U_R.E - U_L.E);

    Flux F_rus;
    F_rus.rho_f = F_central.rho_f - 0.5 * Phi.rho;
    F_rus.rhou_f = F_central.rhou_f - 0.5 * Phi.rhou;
    F_rus.E_f = F_central.E_f - 0.5 * Phi.E;

    return F_rus;
}