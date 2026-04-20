#include "euler_utils.h"
#include <cmath>

Flux solve_rusanov_flux(const State& W_L, const State& W_R, const Config& cfg, char dir) {
    const double gamma = cfg.phys.gamma;
    Flux FL = (dir == 'x') ? fluxX(W_L, gamma) : fluxY(W_L, gamma);
    Flux FR = (dir == 'x') ? fluxX(W_R, gamma) : fluxY(W_R, gamma);

    Flux Fc;
    Fc.rho_f = 0.5 * (FL.rho_f + FR.rho_f);
    Fc.rhou_f = 0.5 * (FL.rhou_f + FR.rhou_f);
    Fc.rhov_f = 0.5 * (FL.rhov_f + FR.rhov_f);
    Fc.E_f = 0.5 * (FL.E_f + FR.E_f);

    Conserved UL = physToCons(W_L, gamma);
    Conserved UR = physToCons(W_R, gamma);
    double aL = soundSpeed(W_L, gamma);
    double aR = soundSpeed(W_R, gamma);
    double uL_n = (dir == 'x') ? W_L.u : W_L.v;
    double uR_n = (dir == 'x') ? W_R.u : W_R.v;
    double lambdaL = std::abs(uL_n) + aL;
    double lambdaR = std::abs(uR_n) + aR;
    double spectral_radius = std::max(lambdaL, lambdaR);
    double omega = (cfg.viscosity_coeff > 0.0) ? cfg.viscosity_coeff : 1.2;
    double alpha = omega * spectral_radius;

    Conserved Phi;
    Phi.rho = alpha * (UR.rho - UL.rho);
    Phi.rhou = alpha * (UR.rhou - UL.rhou);
    Phi.rhov = alpha * (UR.rhov - UL.rhov);
    Phi.E = alpha * (UR.E - UL.E);

    Flux F;
    F.rho_f = Fc.rho_f - 0.5 * Phi.rho;
    F.rhou_f = Fc.rhou_f - 0.5 * Phi.rhou;
    F.rhov_f = Fc.rhov_f - 0.5 * Phi.rhov;
    F.E_f = Fc.E_f - 0.5 * Phi.E;
    return F;
}