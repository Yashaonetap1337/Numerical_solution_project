#include "hll.h"
#include "euler_utils.h"
#include <cmath>

static void hll_wave_speeds(const State& W_L, const State& W_R, double gamma, char dir,
    double& S_L, double& S_R) {
    double aL = soundSpeed(W_L, gamma);
    double aR = soundSpeed(W_R, gamma);
    double uL_n = (dir == 'x') ? W_L.u : W_L.v;
    double uR_n = (dir == 'x') ? W_R.u : W_R.v;
    S_L = std::min(uL_n - aL, uR_n - aR);
    S_R = std::max(uL_n + aL, uR_n + aR);
    const double eps = 1e-8;
    if (S_L >= S_R) {
        S_L = std::min(S_L, -eps);
        S_R = std::max(S_R, eps);
    }
}

Flux hll_flux(const State& W_L, const State& W_R, double gamma, char dir) {
    double S_L, S_R;
    hll_wave_speeds(W_L, W_R, gamma, dir, S_L, S_R);
    Conserved UL = physToCons(W_L, gamma);
    Conserved UR = physToCons(W_R, gamma);
    Flux FL = (dir == 'x') ? fluxX(W_L, gamma) : fluxY(W_L, gamma);
    Flux FR = (dir == 'x') ? fluxX(W_R, gamma) : fluxY(W_R, gamma);
    Flux F;
    if (S_L >= 0.0) {
        F = FL;
    }
    else if (S_R <= 0.0) {
        F = FR;
    }
    else {
        double denom = S_R - S_L;
        if (std::abs(denom) < 1e-12) {
            F = (FL + FR) * 0.5;
        }
        else {
            double inv = 1.0 / denom;
            F.rho_f = inv * (S_R * FL.rho_f - S_L * FR.rho_f + S_L * S_R * (UR.rho - UL.rho));
            F.rhou_f = inv * (S_R * FL.rhou_f - S_L * FR.rhou_f + S_L * S_R * (UR.rhou - UL.rhou));
            F.rhov_f = inv * (S_R * FL.rhov_f - S_L * FR.rhov_f + S_L * S_R * (UR.rhov - UL.rhov));
            F.E_f = inv * (S_R * FL.E_f - S_L * FR.E_f + S_L * S_R * (UR.E - UL.E));
        }
    }
    return F;
}