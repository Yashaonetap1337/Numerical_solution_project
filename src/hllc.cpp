#include "hllc.h"
#include "euler_utils.h"
#include <cmath>

static void hllc_wave_speeds(const State& W_L, const State& W_R, double gamma, char dir,
    double& S_L, double& S_R, double& S_M) {
    double aL = soundSpeed(W_L, gamma);
    double aR = soundSpeed(W_R, gamma);
    double uL_n = (dir == 'x') ? W_L.u : W_L.v;
    double uR_n = (dir == 'x') ? W_R.u : W_R.v;

    double rho_avg = 0.5 * (W_L.rho + W_R.rho);
    double a_avg = 0.5 * (aL + aR);
    double p_star = 0.5 * (W_L.p + W_R.p) - 0.5 * (uR_n - uL_n) * rho_avg * a_avg;
    p_star = std::max(1e-9, p_star);
    double u_star = 0.5 * (uL_n + uR_n) - 0.5 * (W_R.p - W_L.p) / (rho_avg * a_avg);

    auto q_func = [gamma](double p_star, double p_k, double a_k) -> double {
        if (p_star <= p_k) return 1.0;
        return std::sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / p_k - 1.0));
        };
    double qL = q_func(p_star, W_L.p, aL);
    double qR = q_func(p_star, W_R.p, aR);
    S_L = uL_n - aL * qL;
    S_R = uR_n + aR * qR;
    S_M = u_star;
}

Flux hllc_flux(const State& W_L, const State& W_R, double gamma, char dir, WaveSpeedMethod method) {
    double S_L, S_R, S_M;
    hllc_wave_speeds(W_L, W_R, gamma, dir, S_L, S_R, S_M);

    Conserved UL = physToCons(W_L, gamma);
    Conserved UR = physToCons(W_R, gamma);
    Flux FL = (dir == 'x') ? fluxX(W_L, gamma) : fluxY(W_L, gamma);
    Flux FR = (dir == 'x') ? fluxX(W_R, gamma) : fluxY(W_R, gamma);

    double uL_n = (dir == 'x') ? W_L.u : W_L.v;
    double uR_n = (dir == 'x') ? W_R.u : W_R.v;
    double uL_t = (dir == 'x') ? W_L.v : W_L.u;
    double uR_t = (dir == 'x') ? W_R.v : W_R.u;

    double p_star = W_L.p + W_L.rho * (S_L - uL_n) * (S_M - uL_n);
    double rho_starL = W_L.rho * (S_L - uL_n) / (S_L - S_M);
    double rho_starR = W_R.rho * (S_R - uR_n) / (S_R - S_M);
    double E_starL = ((S_L - uL_n) * UL.E - W_L.p * uL_n + p_star * S_M) / (S_L - S_M);
    double E_starR = ((S_R - uR_n) * UR.E - W_R.p * uR_n + p_star * S_M) / (S_R - S_M);

    Conserved U_starL, U_starR;
    if (dir == 'x') {
        U_starL = { rho_starL, rho_starL * S_M, rho_starL * uL_t, E_starL };
        U_starR = { rho_starR, rho_starR * S_M, rho_starR * uR_t, E_starR };
    }
    else {
        U_starL = { rho_starL, rho_starL * uL_t, rho_starL * S_M, E_starL };
        U_starR = { rho_starR, rho_starR * uR_t, rho_starR * S_M, E_starR };
    }

    Flux F;
    if (S_L >= 0.0) {
        F = FL;
    }
    else if (S_M >= 0.0) {
        F.rho_f = FL.rho_f + S_L * (U_starL.rho - UL.rho);
        F.rhou_f = FL.rhou_f + S_L * (U_starL.rhou - UL.rhou);
        F.rhov_f = FL.rhov_f + S_L * (U_starL.rhov - UL.rhov);
        F.E_f = FL.E_f + S_L * (U_starL.E - UL.E);
    }
    else if (S_R > 0.0) {
        F.rho_f = FR.rho_f + S_R * (U_starR.rho - UR.rho);
        F.rhou_f = FR.rhou_f + S_R * (U_starR.rhou - UR.rhou);
        F.rhov_f = FR.rhov_f + S_R * (U_starR.rhov - UR.rhov);
        F.E_f = FR.E_f + S_R * (U_starR.E - UR.E);
    }
    else {
        F = FR;
    }
    return F;
}