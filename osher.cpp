#include "euler_utils.h"
#include "choice_of_riemann_solvers.h" // для solve_rusanov_flux
#include <cmath>
#include <algorithm>

static State calc_sonic_state_1(const State& W, double gamma) {
    double c = std::sqrt(gamma * W.p / W.rho);
    double k = (gamma - 1.0) / (gamma + 1.0);
    double c_sonic = k * (W.u + 2.0 * c / (gamma - 1.0));
    State S;
    S.u = c_sonic;
    S.rho = W.rho * std::pow(c_sonic / c, 2.0 / (gamma - 1.0));
    S.p = W.p * std::pow(c_sonic / c, 2.0 * gamma / (gamma - 1.0));
    S.v = W.v;
    return S;
}

static State calc_sonic_state_3(const State& W, double gamma) {
    double c = std::sqrt(gamma * W.p / W.rho);
    double k = (gamma - 1.0) / (gamma + 1.0);
    double J_minus = W.u - 2.0 * c / (gamma - 1.0);
    double c_sonic = -J_minus * k;
    State S;
    S.u = -c_sonic;
    S.rho = W.rho * std::pow(c_sonic / c, 2.0 / (gamma - 1.0));
    S.p = W.p * std::pow(c_sonic / c, 2.0 * gamma / (gamma - 1.0));
    S.v = W.v;
    return S;
}

Flux solve_osher_flux(const State& W_L, const State& W_R, const Config& cfg, char dir) {
    const double gamma = cfg.phys.gamma;
    const double gm1 = gamma - 1.0;

    // Normalize states: normal velocity always in .u, tangential in .v
    State L = W_L, R = W_R;
    if (dir == 'y') {
        std::swap(L.u, L.v);
        std::swap(R.u, R.v);
    }

    double cL = std::sqrt(gamma * L.p / L.rho);
    double cR = std::sqrt(gamma * R.p / R.rho);

    double HL = L.u + 2.0 * cL / gm1;
    double HR = R.u - 2.0 * cR / gm1;
    double z = gm1 / (2.0 * gamma);
    double KL = std::pow(L.p, z) / cL;
    double KR = std::pow(R.p, z) / cR;
    double u_star = (KL * HL + KR * HR) / (KL + KR);
    double cL_star = cL + 0.5 * gm1 * (L.u - u_star);
    double cR_star = cR + 0.5 * gm1 * (u_star - R.u);

    State W13, W23;
    W13.u = u_star; W23.u = u_star;
    W13.rho = L.rho * std::pow(cL_star / cL, 2.0 / gm1);
    W13.p = L.p * std::pow(cL_star / cL, 2.0 * gamma / gm1);
    W13.v = L.v;
    W23.rho = R.rho * std::pow(cR_star / cR, 2.0 / gm1);
    W23.p = R.p * std::pow(cR_star / cR, 2.0 * gamma / gm1);
    W23.v = R.v;

    auto flux_of = [&](const State& W) -> Flux {
        if (dir == 'x') return fluxX(W, gamma);
        else {
            State Wy = W;
            std::swap(Wy.u, Wy.v);
            return fluxY(Wy, gamma);
        }
        };

    Flux F_total = flux_of(L);

    // Path 1: L -> 1/3
    double lambda_start = L.u - cL;
    double lambda_end = W13.u - cL_star;
    if (!(lambda_start >= 0 && lambda_end >= 0)) {
        if (lambda_start <= 0 && lambda_end <= 0) {
            F_total = F_total + (flux_of(W13) - flux_of(L));
        }
        else {
            State sonic = calc_sonic_state_1(L, gamma);
            if (lambda_start < 0)
                F_total = F_total + (flux_of(sonic) - flux_of(L));
            else
                F_total = F_total + (flux_of(W13) - flux_of(sonic));
        }
    }

    // Path 2: 1/3 -> 2/3
    if (u_star < 0) {
        F_total = F_total + (flux_of(W23) - flux_of(W13));
    }

    // Path 3: 2/3 -> R
    lambda_start = W23.u + cR_star;
    lambda_end = R.u + cR;
    if (!(lambda_start >= 0 && lambda_end >= 0)) {
        if (lambda_start <= 0 && lambda_end <= 0) {
            F_total = F_total + (flux_of(R) - flux_of(W23));
        }
        else {
            State sonic = calc_sonic_state_3(R, gamma);
            if (lambda_start < 0)
                F_total = F_total + (flux_of(sonic) - flux_of(W23));
            else
                F_total = F_total + (flux_of(R) - flux_of(sonic));
        }
    }

    if (!std::isfinite(F_total.rho_f) || !std::isfinite(F_total.rhou_f) ||
        !std::isfinite(F_total.rhov_f) || !std::isfinite(F_total.E_f))
        return solve_rusanov_flux(W_L, W_R, cfg, dir);

    return F_total;
}