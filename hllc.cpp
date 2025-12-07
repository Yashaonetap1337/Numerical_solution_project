#include "hllc.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

static void hllc_wave_speeds(const State& W_L, const State& W_R, double gamma,
    WaveSpeedMethod method, double& S_L, double& S_R, double& S_M) {

    double a_L = soundSpeed(W_L, gamma);
    double a_R = soundSpeed(W_R, gamma);

    switch (method) {
    case WaveSpeedMethod::ISOENTROPIC: {
        double u_star = 0.5 * (W_L.u + W_R.u) + (a_L - a_R) / (gamma - 1.0);
        double a_star = 0.5 * (a_L + a_R) + 0.25 * (gamma - 1.0) * (W_L.u - W_R.u);

        S_L = std::min(W_L.u - a_L, u_star - a_star);
        S_R = std::max(W_R.u + a_R, u_star + a_star);
        S_M = u_star;
        break;
    }

    case WaveSpeedMethod::LINEARIZED: {
     
        double a_bar = 0.5 * (a_L + a_R);
        double rho_bar = 0.5 * (W_L.rho + W_R.rho);

       
        double p_star = 0.5 * (W_L.p + W_R.p) - 0.5 * (W_R.u - W_L.u) * rho_bar * a_bar;
        double u_star = 0.5 * (W_L.u + W_R.u) - 0.5 * (W_R.p - W_L.p) / (rho_bar * a_bar);

        
        double a_L_star = sqrt(gamma * p_star / W_L.rho);
        double a_R_star = sqrt(gamma * p_star / W_R.rho);

        S_L = std::min(W_L.u - a_L, u_star - a_L_star);
        S_R = std::max(W_R.u + a_R, u_star + a_R_star);
        S_M = u_star;
        break;
    }

    case WaveSpeedMethod::HYBRID: {
        double a_bar = 0.5 * (a_L + a_R);
        double rho_bar = 0.5 * (W_L.rho + W_R.rho);
        double p_star = 0.5 * (W_L.p + W_R.p) - 0.5 * (W_R.u - W_L.u) * rho_bar * a_bar;
        double u_star = 0.5 * (W_L.u + W_R.u) - 0.5 * (W_R.p - W_L.p) / (rho_bar * a_bar);

        // Вычисляем H_S = p*/p_S
        double H_L = p_star / W_L.p;
        double H_R = p_star / W_R.p;

        // Вычисляем q_S
        auto q_func = [gamma](double H) -> double {
            if (H <= 1.0) {
                return 1.0;
            }
            else {
                return sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (H - 1.0));
            }
            };

        double q_L = q_func(H_L);
        double q_R = q_func(H_R);

        // Волновые скорости
        S_L = W_L.u - a_L * q_L;
        S_R = W_R.u + a_R * q_R;
        S_M = u_star;
        break;
    }

    default:
        throw std::invalid_argument("Unknown wave speed method");
    }
}


static Flux hllc_flux_with_method(const State& W_L, const State& W_R, double gamma,
    WaveSpeedMethod method = WaveSpeedMethod::HYBRID) {

    double S_L, S_R, S_M;
    hllc_wave_speeds(W_L, W_R, gamma, method, S_L, S_R, S_M);


    double rho_star_L = W_L.rho * (S_L - W_L.u) / (S_L - S_M);
    double rho_star_R = W_R.rho * (S_R - W_R.u) / (S_R - S_M);

    double p_star = 0.5 * (
        W_L.p + W_L.rho * (S_L - W_L.u) * (S_M - W_L.u) +
        W_R.p + W_R.rho * (S_R - W_R.u) * (S_M - W_R.u)
        );

    p_star = std::max(p_star, 1e-9);

    double E_L = W_L.p / (gamma - 1.0) + 0.5 * W_L.rho * W_L.u * W_L.u;
    double E_R = W_R.p / (gamma - 1.0) + 0.5 * W_R.rho * W_R.u * W_R.u;

    double E_star_L = ((S_L - W_L.u) * E_L + p_star * S_M - W_L.p * W_L.u) / (S_L - S_M);
    double E_star_R = ((S_R - W_R.u) * E_R + p_star * S_M - W_R.p * W_R.u) / (S_R - S_M);

    Conserved U_L = physToCons(W_L, gamma);
    Conserved U_R = physToCons(W_R, gamma);

    Conserved U_star_L = { rho_star_L, rho_star_L * S_M, E_star_L };
    Conserved U_star_R = { rho_star_R, rho_star_R * S_M, E_star_R };

    Flux F_L = physToFlux(W_L, gamma);
    Flux F_R = physToFlux(W_R, gamma);

    Flux flux;

    if (S_L >= 0.0) {
        flux = F_L;
    }
    else if (S_M >= 0.0) {
        flux.rho_f = F_L.rho_f + S_L * (U_star_L.rho - U_L.rho);
        flux.rhou_f = F_L.rhou_f + S_L * (U_star_L.rhou - U_L.rhou);
        flux.E_f = F_L.E_f + S_L * (U_star_L.E - U_L.E);
    }
    else if (S_R > 0.0) {
        flux.rho_f = F_R.rho_f + S_R * (U_star_R.rho - U_R.rho);
        flux.rhou_f = F_R.rhou_f + S_R * (U_star_R.rhou - U_R.rhou);
        flux.E_f = F_R.E_f + S_R * (U_star_R.E - U_R.E);
    }
    else {
        flux = F_R;
    }

    return flux;
}

void hllc_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const int num_faces = grid.Nx + 1;
    fluxes.resize(num_faces);

    const double gamma = cfg.phys.gamma;

    // По умолчанию используем гибридную оценку
    WaveSpeedMethod method = WaveSpeedMethod::HYBRID;

    for (int i = 0; i < num_faces; ++i) {
        int left_idx = i + grid.num_fict - 1;
        int right_idx = i + grid.num_fict;

        const State& W_L = grid.W[left_idx];
        const State& W_R = grid.W[right_idx];

        try {
            fluxes[i] = hllc_flux_with_method(W_L, W_R, gamma, method);
        }
        catch (const std::exception& e) {
            // Fallback на upwind
            Flux F_L = physToFlux(W_L, gamma);
            Flux F_R = physToFlux(W_R, gamma);

            if (0.5 * (W_L.u + W_R.u) >= 0.0) {
                fluxes[i] = F_L;
            }
            else {
                fluxes[i] = F_R;
            }
        }
    }
}