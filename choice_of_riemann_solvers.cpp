#pragma once
#include "types.h"
#include "euler_utils.h"
#include "hllc.h"
#include <stdexcept>

State solve_general_riemann_problem(const State& W_L, const State& W_R, double xi, const Config& cfg, const ApproximationType approx_type);
Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg, char dir);
Flux hll_flux(const State& W_L, const State& W_R, double gamma, char dir);
Flux hllc_flux(const State& W_L, const State& W_R, double gamma, char dir, WaveSpeedMethod method);
Flux solve_rusanov_flux(const State& W_L, const State& W_R, const Config& cfg, char dir);
Flux solve_osher_flux(const State& W_L, const State& W_R, const Config& cfg, char dir);
State solve_acoustic_riemann_problem(const State& W_L, const State& W_R, double gamma);

inline Flux compute_interface_flux(const State& W_L, const State& W_R, const Config& cfg, char dir) {
    switch (cfg.riemann_solver_type) {
        case RiemannSolverType::ROE:
            return solve_roe_flux(W_L, W_R, cfg, dir);
        case RiemannSolverType::HLL:
            return hll_flux(W_L, W_R, cfg.phys.gamma, dir);
        case RiemannSolverType::HLLC:
            return hllc_flux(W_L, W_R, cfg.phys.gamma, dir);
        case RiemannSolverType::RUSANOV:
            return solve_rusanov_flux(W_L, W_R, cfg, dir);
        case RiemannSolverType::OSHER:
            return solve_osher_flux(W_L, W_R, cfg, dir);

        case RiemannSolverType::EXACT:
        {
            // ѕоворот: дл€ y-направлени€ v становитс€ нормальной скоростью
            State WL = W_L, WR = W_R;
            if (dir == 'y') {
                std::swap(WL.u, WL.v);
                std::swap(WR.u, WR.v);
            }
            State W_star = solve_general_riemann_problem(WL, WR, 0.0, cfg, cfg.approx_type);
            if (dir == 'y') {
                // ќбратный поворот
                std::swap(W_star.u, W_star.v);
            }
            return (dir == 'x') ? fluxX(W_star, cfg.phys.gamma)
                                 : fluxY(W_star, cfg.phys.gamma);
        }

        case RiemannSolverType::ACOUSTIC:
        {
            State WL = W_L, WR = W_R;
            if (dir == 'y') {
                std::swap(WL.u, WL.v);
                std::swap(WR.u, WR.v);
            }
            State W_star = solve_acoustic_riemann_problem(WL, WR, cfg.phys.gamma);
            if (dir == 'y') {
                std::swap(W_star.u, W_star.v);
            }
            return (dir == 'x') ? fluxX(W_star, cfg.phys.gamma)
                                 : fluxY(W_star, cfg.phys.gamma);
        }

        default:
            throw std::runtime_error("Unknown Riemann solver type");
    }
}