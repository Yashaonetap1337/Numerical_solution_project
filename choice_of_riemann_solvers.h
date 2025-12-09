#pragma once
#include "types.h"
#include "euler_utils.h" 
#include "hllc.h"

State solve_riemann_problem(const State& W_L, const State& W_R, double xi, const Config& cfg, const ApproximationType approx_type);

Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg);
Flux hll_flux(const State& W_L, const State& W_R, double gamma);
Flux hllc_flux_with_method(const State& W_L, const State& W_R, double gamma, WaveSpeedMethod method = WaveSpeedMethod::HYBRID);
Flux solve_rusanov_flux(const State& W_L, const State& W_R, const Config& cfg);
Flux solve_osher_flux(const State& W_L, const State& W_R, const Config& cfg); 


inline Flux compute_interface_flux(const State& W_L, const State& W_R, const Config& cfg) {
    if (cfg.riemann_solver_type == RiemannSolverType::ROE) {
        return solve_roe_flux(W_L, W_R, cfg);
    }
    else if (cfg.riemann_solver_type == RiemannSolverType::HLL) {
        return hll_flux(W_L, W_R, cfg.phys.gamma);
    }
    else if (cfg.riemann_solver_type == RiemannSolverType::HLLC) {
        return hllc_flux_with_method(W_L, W_R, cfg.phys.gamma);
    }
    else if (cfg.riemann_solver_type == RiemannSolverType::RUSANOV) {
        return solve_rusanov_flux(W_L, W_R, cfg);
    }
    else if (cfg.riemann_solver_type == RiemannSolverType::OSHER) {
        return solve_osher_flux(W_L, W_R, cfg);
    }
    else {
        State W_star = solve_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
        return physToFlux(W_star, cfg.phys.gamma);
    }
}
