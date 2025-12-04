#pragma once
#include "types.h"
#include "euler_utils.h" 

State solve_riemann_problem(const State& W_L, const State& W_R, double xi, const Config& cfg, const ApproximationType approx_type);

Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg);


inline Flux compute_interface_flux(const State& W_L, const State& W_R, const Config& cfg) {
    if (cfg.riemann_solver_type == RiemannSolverType::ROE) {
        return solve_roe_flux(W_L, W_R, cfg);
    }
    else {
        State W_star = solve_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
        return physToFlux(W_star, cfg.phys.gamma);
    }
}
