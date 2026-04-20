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

// choice_of_riemann_solvers.h Ц добавить после существующих объ€влений

// ѕерегрузка дл€ 2D Ц принимает нормаль (nx, ny)
inline Flux compute_interface_flux(const State& W_L, const State& W_R,
    const Config& cfg, double nx, double ny) {
    // ѕоворот в локальную систему координат грани
    double unL = W_L.u * nx + W_L.v * ny;
    double utL = -W_L.u * ny + W_L.v * nx;
    double unR = W_R.u * nx + W_R.v * ny;
    double utR = -W_R.u * ny + W_R.v * nx;

    State W_rotL = { W_L.rho, unL, utL, W_L.p };
    State W_rotR = { W_R.rho, unR, utR, W_R.p };

    // ¬ычисл€ем 1D поток в локальной системе (ось x)
    Flux flux_rot;
    switch (cfg.riemann_solver_type) {
    case RiemannSolverType::EXACT:
        flux_rot = fluxX(solve_general_riemann_problem(W_rotL, W_rotR, 0.0, cfg, cfg.approx_type), cfg.phys.gamma);
        break;
    case RiemannSolverType::ACOUSTIC:
        flux_rot = fluxX(solve_acoustic_riemann_problem(W_rotL, W_rotR, cfg.phys.gamma), cfg.phys.gamma);
        break;
    case RiemannSolverType::ROE:
        flux_rot = solve_roe_flux(W_rotL, W_rotR, cfg, 'x');
        break;
    case RiemannSolverType::HLL:
        flux_rot = hll_flux(W_rotL, W_rotR, cfg.phys.gamma, 'x');
        break;
    case RiemannSolverType::HLLC:
        flux_rot = hllc_flux(W_rotL, W_rotR, cfg.phys.gamma, 'x');
        break;
    case RiemannSolverType::RUSANOV:
        flux_rot = solve_rusanov_flux(W_rotL, W_rotR, cfg, 'x');
        break;
    case RiemannSolverType::OSHER:
        flux_rot = solve_osher_flux(W_rotL, W_rotR, cfg, 'x');
        break;
    default:
        throw std::runtime_error("Unknown Riemann solver in 2D flux");
    }

    // ќбратный поворот потока в глобальную систему
    Flux flux;
    flux.rho_f = flux_rot.rho_f;
    flux.rhou_f = flux_rot.rhou_f * nx - flux_rot.rhov_f * ny;
    flux.rhov_f = flux_rot.rhou_f * ny + flux_rot.rhov_f * nx;
    flux.E_f = flux_rot.E_f;

    return flux;
}