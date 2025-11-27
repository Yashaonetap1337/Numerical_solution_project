#include "riemann_solver.h"
#include "euler_utils.h"
#include "choice_of_riemann_solvers.h" 
#include <stdexcept>


// акустический решатель Римана
static State solve_acoustic_riemann_problem(const State& W_L, const State& W_R, const double gamma) {
    const double a_L = soundSpeed(W_L, gamma);
    const double a_R = soundSpeed(W_R, gamma);

    // используем простые средние значения для всех параметров
    const double rho_avg = 0.5 * (W_L.rho + W_R.rho);
    const double a_avg = 0.5 * (a_L + a_R);

    // простые формулы для p* и u* (PVRS приближение)
    const double p_star = 0.5 * (W_L.p + W_R.p) - 0.5 * (W_R.u - W_L.u) * rho_avg * a_avg;
    const double u_star = 0.5 * (W_L.u + W_R.u) - 0.5 * (W_R.p - W_L.p) / (rho_avg * a_avg);

    // отбор решения на границе (x/t = 0)
    State solution_state;
    const double S_L = W_L.u - a_L; // приближенная скорость левой волны
    const double S_R = W_R.u + a_R; // приближенная скорость правой волны

    if (S_L > 0.0) {
        solution_state = W_L;
    }
    else if (S_R < 0.0) {
        solution_state = W_R;
    }
    else {

        solution_state.rho = rho_avg;
        solution_state.u = u_star;
        solution_state.p = std::max(1e-9, p_star);
    }

    return solution_state;
}

State solve_riemann_problem(const State& W_L, const State& W_R, double xi, const Config& cfg, const ApproximationType approx_type) {
    switch (cfg.riemann_solver_type) {
    case RiemannSolverType::EXACT:
        return solve_general_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
    case RiemannSolverType::ACOUSTIC:
        return solve_acoustic_riemann_problem(W_L, W_R, cfg.phys.gamma);
    default:
        throw std::runtime_error("Unknown Riemann solver type selected in config!");
    }
}