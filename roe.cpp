#include "roe.h"
#include "euler_utils.h"
#include "choice_of_riemann_solvers.h"
#include <cmath>

static void roe_average(const State& W_L, const State& W_R, double gamma,
    double& rho_tilde, double& u_tilde, double& v_tilde, double& H_tilde) {
    double sqrtL = std::sqrt(W_L.rho);
    double sqrtR = std::sqrt(W_R.rho);
    double sum = sqrtL + sqrtR;
    rho_tilde = sqrtL * sqrtR;
    u_tilde = (sqrtL * W_L.u + sqrtR * W_R.u) / sum;
    v_tilde = (sqrtL * W_L.v + sqrtR * W_R.v) / sum;
    double EL = W_L.p / (gamma - 1.0) + 0.5 * W_L.rho * (W_L.u * W_L.u + W_L.v * W_L.v);
    double ER = W_R.p / (gamma - 1.0) + 0.5 * W_R.rho * (W_R.u * W_R.u + W_R.v * W_R.v);
    double HL = (EL + W_L.p) / W_L.rho;
    double HR = (ER + W_R.p) / W_R.rho;
    H_tilde = (sqrtL * HL + sqrtR * HR) / sum;
}

static double roe_sound_speed(double u_tilde, double v_tilde, double H_tilde, double gamma) {
    double q2 = u_tilde * u_tilde + v_tilde * v_tilde;
    double a2 = (gamma - 1.0) * (H_tilde - 0.5 * q2);
    return std::sqrt(std::max(1e-12, a2));
}

static double entropy_fix(double lambda, double delta) {
    if (std::abs(lambda) < delta)
        return (lambda * lambda + delta * delta) / (2.0 * delta);
    return std::abs(lambda);
}

Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg, char dir) {
    const double gamma = cfg.phys.gamma;
    State safeL = W_L, safeR = W_R;
    safeL.rho = std::max(1e-9, safeL.rho);
    safeR.rho = std::max(1e-9, safeR.rho);
    safeL.p = std::max(1e-9, safeL.p);
    safeR.p = std::max(1e-9, safeR.p);

    double uL_n = (dir == 'x') ? safeL.u : safeL.v;
    double uR_n = (dir == 'x') ? safeR.u : safeR.v;
    double uL_t = (dir == 'x') ? safeL.v : safeL.u;
    double uR_t = (dir == 'x') ? safeR.v : safeR.u;

    double aL = soundSpeed(safeL, gamma);
    double aR = soundSpeed(safeR, gamma);

    bool low_density = (safeL.rho < 1e-4) || (safeR.rho < 1e-4);
    double escape_vel = 2.0 * (aL + aR) / (gamma - 1.0);
    bool vacuum = (uR_n - uL_n) > escape_vel;
    double p_ratio = std::max(safeL.p, safeR.p) / std::min(safeL.p, safeR.p);
    bool strong_shock = p_ratio > 100.0;

    if (low_density || vacuum || strong_shock)
        return solve_rusanov_flux(safeL, safeR, cfg, dir);

    double rho_tilde, u_tilde, v_tilde, H_tilde;
    roe_average(safeL, safeR, gamma, rho_tilde, u_tilde, v_tilde, H_tilde);
    double a_tilde = roe_sound_speed(u_tilde, v_tilde, H_tilde, gamma);

    double un_tilde = (dir == 'x') ? u_tilde : v_tilde;
    double ut_tilde = (dir == 'x') ? v_tilde : u_tilde;

    double lambda[4] = { un_tilde - a_tilde, un_tilde, un_tilde, un_tilde + a_tilde };
    double delta_fix = 0.2 * a_tilde;
    double vel_diff = uR_n - uL_n;
    if (vel_diff > 0) delta_fix += 0.4 * vel_diff;
    delta_fix = std::max(delta_fix, 1e-5);

    double drho = safeR.rho - safeL.rho;
    double du_n = uR_n - uL_n;
    double du_t = (dir == 'x') ? (safeR.v - safeL.v) : (safeR.u - safeL.u);
    double dp = safeR.p - safeL.p;

    double c = a_tilde;
    double alpha1 = 0.5 * (dp - rho_tilde * c * du_n) / (c * c);
    double alpha2 = drho - dp / (c * c);
    double alpha3 = rho_tilde * du_t;
    double alpha4 = 0.5 * (dp + rho_tilde * c * du_n) / (c * c);

    // Eigenvectors for conservative variables
    Conserved r[4];
    r[0] = { 1.0, (dir == 'x') ? un_tilde - c : ut_tilde, (dir == 'x') ? ut_tilde : un_tilde - c, H_tilde - un_tilde * c };
    r[1] = { 1.0, (dir == 'x') ? un_tilde : ut_tilde,     (dir == 'x') ? ut_tilde : un_tilde,     0.5 * (un_tilde * un_tilde + ut_tilde * ut_tilde) };
    r[2] = { 0.0, (dir == 'x') ? 0.0 : 1.0,               (dir == 'x') ? 1.0 : 0.0,               ut_tilde };
    r[3] = { 1.0, (dir == 'x') ? un_tilde + c : ut_tilde, (dir == 'x') ? ut_tilde : un_tilde + c, H_tilde + un_tilde * c };

    Flux FL = (dir == 'x') ? fluxX(safeL, gamma) : fluxY(safeL, gamma);
    Flux FR = (dir == 'x') ? fluxX(safeR, gamma) : fluxY(safeR, gamma);

    Flux F;
    F.rho_f = 0.5 * (FL.rho_f + FR.rho_f);
    F.rhou_f = 0.5 * (FL.rhou_f + FR.rhou_f);
    F.rhov_f = 0.5 * (FL.rhov_f + FR.rhov_f);
    F.E_f = 0.5 * (FL.E_f + FR.E_f);

    double alpha[4] = { alpha1, alpha2, alpha3, alpha4 };
    for (int k = 0; k < 4; ++k) {
        double abs_lam = entropy_fix(lambda[k], delta_fix);
        double coeff = 0.5 * alpha[k] * abs_lam;
        F.rho_f -= coeff * r[k].rho;
        F.rhou_f -= coeff * r[k].rhou;
        F.rhov_f -= coeff * r[k].rhov;
        F.E_f -= coeff * r[k].E;
    }

    if (!std::isfinite(F.rho_f) || !std::isfinite(F.rhou_f) ||
        !std::isfinite(F.rhov_f) || !std::isfinite(F.E_f))
        return solve_rusanov_flux(safeL, safeR, cfg, dir);

    return F;
}