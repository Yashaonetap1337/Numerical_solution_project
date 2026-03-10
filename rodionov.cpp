#include "rodionov.h"
#include "choice_of_riemann_solvers.h"
#include "euler_utils.h"
#include "slope_limiter.h"
#include "limiters_math.h"
#include <vector>
#include <cmath>

void rodionov_step(Grid& grid, double dt, const Config& cfg) {
    const double dx = grid.dx, dy = grid.dy;
    const double gamma = cfg.phys.gamma;
    const int ng = grid.num_fict;
    const int Nx = grid.Nx, Ny = grid.Ny;

    // Save current state
    auto U0 = grid.U;

    // Arrays for slopes (simplified: just use math_limiter)
    std::vector<std::vector<double>> rho(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    std::vector<std::vector<double>> u(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    std::vector<std::vector<double>> v(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    std::vector<std::vector<double>> p(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));

    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            rho[i][j] = grid.W[i][j].rho;
            u[i][j] = grid.W[i][j].u;
            v[i][j] = grid.W[i][j].v;
            p[i][j] = grid.W[i][j].p;
        }

    // Slopes in x direction
    std::vector<std::vector<double>> drho_x(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> du_x(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dv_x(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dp_x(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));

    for (int j = 0; j < Ny + 2 * ng; ++j)
        for (int i = 1; i < Nx + 2 * ng - 1; ++i) {
            drho_x[i][j] = math_limiter(rho[i][j] - rho[i - 1][j], rho[i + 1][j] - rho[i][j], cfg.slope_limiter);
            du_x[i][j] = math_limiter(u[i][j] - u[i - 1][j], u[i + 1][j] - u[i][j], cfg.slope_limiter);
            dv_x[i][j] = math_limiter(v[i][j] - v[i - 1][j], v[i + 1][j] - v[i][j], cfg.slope_limiter);
            dp_x[i][j] = math_limiter(p[i][j] - p[i - 1][j], p[i + 1][j] - p[i][j], cfg.slope_limiter);
        }

    // Slopes in y direction
    std::vector<std::vector<double>> drho_y(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> du_y(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dv_y(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dp_y(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));

    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 1; j < Ny + 2 * ng - 1; ++j) {
            drho_y[i][j] = math_limiter(rho[i][j] - rho[i][j - 1], rho[i][j + 1] - rho[i][j], cfg.slope_limiter);
            du_y[i][j] = math_limiter(u[i][j] - u[i][j - 1], u[i][j + 1] - u[i][j], cfg.slope_limiter);
            dv_y[i][j] = math_limiter(v[i][j] - v[i][j - 1], v[i][j + 1] - v[i][j], cfg.slope_limiter);
            dp_y[i][j] = math_limiter(p[i][j] - p[i][j - 1], p[i][j + 1] - p[i][j], cfg.slope_limiter);
        }

    // Predictor fluxes
    std::vector<std::vector<Flux>> Fx_pred(Nx + 1, std::vector<Flux>(Ny + 2 * ng));
    std::vector<std::vector<Flux>> Fy_pred(Nx + 2 * ng, std::vector<Flux>(Ny + 1));

    for (int i = 0; i <= Nx; ++i) {
        int idx = i + ng;
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            double rhoL = rho[idx - 1][j] + 0.5 * drho_x[idx - 1][j];
            double uL = u[idx - 1][j] + 0.5 * du_x[idx - 1][j];
            double vL = v[idx - 1][j] + 0.5 * dv_x[idx - 1][j];
            double pL = p[idx - 1][j] + 0.5 * dp_x[idx - 1][j];
            double rhoR = rho[idx][j] - 0.5 * drho_x[idx][j];
            double uR = u[idx][j] - 0.5 * du_x[idx][j];
            double vR = v[idx][j] - 0.5 * dv_x[idx][j];
            double pR = p[idx][j] - 0.5 * dp_x[idx][j];
            State WL = { std::max(1e-9,rhoL), uL, vL, std::max(1e-9,pL) };
            State WR = { std::max(1e-9,rhoR), uR, vR, std::max(1e-9,pR) };
            Fx_pred[i][j] = compute_interface_flux(WL, WR, cfg, 'x');
        }
    }

    for (int j = 0; j <= Ny; ++j) {
        int idx = j + ng;
        for (int i = 0; i < Nx + 2 * ng; ++i) {
            double rhoB = rho[i][idx - 1] + 0.5 * drho_y[i][idx - 1];
            double uB = u[i][idx - 1] + 0.5 * du_y[i][idx - 1];
            double vB = v[i][idx - 1] + 0.5 * dv_y[i][idx - 1];
            double pB = p[i][idx - 1] + 0.5 * dp_y[i][idx - 1];
            double rhoT = rho[i][idx] - 0.5 * drho_y[i][idx];
            double uT = u[i][idx] - 0.5 * du_y[i][idx];
            double vT = v[i][idx] - 0.5 * dv_y[i][idx];
            double pT = p[i][idx] - 0.5 * dp_y[i][idx];
            State WB = { std::max(1e-9,rhoB), uB, vB, std::max(1e-9,pB) };
            State WT = { std::max(1e-9,rhoT), uT, vT, std::max(1e-9,pT) };
            Fy_pred[i][j] = compute_interface_flux(WB, WT, cfg, 'y');
        }
    }

    // Predictor step
    std::vector<std::vector<Conserved>> U_pred = U0;
    for (int i = ng; i < Nx + ng; ++i) {
        for (int j = ng; j < Ny + ng; ++j) {
            double div_x_rho = (Fx_pred[i - ng + 1][j].rho_f - Fx_pred[i - ng][j].rho_f) / dx;
            double div_y_rho = (Fy_pred[i][j - ng + 1].rho_f - Fy_pred[i][j - ng].rho_f) / dy;
            U_pred[i][j].rho = U0[i][j].rho - 0.5 * dt * (div_x_rho + div_y_rho);
            double div_x_rhou = (Fx_pred[i - ng + 1][j].rhou_f - Fx_pred[i - ng][j].rhou_f) / dx;
            double div_y_rhou = (Fy_pred[i][j - ng + 1].rhou_f - Fy_pred[i][j - ng].rhou_f) / dy;
            U_pred[i][j].rhou = U0[i][j].rhou - 0.5 * dt * (div_x_rhou + div_y_rhou);
            double div_x_rhov = (Fx_pred[i - ng + 1][j].rhov_f - Fx_pred[i - ng][j].rhov_f) / dx;
            double div_y_rhov = (Fy_pred[i][j - ng + 1].rhov_f - Fy_pred[i][j - ng].rhov_f) / dy;
            U_pred[i][j].rhov = U0[i][j].rhov - 0.5 * dt * (div_x_rhov + div_y_rhov);
            double div_x_E = (Fx_pred[i - ng + 1][j].E_f - Fx_pred[i - ng][j].E_f) / dx;
            double div_y_E = (Fy_pred[i][j - ng + 1].E_f - Fy_pred[i][j - ng].E_f) / dy;
            U_pred[i][j].E = U0[i][j].E - 0.5 * dt * (div_x_E + div_y_E);
        }
    }

    // Apply boundary conditions to predictor (copy from U0 for simplicity)
    for (int i = 0; i < ng; ++i) {
        for (int j = 0; j < Ny + 2 * ng; ++j) U_pred[i][j] = U0[i][j];
        for (int j = 0; j < Ny + 2 * ng; ++j) U_pred[Nx + 2 * ng - 1 - i][j] = U0[Nx + 2 * ng - 1 - i][j];
    }
    for (int j = 0; j < ng; ++j) {
        for (int i = 0; i < Nx + 2 * ng; ++i) U_pred[i][j] = U0[i][j];
        for (int i = 0; i < Nx + 2 * ng; ++i) U_pred[i][Ny + 2 * ng - 1 - j] = U0[i][Ny + 2 * ng - 1 - j];
    }

    // Convert predictor to primitive
    std::vector<std::vector<State>> W_pred(Nx + 2 * ng, std::vector<State>(Ny + 2 * ng));
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            W_pred[i][j] = consToPhys(U_pred[i][j], gamma);

    // Extract arrays for predictor
    std::vector<std::vector<double>> rho_pred(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    std::vector<std::vector<double>> u_pred(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    std::vector<std::vector<double>> v_pred(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    std::vector<std::vector<double>> p_pred(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng));
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            rho_pred[i][j] = W_pred[i][j].rho;
            u_pred[i][j] = W_pred[i][j].u;
            v_pred[i][j] = W_pred[i][j].v;
            p_pred[i][j] = W_pred[i][j].p;
        }

    // Slopes for predictor
    std::vector<std::vector<double>> drho_xp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> du_xp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dv_xp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dp_xp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    for (int j = 0; j < Ny + 2 * ng; ++j)
        for (int i = 1; i < Nx + 2 * ng - 1; ++i) {
            drho_xp[i][j] = math_limiter(rho_pred[i][j] - rho_pred[i - 1][j], rho_pred[i + 1][j] - rho_pred[i][j], cfg.slope_limiter);
            du_xp[i][j] = math_limiter(u_pred[i][j] - u_pred[i - 1][j], u_pred[i + 1][j] - u_pred[i][j], cfg.slope_limiter);
            dv_xp[i][j] = math_limiter(v_pred[i][j] - v_pred[i - 1][j], v_pred[i + 1][j] - v_pred[i][j], cfg.slope_limiter);
            dp_xp[i][j] = math_limiter(p_pred[i][j] - p_pred[i - 1][j], p_pred[i + 1][j] - p_pred[i][j], cfg.slope_limiter);
        }

    std::vector<std::vector<double>> drho_yp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> du_yp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dv_yp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    std::vector<std::vector<double>> dp_yp(Nx + 2 * ng, std::vector<double>(Ny + 2 * ng, 0));
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 1; j < Ny + 2 * ng - 1; ++j) {
            drho_yp[i][j] = math_limiter(rho_pred[i][j] - rho_pred[i][j - 1], rho_pred[i][j + 1] - rho_pred[i][j], cfg.slope_limiter);
            du_yp[i][j] = math_limiter(u_pred[i][j] - u_pred[i][j - 1], u_pred[i][j + 1] - u_pred[i][j], cfg.slope_limiter);
            dv_yp[i][j] = math_limiter(v_pred[i][j] - v_pred[i][j - 1], v_pred[i][j + 1] - v_pred[i][j], cfg.slope_limiter);
            dp_yp[i][j] = math_limiter(p_pred[i][j] - p_pred[i][j - 1], p_pred[i][j + 1] - p_pred[i][j], cfg.slope_limiter);
        }

    // Corrector fluxes
    std::vector<std::vector<Flux>> Fx_corr(Nx + 1, std::vector<Flux>(Ny + 2 * ng));
    std::vector<std::vector<Flux>> Fy_corr(Nx + 2 * ng, std::vector<Flux>(Ny + 1));

    for (int i = 0; i <= Nx; ++i) {
        int idx = i + ng;
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            double rhoL = rho_pred[idx - 1][j] + 0.5 * drho_xp[idx - 1][j];
            double uL = u_pred[idx - 1][j] + 0.5 * du_xp[idx - 1][j];
            double vL = v_pred[idx - 1][j] + 0.5 * dv_xp[idx - 1][j];
            double pL = p_pred[idx - 1][j] + 0.5 * dp_xp[idx - 1][j];
            double rhoR = rho_pred[idx][j] - 0.5 * drho_xp[idx][j];
            double uR = u_pred[idx][j] - 0.5 * du_xp[idx][j];
            double vR = v_pred[idx][j] - 0.5 * dv_xp[idx][j];
            double pR = p_pred[idx][j] - 0.5 * dp_xp[idx][j];
            State WL = { std::max(1e-9,rhoL), uL, vL, std::max(1e-9,pL) };
            State WR = { std::max(1e-9,rhoR), uR, vR, std::max(1e-9,pR) };
            Fx_corr[i][j] = compute_interface_flux(WL, WR, cfg, 'x');
        }
    }

    for (int j = 0; j <= Ny; ++j) {
        int idx = j + ng;
        for (int i = 0; i < Nx + 2 * ng; ++i) {
            double rhoB = rho_pred[i][idx - 1] + 0.5 * drho_yp[i][idx - 1];
            double uB = u_pred[i][idx - 1] + 0.5 * du_yp[i][idx - 1];
            double vB = v_pred[i][idx - 1] + 0.5 * dv_yp[i][idx - 1];
            double pB = p_pred[i][idx - 1] + 0.5 * dp_yp[i][idx - 1];
            double rhoT = rho_pred[i][idx] - 0.5 * drho_yp[i][idx];
            double uT = u_pred[i][idx] - 0.5 * du_yp[i][idx];
            double vT = v_pred[i][idx] - 0.5 * dv_yp[i][idx];
            double pT = p_pred[i][idx] - 0.5 * dp_yp[i][idx];
            State WB = { std::max(1e-9,rhoB), uB, vB, std::max(1e-9,pB) };
            State WT = { std::max(1e-9,rhoT), uT, vT, std::max(1e-9,pT) };
            Fy_corr[i][j] = compute_interface_flux(WB, WT, cfg, 'y');
        }
    }

    // Corrector step
    for (int i = ng; i < Nx + ng; ++i) {
        for (int j = ng; j < Ny + ng; ++j) {
            double div_x_rho = (Fx_corr[i - ng + 1][j].rho_f - Fx_corr[i - ng][j].rho_f) / dx;
            double div_y_rho = (Fy_corr[i][j - ng + 1].rho_f - Fy_corr[i][j - ng].rho_f) / dy;
            grid.U[i][j].rho = U0[i][j].rho - dt * (div_x_rho + div_y_rho);
            double div_x_rhou = (Fx_corr[i - ng + 1][j].rhou_f - Fx_corr[i - ng][j].rhou_f) / dx;
            double div_y_rhou = (Fy_corr[i][j - ng + 1].rhou_f - Fy_corr[i][j - ng].rhou_f) / dy;
            grid.U[i][j].rhou = U0[i][j].rhou - dt * (div_x_rhou + div_y_rhou);
            double div_x_rhov = (Fx_corr[i - ng + 1][j].rhov_f - Fx_corr[i - ng][j].rhov_f) / dx;
            double div_y_rhov = (Fy_corr[i][j - ng + 1].rhov_f - Fy_corr[i][j - ng].rhov_f) / dy;
            grid.U[i][j].rhov = U0[i][j].rhov - dt * (div_x_rhov + div_y_rhov);
            double div_x_E = (Fx_corr[i - ng + 1][j].E_f - Fx_corr[i - ng][j].E_f) / dx;
            double div_y_E = (Fy_corr[i][j - ng + 1].E_f - Fy_corr[i][j - ng].E_f) / dy;
            grid.U[i][j].E = U0[i][j].E - dt * (div_x_E + div_y_E);
        }
    }
}