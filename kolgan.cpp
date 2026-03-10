#include "kolgan.h"
#include "choice_of_riemann_solvers.h"
#include "slope_limiter.h"
#include "euler_utils.h"
#include <vector>
#include <cmath>

void kolgan_flux_computation(const Grid& grid, const Config& cfg,
    std::vector<std::vector<Flux>>& fluxes_x,
    std::vector<std::vector<Flux>>& fluxes_y) {
    const double gamma = cfg.phys.gamma;
    const int ng = grid.num_fict;
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const double EPS = 1e-9;

    // Fluxes in x-direction
    fluxes_x.assign(Nx + 1, std::vector<Flux>(Ny + 2 * ng));
    for (int j = 0; j < Ny + 2 * ng; ++j) {
        std::vector<double> rho(Nx + 2 * ng), u(Nx + 2 * ng), v(Nx + 2 * ng), p(Nx + 2 * ng);
        for (int i = 0; i < Nx + 2 * ng; ++i) {
            rho[i] = grid.W[i][j].rho;
            u[i] = grid.W[i][j].u;
            v[i] = grid.W[i][j].v;
            p[i] = grid.W[i][j].p;
        }
        for (int i = 0; i <= Nx; ++i) {
            int idx = i + ng;
            auto rec_rhoL = reconstruct_slope_1d(rho, idx - 1, cfg.slope_limiter);
            auto rec_rhoR = reconstruct_slope_1d(rho, idx, cfg.slope_limiter);
            auto rec_uL = reconstruct_slope_1d(u, idx - 1, cfg.slope_limiter);
            auto rec_uR = reconstruct_slope_1d(u, idx, cfg.slope_limiter);
            auto rec_vL = reconstruct_slope_1d(v, idx - 1, cfg.slope_limiter);
            auto rec_vR = reconstruct_slope_1d(v, idx, cfg.slope_limiter);
            auto rec_pL = reconstruct_slope_1d(p, idx - 1, cfg.slope_limiter);
            auto rec_pR = reconstruct_slope_1d(p, idx, cfg.slope_limiter);

            State WL, WR;
            WL.rho = rec_rhoL.second; WR.rho = rec_rhoR.first;
            WL.u = rec_uL.second;   WR.u = rec_uR.first;
            WL.v = rec_vL.second;   WR.v = rec_vR.first;
            WL.p = rec_pL.second;   WR.p = rec_pR.first;
            WL.rho = std::max(EPS, WL.rho); WR.rho = std::max(EPS, WR.rho);
            WL.p = std::max(EPS, WL.p);   WR.p = std::max(EPS, WR.p);
            fluxes_x[i][j] = compute_interface_flux(WL, WR, cfg, 'x');
        }
    }

    // Fluxes in y-direction
    fluxes_y.assign(Nx + 2 * ng, std::vector<Flux>(Ny + 1));
    for (int i = 0; i < Nx + 2 * ng; ++i) {
        std::vector<double> rho(Ny + 2 * ng), u(Ny + 2 * ng), v(Ny + 2 * ng), p(Ny + 2 * ng);
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            rho[j] = grid.W[i][j].rho;
            u[j] = grid.W[i][j].u;
            v[j] = grid.W[i][j].v;
            p[j] = grid.W[i][j].p;
        }
        for (int j = 0; j <= Ny; ++j) {
            int idx = j + ng;
            auto rec_rhoB = reconstruct_slope_1d(rho, idx - 1, cfg.slope_limiter);
            auto rec_rhoT = reconstruct_slope_1d(rho, idx, cfg.slope_limiter);
            auto rec_uB = reconstruct_slope_1d(u, idx - 1, cfg.slope_limiter);
            auto rec_uT = reconstruct_slope_1d(u, idx, cfg.slope_limiter);
            auto rec_vB = reconstruct_slope_1d(v, idx - 1, cfg.slope_limiter);
            auto rec_vT = reconstruct_slope_1d(v, idx, cfg.slope_limiter);
            auto rec_pB = reconstruct_slope_1d(p, idx - 1, cfg.slope_limiter);
            auto rec_pT = reconstruct_slope_1d(p, idx, cfg.slope_limiter);

            State WB, WT;
            WB.rho = rec_rhoB.second; WT.rho = rec_rhoT.first;
            WB.u = rec_uB.second;   WT.u = rec_uT.first;
            WB.v = rec_vB.second;   WT.v = rec_vT.first;
            WB.p = rec_pB.second;   WT.p = rec_pT.first;
            WB.rho = std::max(EPS, WB.rho); WT.rho = std::max(EPS, WT.rho);
            WB.p = std::max(EPS, WB.p);   WT.p = std::max(EPS, WT.p);
            fluxes_y[i][j] = compute_interface_flux(WB, WT, cfg, 'y');
        }
    }
}