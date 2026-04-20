#include "acoustic.h"
#include "choice_of_riemann_solvers.h"
#include "euler_utils.h"

void acoustic_flux_computation(const Grid& grid, const Config& cfg,
    std::vector<std::vector<Flux>>& fluxes_x,
    std::vector<std::vector<Flux>>& fluxes_y) {
    const double gamma = cfg.phys.gamma;
    const int ng = grid.num_fict;
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const double EPS = 1e-9;

    fluxes_x.assign(Nx + 1, std::vector<Flux>(Ny + 2 * ng));
    for (int i = 0; i <= Nx; ++i) {
        int iL = i + ng - 1;
        int iR = i + ng;
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            State WL = grid.W[iL][j];
            State WR = grid.W[iR][j];
            if (WL.rho < EPS) WL.rho = EPS;
            if (WR.rho < EPS) WR.rho = EPS;
            if (WL.p < EPS) WL.p = EPS;
            if (WR.p < EPS) WR.p = EPS;
            fluxes_x[i][j] = compute_interface_flux(WL, WR, cfg, 'x');
        }
    }

    fluxes_y.assign(Nx + 2 * ng, std::vector<Flux>(Ny + 1));
    for (int j = 0; j <= Ny; ++j) {
        int jB = j + ng - 1;
        int jT = j + ng;
        for (int i = 0; i < Nx + 2 * ng; ++i) {
            State WB = grid.W[i][jB];
            State WT = grid.W[i][jT];
            if (WB.rho < EPS) WB.rho = EPS;
            if (WT.rho < EPS) WT.rho = EPS;
            if (WB.p < EPS) WB.p = EPS;
            if (WT.p < EPS) WT.p = EPS;
            fluxes_y[i][j] = compute_interface_flux(WB, WT, cfg, 'y');
        }
    }
}