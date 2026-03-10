#include "maccormack.h"
#include "euler_utils.h"
#include <vector>

void maccormack_step(Grid& grid, double dt, const Config& cfg) {
    const double dx = grid.dx, dy = grid.dy;
    const double gamma = cfg.phys.gamma;
    const int ng = grid.num_fict;
    const int Nx = grid.Nx, Ny = grid.Ny;
    auto U0 = grid.U;

    std::vector<std::vector<Flux>> F(Nx + 2 * ng, std::vector<Flux>(Ny + 2 * ng));
    std::vector<std::vector<Flux>> G(Nx + 2 * ng, std::vector<Flux>(Ny + 2 * ng));
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            F[i][j] = fluxX(grid.W[i][j], gamma);
            G[i][j] = fluxY(grid.W[i][j], gamma);
        }

    // Predictor (forward differences)
    auto U_pred = U0;
    for (int i = ng; i < Nx + ng; ++i)
        for (int j = ng; j < Ny + ng; ++j) {
            U_pred[i][j].rho = U0[i][j].rho - dt / dx * (F[i + 1][j].rho_f - F[i][j].rho_f) - dt / dy * (G[i][j + 1].rho_f - G[i][j].rho_f);
            U_pred[i][j].rhou = U0[i][j].rhou - dt / dx * (F[i + 1][j].rhou_f - F[i][j].rhou_f) - dt / dy * (G[i][j + 1].rhou_f - G[i][j].rhou_f);
            U_pred[i][j].rhov = U0[i][j].rhov - dt / dx * (F[i + 1][j].rhov_f - F[i][j].rhov_f) - dt / dy * (G[i][j + 1].rhov_f - G[i][j].rhov_f);
            U_pred[i][j].E = U0[i][j].E - dt / dx * (F[i + 1][j].E_f - F[i][j].E_f) - dt / dy * (G[i][j + 1].E_f - G[i][j].E_f);
        }

    // Simple boundary for predictor (copy from U0)
    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            U_pred[i][j] = U0[i][j];
            U_pred[Nx + 2 * ng - 1 - i][j] = U0[Nx + 2 * ng - 1 - i][j];
        }
    for (int j = 0; j < ng; ++j)
        for (int i = 0; i < Nx + 2 * ng; ++i) {
            U_pred[i][j] = U0[i][j];
            U_pred[i][Ny + 2 * ng - 1 - j] = U0[i][Ny + 2 * ng - 1 - j];
        }

    // Convert predictor to primitive
    std::vector<std::vector<State>> W_pred(Nx + 2 * ng, std::vector<State>(Ny + 2 * ng));
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            W_pred[i][j] = consToPhys(U_pred[i][j], gamma);

    // Fluxes for predictor
    std::vector<std::vector<Flux>> Fp(Nx + 2 * ng, std::vector<Flux>(Ny + 2 * ng));
    std::vector<std::vector<Flux>> Gp(Nx + 2 * ng, std::vector<Flux>(Ny + 2 * ng));
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j) {
            Fp[i][j] = fluxX(W_pred[i][j], gamma);
            Gp[i][j] = fluxY(W_pred[i][j], gamma);
        }

    // Corrector (backward differences)
    for (int i = ng; i < Nx + ng; ++i)
        for (int j = ng; j < Ny + ng; ++j) {
            grid.U[i][j].rho = 0.5 * (U0[i][j].rho + U_pred[i][j].rho - dt / dx * (Fp[i][j].rho_f - Fp[i - 1][j].rho_f) - dt / dy * (Gp[i][j].rho_f - Gp[i][j - 1].rho_f));
            grid.U[i][j].rhou = 0.5 * (U0[i][j].rhou + U_pred[i][j].rhou - dt / dx * (Fp[i][j].rhou_f - Fp[i - 1][j].rhou_f) - dt / dy * (Gp[i][j].rhou_f - Gp[i][j - 1].rhou_f));
            grid.U[i][j].rhov = 0.5 * (U0[i][j].rhov + U_pred[i][j].rhov - dt / dx * (Fp[i][j].rhov_f - Fp[i - 1][j].rhov_f) - dt / dy * (Gp[i][j].rhov_f - Gp[i][j - 1].rhov_f));
            grid.U[i][j].E = 0.5 * (U0[i][j].E + U_pred[i][j].E - dt / dx * (Fp[i][j].E_f - Fp[i - 1][j].E_f) - dt / dy * (Gp[i][j].E_f - Gp[i][j - 1].E_f));
        }
}