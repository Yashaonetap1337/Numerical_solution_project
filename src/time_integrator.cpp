#include "time_integrator.h"
#include "euler_utils.h"
#include "boundary_conditions.h"
#include "godunov.h"
#include "acoustic.h"
#include "kolgan.h"
#include <vector>
#include <stdexcept>

static void compute_rhs(Grid& grid, const Config& cfg,
    std::vector<std::vector<Conserved>>& rhs) {
    const double dx = grid.dx, dy = grid.dy;
    const int ng = grid.num_fict;
    const int Nx = grid.Nx, Ny = grid.Ny;

    std::vector<std::vector<Flux>> Fx, Fy;
    switch (cfg.method) {
    case NumericalMethod::GODUNOV:
        godunov_flux_computation(grid, cfg, Fx, Fy);
        break;
    case NumericalMethod::ACOUSTIC:
        acoustic_flux_computation(grid, cfg, Fx, Fy);
        break;
    case NumericalMethod::KOLGAN:
        kolgan_flux_computation(grid, cfg, Fx, Fy);
        break;
    default:
        throw std::runtime_error("Unsupported method in compute_rhs");
    }

    rhs.assign(Nx + 2 * ng, std::vector<Conserved>(Ny + 2 * ng, { 0,0,0,0 }));
    for (int i = ng; i < Nx + ng; ++i) {
        for (int j = ng; j < Ny + ng; ++j) {
            const Flux& FL = Fx[i - ng][j];
            const Flux& FR = Fx[i - ng + 1][j];
            const Flux& GB = Fy[i][j - ng];
            const Flux& GT = Fy[i][j - ng + 1];
            rhs[i][j].rho = -((FR.rho_f - FL.rho_f) / dx + (GT.rho_f - GB.rho_f) / dy);
            rhs[i][j].rhou = -((FR.rhou_f - FL.rhou_f) / dx + (GT.rhou_f - GB.rhou_f) / dy);
            rhs[i][j].rhov = -((FR.rhov_f - FL.rhov_f) / dx + (GT.rhov_f - GB.rhov_f) / dy);
            rhs[i][j].E = -((FR.E_f - FL.E_f) / dx + (GT.E_f - GB.E_f) / dy);
        }
    }
}

static void euler_step(Grid& grid, double dt, const Config& cfg) {
    std::vector<std::vector<Conserved>> rhs;
    compute_rhs(grid, cfg, rhs);
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
            grid.U[i][j] = grid.U[i][j] + dt * rhs[i][j];
}

static void tvd_rk3_step(Grid& grid, double dt, const Config& cfg) {
    const int ng = grid.num_fict;
    const int Nx = grid.Nx, Ny = grid.Ny;
    auto U0 = grid.U;

    std::vector<std::vector<Conserved>> rhs1;
    compute_rhs(grid, cfg, rhs1);
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            grid.U[i][j] = U0[i][j] + dt * rhs1[i][j];
    apply_boundary_conditions(grid, cfg);
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            grid.W[i][j] = consToPhys(grid.U[i][j], cfg.phys.gamma);

    std::vector<std::vector<Conserved>> rhs2;
    compute_rhs(grid, cfg, rhs2);
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            grid.U[i][j] = 0.75 * U0[i][j] + 0.25 * grid.U[i][j] + 0.25 * dt * rhs2[i][j];
    apply_boundary_conditions(grid, cfg);
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            grid.W[i][j] = consToPhys(grid.U[i][j], cfg.phys.gamma);

    std::vector<std::vector<Conserved>> rhs3;
    compute_rhs(grid, cfg, rhs3);
    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            grid.U[i][j] = (1.0 / 3.0) * U0[i][j] + (2.0 / 3.0) * grid.U[i][j] + (2.0 / 3.0) * dt * rhs3[i][j];
}

void time_step(Grid& grid, double dt, const Config& cfg) {
    switch (cfg.time_integrator) {
    case TimeIntegrator::EULER:   euler_step(grid, dt, cfg); break;
    case TimeIntegrator::TVD_RK3: tvd_rk3_step(grid, dt, cfg); break;
    default: throw std::runtime_error("Unknown time integrator");
    }
    apply_boundary_conditions(grid, cfg);
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
            grid.W[i][j] = consToPhys(grid.U[i][j], cfg.phys.gamma);
}