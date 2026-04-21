#include "flic.h"
#include "euler_utils.h"
#include "boundary_conditions.h"
#include <vector>
#include <algorithm>
#include <cmath>

// Псевдоним для поля примитивных переменных (как W_tilde в эталоне)
using StateField = std::vector<std::vector<State>>;

static StateField make_state_field(int nx, int ny) {
    return StateField(nx, std::vector<State>(ny));
}

static void apply_bc_to_tilde(StateField& W_tilde, Grid& grid, const Config& cfg) {
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
            grid.W[i][j] = W_tilde[i][j];

    apply_boundary_conditions(grid, cfg);

    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
            W_tilde[i][j] = grid.W[i][j];
}

// CHECK: FLIC_LAGRANGE
static void flic_lagrange(const Grid& grid, const Config& cfg,
    StateField& W_tilde, double dt, int dir)
{
    const int ng  = grid.num_fict;
    const int Nx  = grid.Nx;
    const int Ny  = grid.Ny;
    const double gamma = cfg.phys.gamma;

    for (int i = 0; i < Nx + 2 * ng; ++i)
        for (int j = 0; j < Ny + 2 * ng; ++j)
            W_tilde[i][j] = grid.W[i][j];

    if (dir == 0) {
        for (int i = ng; i < Nx + ng - 1; ++i) {
            for (int j = ng; j < Ny + ng; ++j) {
                const double dx = grid.dx;

                const double rho = grid.W[i][j].rho;
                const double u   = grid.W[i][j].u;
                const double v   = grid.W[i][j].v;
                const double p   = grid.W[i][j].p;

                const double p_ip = 0.5 * (p + grid.W[i + 1][j].p);
                const double p_im = 0.5 * (grid.W[i - 1][j].p + p);

                const double u_ip = 0.5 * (u + grid.W[i + 1][j].u);
                const double u_im = 0.5 * (grid.W[i - 1][j].u + u);

                W_tilde[i][j].rho = rho;

                const double u_tilde = u - dt / rho * (p_ip - p_im) / dx;
                W_tilde[i][j].u = u_tilde;

                W_tilde[i][j].v = v;

                const double E = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
                const double E_tilde = E - dt * (p_ip * u_ip - p_im * u_im) / dx;

                const double kinetic_tilde = 0.5 * rho * (u_tilde * u_tilde + v * v);
                const double p_tilde = (gamma - 1.0) * (E_tilde - kinetic_tilde);
                W_tilde[i][j].p = std::max(1e-8, p_tilde);
            }
        }
    }

    if (dir == 1) {
        for (int i = ng; i < Nx + ng; ++i) {
            for (int j = ng; j < Ny + ng - 1; ++j) {
                const double dy = grid.dy;

                const double rho = grid.W[i][j].rho;
                const double u   = grid.W[i][j].u;
                const double v   = grid.W[i][j].v;
                const double p   = grid.W[i][j].p;

                const double p_jp = 0.5 * (p + grid.W[i][j + 1].p);
                const double p_jm = 0.5 * (grid.W[i][j - 1].p + p);

                const double v_jp = 0.5 * (v + grid.W[i][j + 1].v);
                const double v_jm = 0.5 * (grid.W[i][j - 1].v + v);

                W_tilde[i][j].rho = rho;

                W_tilde[i][j].u = u;

                const double dPdy = (p_jp - p_jm) / dy;
                const double v_tilde = v - dt / rho * dPdy;
                W_tilde[i][j].v = v_tilde;

                const double E = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
                const double dPVdy = (p_jp * v_jp - p_jm * v_jm) / dy;
                const double E_tilde = E - dt * dPVdy;

                const double kinetic_tilde = 0.5 * rho * (u * u + v_tilde * v_tilde);
                const double p_tilde = (gamma - 1.0) * (E_tilde - kinetic_tilde);
                W_tilde[i][j].p = std::max(1e-8, p_tilde);
            }
        }
    }
}

// CHECK: FLIC_EULER
static void flic_euler(const StateField& W_tilde, Grid& grid,
    const Config& cfg, double dt, int dir)
{
    const int ng  = grid.num_fict;
    const int Nx  = grid.Nx;
    const int Ny  = grid.Ny;
    const double gamma    = cfg.phys.gamma;
    const double vel_eps  = 1e-12;

    if (dir == 0) {
        for (int i = ng; i < Nx + ng - 1; ++i) {
            for (int j = ng; j < Ny + ng; ++j) {
                const double dx = grid.dx;

                double rho_flux_R = 0.0, rhou_flux_R = 0.0;
                double rhov_flux_R = 0.0, E_flux_R = 0.0;

                const double u_face_R = 0.5 * (W_tilde[i][j].u + W_tilde[i + 1][j].u);
                if (std::abs(u_face_R) > vel_eps) {
                    const int up = (u_face_R > 0.0) ? i : i + 1;
                    const double rho_d = W_tilde[up][j].rho;
                    const double u_d   = W_tilde[up][j].u;
                    const double v_d   = W_tilde[up][j].v;
                    const double p_d   = W_tilde[up][j].p;
                    const double E_d   = p_d / (gamma - 1.0)
                                       + 0.5 * rho_d * (u_d * u_d + v_d * v_d);
                    rho_flux_R  = rho_d * u_d;
                    rhou_flux_R = rho_flux_R * u_d;
                    rhov_flux_R = rho_flux_R * v_d;
                    E_flux_R    = u_d * E_d;
                }

                double rho_flux_L = 0.0, rhou_flux_L = 0.0;
                double rhov_flux_L = 0.0, E_flux_L = 0.0;

                const double u_face_L = 0.5 * (W_tilde[i - 1][j].u + W_tilde[i][j].u);
                if (std::abs(u_face_L) > vel_eps) {
                    const int up = (u_face_L > 0.0) ? i - 1 : i;
                    const double rho_d = W_tilde[up][j].rho;
                    const double u_d   = W_tilde[up][j].u;
                    const double v_d   = W_tilde[up][j].v;
                    const double p_d   = W_tilde[up][j].p;
                    const double E_d   = p_d / (gamma - 1.0)
                                       + 0.5 * rho_d * (u_d * u_d + v_d * v_d);
                    rho_flux_L  = rho_d * u_d;
                    rhou_flux_L = rho_flux_L * u_d;
                    rhov_flux_L = rho_flux_L * v_d;
                    E_flux_L    = u_d * E_d;
                }

                const double rho_t = W_tilde[i][j].rho;
                const double u_t   = W_tilde[i][j].u;
                const double v_t   = W_tilde[i][j].v;
                const double p_t   = W_tilde[i][j].p;
                const double E_t   = p_t / (gamma - 1.0)
                                   + 0.5 * rho_t * (u_t * u_t + v_t * v_t);

                // CHECK: FLIC_CONSERV
                double rho_new  = rho_t             - dt / dx * (rho_flux_R  - rho_flux_L);
                double rhou_new = rho_t * u_t       - dt / dx * (rhou_flux_R - rhou_flux_L);
                double rhov_new = rho_t * v_t       - dt / dx * (rhov_flux_R - rhov_flux_L);
                double E_new    = E_t               - dt / dx * (E_flux_R    - E_flux_L);

                rho_new = std::max(1e-8, rho_new);

                double u_new = rhou_new / rho_new;
                double v_new = rhov_new / rho_new;
                if (std::abs(u_new) < 1e-9) u_new = 0.0;
                if (std::abs(v_new) < 1e-9) v_new = 0.0;

                const double kinetic = 0.5 * rho_new * (u_new * u_new + v_new * v_new);
                const double p_new   = (gamma - 1.0) * (E_new - kinetic);

                grid.W[i][j].rho = rho_new;
                grid.W[i][j].u   = u_new;
                grid.W[i][j].v   = v_new;
                grid.W[i][j].p   = std::max(1e-8, p_new);
                grid.U[i][j]     = physToCons(grid.W[i][j], gamma);
            }
        }
    }

    if (dir == 1) {
        for (int i = ng; i < Nx + ng; ++i) {
            for (int j = ng; j < Ny + ng - 1; ++j) {
                const double dy = grid.dy;

                double rho_flux_T = 0.0, rhou_flux_T = 0.0;
                double rhov_flux_T = 0.0, E_flux_T = 0.0;

                const double v_face_T = 0.5 * (W_tilde[i][j].v + W_tilde[i][j + 1].v);
                if (std::abs(v_face_T) > vel_eps) {
                    const int up = (v_face_T > 0.0) ? j : j + 1;
                    const double rho_d = W_tilde[i][up].rho;
                    const double u_d   = W_tilde[i][up].u;
                    const double v_d   = W_tilde[i][up].v;
                    const double p_d   = W_tilde[i][up].p;
                    const double E_d   = p_d / (gamma - 1.0)
                                       + 0.5 * rho_d * (u_d * u_d + v_d * v_d);
                    rho_flux_T  = rho_d * v_d;
                    rhou_flux_T = rho_flux_T * u_d;
                    rhov_flux_T = rho_flux_T * v_d;
                    E_flux_T    = v_d * E_d;
                }

                double rho_flux_B = 0.0, rhou_flux_B = 0.0;
                double rhov_flux_B = 0.0, E_flux_B = 0.0;

                const double v_face_B = 0.5 * (W_tilde[i][j - 1].v + W_tilde[i][j].v);
                if (std::abs(v_face_B) > vel_eps) {
                    const int up = (v_face_B > 0.0) ? j - 1 : j;
                    const double rho_d = W_tilde[i][up].rho;
                    const double u_d   = W_tilde[i][up].u;
                    const double v_d   = W_tilde[i][up].v;
                    const double p_d   = W_tilde[i][up].p;
                    const double E_d   = p_d / (gamma - 1.0)
                                       + 0.5 * rho_d * (u_d * u_d + v_d * v_d);
                    rho_flux_B  = rho_d * v_d;
                    rhou_flux_B = rho_flux_B * u_d;
                    rhov_flux_B = rho_flux_B * v_d;
                    E_flux_B    = v_d * E_d;
                }

                const double rho_t = W_tilde[i][j].rho;
                const double u_t   = W_tilde[i][j].u;
                const double v_t   = W_tilde[i][j].v;
                const double p_t   = W_tilde[i][j].p;
                const double E_t   = p_t / (gamma - 1.0)
                                   + 0.5 * rho_t * (u_t * u_t + v_t * v_t);

                double rho_new  = rho_t       - dt / dy * (rho_flux_T  - rho_flux_B);
                double rhou_new = rho_t * u_t - dt / dy * (rhou_flux_T - rhou_flux_B);
                double rhov_new = rho_t * v_t - dt / dy * (rhov_flux_T - rhov_flux_B);
                double E_new    = E_t         - dt / dy * (E_flux_T    - E_flux_B);

                rho_new = std::max(1e-8, rho_new);

                double u_new = rhou_new / rho_new;
                double v_new = rhov_new / rho_new;
                if (std::abs(u_new) < 1e-9) u_new = 0.0;
                if (std::abs(v_new) < 1e-9) v_new = 0.0;

                const double kinetic = 0.5 * rho_new * (u_new * u_new + v_new * v_new);
                const double p_new   = (gamma - 1.0) * (E_new - kinetic);

                grid.W[i][j].rho = rho_new;
                grid.W[i][j].u   = u_new;
                grid.W[i][j].v   = v_new;
                grid.W[i][j].p   = std::max(1e-8, p_new);
                grid.U[i][j]     = physToCons(grid.W[i][j], gamma);
            }
        }
    }
}

static void flic_substep(Grid& grid, const Config& cfg,
    StateField& W_tilde, double dt, int dir)
{
    flic_lagrange(grid, cfg, W_tilde, dt, dir);
    apply_bc_to_tilde(W_tilde, grid, cfg);
    flic_euler(W_tilde, grid, cfg, dt, dir);
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
            grid.U[i][j] = physToCons(grid.W[i][j], cfg.phys.gamma);
    apply_boundary_conditions(grid, cfg);
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
            grid.W[i][j] = consToPhys(grid.U[i][j], cfg.phys.gamma);
}

void flic_step(Grid& grid, double dt, const Config& cfg) {
    static bool swap = false;

    StateField W_tilde = make_state_field(
        grid.Nx + 2 * grid.num_fict,
        grid.Ny + 2 * grid.num_fict);

    const double dt_half = 0.5 * dt;

    if (!swap) {
        flic_substep(grid, cfg, W_tilde, dt_half, 0);
        flic_substep(grid, cfg, W_tilde, dt,      1);
        flic_substep(grid, cfg, W_tilde, dt_half, 0);
    } else {
        flic_substep(grid, cfg, W_tilde, dt_half, 1);
        flic_substep(grid, cfg, W_tilde, dt,      0);
        flic_substep(grid, cfg, W_tilde, dt_half, 1);
    }

    const int ng = grid.num_fict;
    for (int i = ng; i < grid.Nx + ng; ++i) {
        for (int j = ng; j < grid.Ny + ng; ++j) {
            if (std::abs(grid.W[i][j].u) < 1e-9) grid.W[i][j].u = 0.0;
            if (std::abs(grid.W[i][j].v) < 1e-9) grid.W[i][j].v = 0.0;
            grid.U[i][j] = physToCons(grid.W[i][j], cfg.phys.gamma);
        }
    }

    swap = !swap;
}
