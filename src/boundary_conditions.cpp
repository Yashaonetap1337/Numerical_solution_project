#include "boundary_conditions.h"
#include "euler_utils.h"

void apply_boundary_conditions(Grid& grid, const Config& cfg) {
    const int ng = grid.num_fict;
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;

    // Left boundary
    for (int i = 0; i < ng; ++i) {
        int fict = i;
        int real = ng + (ng - 1 - i);
        switch (cfg.left_boundary) {
        case BoundaryType::WALL:
            for (int j = 0; j < Ny + 2 * ng; ++j) {
                grid.W[fict][j].rho = grid.W[real][j].rho;
                grid.W[fict][j].p = grid.W[real][j].p;
                grid.W[fict][j].u = -grid.W[real][j].u;
                grid.W[fict][j].v = grid.W[real][j].v;
                grid.U[fict][j] = physToCons(grid.W[fict][j], cfg.phys.gamma);
            }
            break;
        case BoundaryType::FREE:
            for (int j = 0; j < Ny + 2 * ng; ++j) {
                grid.W[fict][j] = grid.W[ng][j];
                grid.U[fict][j] = grid.U[ng][j];
            }
            break;
        case BoundaryType::PERIODIC:
            for (int j = 0; j < Ny + 2 * ng; ++j) {
                int src = Nx + ng + i;
                grid.W[fict][j] = grid.W[src][j];
                grid.U[fict][j] = grid.U[src][j];
            }
            break;
        }
    }

    // Right boundary
    for (int i = 0; i < ng; ++i) {
        int fict = Nx + ng + i;
        int real = Nx + ng - 1 - i;
        switch (cfg.right_boundary) {
        case BoundaryType::WALL:
            for (int j = 0; j < Ny + 2 * ng; ++j) {
                grid.W[fict][j].rho = grid.W[real][j].rho;
                grid.W[fict][j].p = grid.W[real][j].p;
                grid.W[fict][j].u = -grid.W[real][j].u;
                grid.W[fict][j].v = grid.W[real][j].v;
                grid.U[fict][j] = physToCons(grid.W[fict][j], cfg.phys.gamma);
            }
            break;
        case BoundaryType::FREE:
            for (int j = 0; j < Ny + 2 * ng; ++j) {
                grid.W[fict][j] = grid.W[Nx + ng - 1][j];
                grid.U[fict][j] = grid.U[Nx + ng - 1][j];
            }
            break;
        case BoundaryType::PERIODIC:
            for (int j = 0; j < Ny + 2 * ng; ++j) {
                int src = i;
                grid.W[fict][j] = grid.W[src][j];
                grid.U[fict][j] = grid.U[src][j];
            }
            break;
        }
    }

    // Bottom boundary
    for (int j = 0; j < ng; ++j) {
        int fict = j;
        int real = ng + (ng - 1 - j);
        switch (cfg.bottom_boundary) {
        case BoundaryType::WALL:
            for (int i = 0; i < Nx + 2 * ng; ++i) {
                grid.W[i][fict].rho = grid.W[i][real].rho;
                grid.W[i][fict].p = grid.W[i][real].p;
                grid.W[i][fict].v = -grid.W[i][real].v;
                grid.W[i][fict].u = grid.W[i][real].u;
                grid.U[i][fict] = physToCons(grid.W[i][fict], cfg.phys.gamma);
            }
            break;
        case BoundaryType::FREE:
            for (int i = 0; i < Nx + 2 * ng; ++i) {
                grid.W[i][fict] = grid.W[i][ng];
                grid.U[i][fict] = grid.U[i][ng];
            }
            break;
        case BoundaryType::PERIODIC:
            for (int i = 0; i < Nx + 2 * ng; ++i) {
                int src = Ny + ng + j;
                grid.W[i][fict] = grid.W[i][src];
                grid.U[i][fict] = grid.U[i][src];
            }
            break;
        }
    }

    // Top boundary
    for (int j = 0; j < ng; ++j) {
        int fict = Ny + ng + j;
        int real = Ny + ng - 1 - j;
        switch (cfg.top_boundary) {
        case BoundaryType::WALL:
            for (int i = 0; i < Nx + 2 * ng; ++i) {
                grid.W[i][fict].rho = grid.W[i][real].rho;
                grid.W[i][fict].p = grid.W[i][real].p;
                grid.W[i][fict].v = -grid.W[i][real].v;
                grid.W[i][fict].u = grid.W[i][real].u;
                grid.U[i][fict] = physToCons(grid.W[i][fict], cfg.phys.gamma);
            }
            break;
        case BoundaryType::FREE:
            for (int i = 0; i < Nx + 2 * ng; ++i) {
                grid.W[i][fict] = grid.W[i][Ny + ng - 1];
                grid.U[i][fict] = grid.U[i][Ny + ng - 1];
            }
            break;
        case BoundaryType::PERIODIC:
            for (int i = 0; i < Nx + 2 * ng; ++i) {
                int src = j;
                grid.W[i][fict] = grid.W[i][src];
                grid.U[i][fict] = grid.U[i][src];
            }
            break;
        }
    }
}