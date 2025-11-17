#include "boundary_conditions.h"
#include "euler_utils.h"
#include <iostream>

void apply_boundary_conditions(Grid& grid, const Config& cfg) {
    const double gamma = cfg.phys.gamma;

    // Левая граница
    switch (cfg.left_boundary) {
    case BoundaryType::WALL: {
        // Отражающая стенка
        for (int i = 0; i < grid.num_fict; ++i) {
            int fict_idx = grid.num_fict - 1 - i;
            int real_idx = grid.num_fict + i;

            grid.W[fict_idx].rho = grid.W[real_idx].rho;
            grid.W[fict_idx].p = grid.W[real_idx].p;
            grid.W[fict_idx].u = -grid.W[real_idx].u;

            grid.U[fict_idx] = physToCons(grid.W[fict_idx], gamma);
        }
        break;
    }
    case BoundaryType::FREE: {
        // Свободная граница (транслирующая)
        for (int i = 0; i < grid.num_fict; ++i) {
            int fict_idx = i;
            int real_idx = grid.num_fict;

            grid.W[fict_idx] = grid.W[real_idx];
            grid.U[fict_idx] = grid.U[real_idx];
        }
        break;
    }
    case BoundaryType::PERIODIC: {
        // Периодическая граница
        for (int i = 0; i < grid.num_fict; ++i) {
            grid.W[i] = grid.W[grid.Nx + i];
            grid.U[i] = grid.U[grid.Nx + i];
        }
        break;
    }
    }

    // Правая граница
    switch (cfg.right_boundary) {
    case BoundaryType::WALL: {
        // Отражающая стенка
        for (int i = 0; i < grid.num_fict; ++i) {
            int fict_idx = grid.num_fict + grid.Nx + i;
            int real_idx = grid.num_fict + grid.Nx - 1 - i;

            grid.W[fict_idx].rho = grid.W[real_idx].rho;
            grid.W[fict_idx].p = grid.W[real_idx].p;
            grid.W[fict_idx].u = -grid.W[real_idx].u;

            grid.U[fict_idx] = physToCons(grid.W[fict_idx], gamma);
        }
        break;
    }
    case BoundaryType::FREE: {
        // Свободная граница (транслирующая)
        for (int i = 0; i < grid.num_fict; ++i) {
            int fict_idx = grid.num_fict + grid.Nx + i;
            int real_idx = grid.num_fict + grid.Nx - 1;

            grid.W[fict_idx] = grid.W[real_idx];
            grid.U[fict_idx] = grid.U[real_idx];
        }
        break;
    }
    case BoundaryType::PERIODIC: {
        // Периодическая граница
        for (int i = 0; i < grid.num_fict; ++i) {
            grid.W[grid.num_fict + grid.Nx + i] = grid.W[grid.num_fict + i];
            grid.U[grid.num_fict + grid.Nx + i] = grid.U[grid.num_fict + i];
        }
        break;
    }
    }
}