#include "grid.h"
#include "euler_utils.h"

void initialize_grid(Grid& grid, const Config& cfg) {
    grid.dx = (cfg.grid.x_max - cfg.grid.x_min) / cfg.grid.Nx;
    grid.dy = (cfg.grid.y_max - cfg.grid.y_min) / cfg.grid.Ny;

    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        grid.x_centers[i] = cfg.grid.x_min + (i - grid.num_fict + 0.5) * grid.dx;
    for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
        grid.y_centers[j] = cfg.grid.y_min + (j - grid.num_fict + 0.5) * grid.dy;

    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i) {
        double x = grid.x_centers[i];
        for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j) {
            double y = grid.y_centers[j];
            State state;
            if (x < cfg.grid.x_diaphragm && y < cfg.grid.y_diaphragm)
                state = cfg.phys.left_bottom;
            else if (x < cfg.grid.x_diaphragm && y >= cfg.grid.y_diaphragm)
                state = cfg.phys.left_top;
            else if (x >= cfg.grid.x_diaphragm && y < cfg.grid.y_diaphragm)
                state = cfg.phys.right_bottom;
            else
                state = cfg.phys.right_top;
            grid.W[i][j] = state;
            grid.U[i][j] = physToCons(state, cfg.phys.gamma);
        }
    }
}