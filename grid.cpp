#include "grid.h"
#include "euler_utils.h"

#include <iostream>


void initialize_grid(Grid& grid, const Config& cfg) {
    grid.dx = (cfg.grid.x_max - cfg.grid.x_min) / cfg.grid.Nx;
;

    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i) {
        grid.x_centers[i] = cfg.grid.x_min + (i - grid.num_fict + 0.5) * grid.dx;

        if (grid.x_centers[i] < cfg.grid.x_diaphragm) {
            grid.W[i] = cfg.phys.left;
        }
        else {
            grid.W[i] = cfg.phys.right;
        }

        grid.U[i] = physToCons(grid.W[i], cfg.phys.gamma);
    }
}