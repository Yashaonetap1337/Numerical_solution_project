#include "grid.h"
#include "euler_utils.h"

#include <iostream>


void initialize_grid(Grid& grid, const Config& cfg) {
    grid.dx = (cfg.grid.x_max - cfg.grid.x_min) / cfg.grid.Nx;

    // цикл по всем ячейкам, включая фиктивные
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i) {
        // координата центра ячейки i (с учетом смещения из-за фиктивных ячеек)
        grid.x_centers[i] = cfg.grid.x_min + (i - grid.num_fict + 0.5) * grid.dx;

        // задаем начальные условия внутри основной сетки (фиктивные пока не трогаем)
        if (grid.x_centers[i] < cfg.grid.x_diaphragm) {
            grid.W[i] = cfg.phys.left;
        }
        else {
            grid.W[i] = cfg.phys.right;
        }
        // переводим начальные условия в консервативные переменные
        grid.U[i] = physToCons(grid.W[i], cfg.phys.gamma);
    }
    std::cout << "Grid initialized. Nx = " << grid.Nx
        << ", fict = " << grid.num_fict
        << ", dx = " << grid.dx << std::endl;
    std::cout << std::endl;
}