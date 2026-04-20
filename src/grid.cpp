#include "grid.h"
#include "euler_utils.h"

void initialize_grid(Grid& grid, const Config& cfg) {
    grid.dx = (cfg.grid.x_max - cfg.grid.x_min) / cfg.grid.Nx;
    grid.dy = (cfg.grid.y_max - cfg.grid.y_min) / cfg.grid.Ny;

    // Сначала заполняем координаты центров ячеек (включая фиктивные)
    for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i)
        grid.x_centers[i] = cfg.grid.x_min + (i - grid.num_fict + 0.5) * grid.dx;
    for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j)
        grid.y_centers[j] = cfg.grid.y_min + (j - grid.num_fict + 0.5) * grid.dy;

    // ----- Тест Mader3 (обтекание угла) -----
    if (cfg.phys.test_case == 6) {
        const double x_corner = 2.0;
        const double y_corner = 3.0;
        const double x_det = 0.15;
        const double rho_wall = 1000.0;

        for (int i = 0; i < grid.Nx + 2 * grid.num_fict; ++i) {
            double x = grid.x_centers[i];
            for (int j = 0; j < grid.Ny + 2 * grid.num_fict; ++j) {
                double y = grid.y_centers[j];
                State state;

                // Стенка
                if (x < x_corner && y < y_corner) {
                    state.rho = rho_wall;
                    state.u = 0.0;
                    state.v = 0.0;
                    state.p = 0.0;
                }
                // Детонатор (CJ-продукты)
                else if (x < x_det && y >= y_corner) {
                    state = cfg.phys.left_bottom;   // rho=2.4533, u=0.22, p=0.3562
                }
                // Непрореагировавшее ВВ
                else {
                    state = cfg.phys.right_bottom;   // rho=1.84, p=0
                    state.u = 0.0;
                    state.v = 0.0;
                }

                grid.W[i][j] = state;
                grid.U[i][j] = physToCons(state, cfg.phys.gamma);
            }
        }
        return;
    }

    // ----- Стандартная инициализация с диафрагмой -----
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