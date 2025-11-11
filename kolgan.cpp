#include "kolgan.h"
#include "riemann_solver.h"
#include "euler_utils.h"
#include <vector>
#include <cmath>
#include <algorithm>


// ограничитель наклона minmod
// задача этой функции - выбрать самый подходящий наклон, чтобы предотвратить появление новых пиков (осцилляций) на решении.
static double minmod(double a, double b) {
    if (a * b <= 0) {
        return 0.0; 
    }
    if (a * a < a * b) {
        return a;
    }
    return b;
}


void kolgan_step(Grid& grid, double dt, const Config& cfg, std::vector<Flux>& fluxes) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    // вычисление ограниченных наклонов (ΔQ) для каждой ячейки

    
    std::vector<double> delta_rho(total_cells), delta_u(total_cells), delta_p(total_cells); // временные векторы для хранения наклонов для каждой переменной


    for (int i = 1; i < total_cells - 1; ++i) {
        // вычисляем наклоны для каждой переменной
        double drho_L = grid.W[i].rho - grid.W[i - 1].rho;
        double drho_R = grid.W[i + 1].rho - grid.W[i].rho;

        double du_L = grid.W[i].u - grid.W[i - 1].u;
        double du_R = grid.W[i + 1].u - grid.W[i].u;

        double dp_L = grid.W[i].p - grid.W[i - 1].p;
        double dp_R = grid.W[i + 1].p - grid.W[i].p;

        // применяем ограничитель minmod к каждой паре наклонов
        delta_rho[i] = minmod(drho_L, drho_R);
        delta_u[i] = minmod(du_L, du_R);
        delta_p[i] = minmod(dp_L, dp_R);
    }

    // цикл по границам для вычисления потоков
    for (int i = 0; i <= grid.Nx; ++i) {
        // глобальные индексы ячеек слева и справа от текущей границы i
        int cell_i = i + grid.num_fict - 1; // это ячейка "i"
        int cell_i_plus_1 = i + grid.num_fict;   // это ячейка "i+1" 

        // собираем состояния W_L и W_R для Задачи Римана на границе i+1/2 

        State W_L_interface;
        W_L_interface.rho = grid.W[cell_i].rho + 0.5 * delta_rho[cell_i];
        W_L_interface.u = grid.W[cell_i].u + 0.5 * delta_u[cell_i];
        W_L_interface.p = grid.W[cell_i].p + 0.5 * delta_p[cell_i];

        State W_R_interface;
        W_R_interface.rho = grid.W[cell_i_plus_1].rho - 0.5 * delta_rho[cell_i_plus_1];
        W_R_interface.u = grid.W[cell_i_plus_1].u - 0.5 * delta_u[cell_i_plus_1];
        W_R_interface.p = grid.W[cell_i_plus_1].p - 0.5 * delta_p[cell_i_plus_1];

        // решаем задачу Римана с этими реконструированными значениями
        const State state_at_interface = solve_general_riemann_problem(W_L_interface, W_R_interface, 0.0, cfg, cfg.approx_type);

        // вычисляем поток
        fluxes[i] = physToFlux(state_at_interface, gamma);
    }

    // обновление консервативных переменных в ячейках
    for (int i = 0; i < grid.Nx; ++i) {
        const int cell_idx = i + grid.num_fict;
        const Flux& F_left = fluxes[i];
        const Flux& F_right = fluxes[i + 1];

        grid.U[cell_idx].rho -= (dt / dx) * (F_right.rho_f - F_left.rho_f);
        grid.U[cell_idx].rhou -= (dt / dx) * (F_right.rhou_f - F_left.rhou_f);
        grid.U[cell_idx].E -= (dt / dx) * (F_right.E_f - F_left.E_f);
    }
}
