#include "types.h"
#include <stdexcept>

static double minmod(double a, double b) {
    if (a * b <= 0) {
        return 0.0;
    }
    if (a * a < a * b) {
        return a;
    }
    return b;
}

void reconstruction(Grid& grid, const Config& cfg){
    const int total_cells = grid.Nx + 2 * grid.num_fict;
    std::vector<double> delta_rho(total_cells), delta_u(total_cells), delta_p(total_cells), delta_rhou(total_cells), delta_e(total_cells); // временные векторы для хранения наклонов для каждой переменной
    switch (cfg.var_type) {
    case TypesOfVarForReconstruction::CONSERVATIVE:
        for (int i = 1; i < total_cells - 1; ++i) {
            // вычисляем наклоны для каждой переменной
            double drho_L = grid.U[i].rho - grid.U[i - 1].rho;
            double drho_R = grid.U[i + 1].rho - grid.U[i].rho;

            double drhou_L = grid.U[i].rhou - grid.U[i - 1].rhou;
            double drhou_R = grid.U[i + 1].rhou - grid.U[i].rhou;

            double de_L = grid.U[i].E - grid.U[i - 1].E;
            double de_R = grid.U[i + 1].E - grid.U[i].E;

            // применяем ограничитель minmod к каждой паре наклонов
            delta_rho[i] = minmod(drho_L, drho_R);
            delta_rhou[i] = minmod(drhou_L, drhou_R);
            delta_e[i] = minmod(de_L, de_R);
        }

        // цикл по границам для вычисления потоков
        for (int i = 0; i <= grid.Nx; ++i) {
            // глобальные индексы ячеек слева и справа от текущей границы i
            int cell_i = i + grid.num_fict - 1; // это ячейка "i"
            int cell_i_plus_1 = i + grid.num_fict;   // это ячейка "i+1" 

            // собираем состояния W_L и W_R для Задачи Римана на границе i+1/2 

            Conserved U_L_interface;
            U_L_interface.rho = grid.U[cell_i].rho + 0.5 * delta_rho[cell_i];
            U_L_interface.rhou = grid.U[cell_i].rhou + 0.5 * delta_rhou[cell_i];
            U_L_interface.E = grid.U[cell_i].E + 0.5 * delta_e[cell_i];

            Conserved U_R_interface;
            U_R_interface.rho = grid.U[cell_i_plus_1].rho - 0.5 * delta_rho[cell_i_plus_1];
            U_R_interface.rhou = grid.U[cell_i_plus_1].rhou - 0.5 * delta_rhou[cell_i_plus_1];
            U_R_interface.E = grid.U[cell_i_plus_1].E - 0.5 * delta_e[cell_i_plus_1];

        }
    case TypesOfVarForReconstruction::NONCONSERVATIVE:
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

        }
    default:
        throw std::runtime_error("This type of variables is not defined!");
    }
    
}
