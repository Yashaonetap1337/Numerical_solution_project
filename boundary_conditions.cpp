#include "boundary_conditions.h"
#include "euler_utils.h"
#include <iostream>

void apply_boundary_conditions(Grid& grid, const Config& cfg) {
    const double gamma = cfg.phys.gamma;
    const int Nf = grid.num_fict;
    const int Nx = grid.Nx;

    switch (cfg.left_boundary) {
    case BoundaryType::WALL: {
        for (int i = 0; i < Nf; ++i) {
            int fict_idx = Nf - 1 - i; // Индекс фиктивной ячейки (снаружи)
            int real_idx = Nf + i;     // Индекс "зеркальной" реальной ячейки (внутри)


            grid.W[fict_idx].rho = grid.W[real_idx].rho;
            grid.W[fict_idx].p = grid.W[real_idx].p;
            grid.W[fict_idx].u = -grid.W[real_idx].u;

            // Обновляем консервативные переменные
            grid.U[fict_idx].rho = grid.U[real_idx].rho;
            grid.U[fict_idx].rhou = -grid.U[real_idx].rhou; 
            grid.U[fict_idx].E = grid.U[real_idx].E;    
        }
        break;
    }
    case BoundaryType::FREE: {
        for (int i = 0; i < Nf; ++i) {
            int fict_idx = i;
            int real_idx = Nf; 

            grid.W[fict_idx] = grid.W[real_idx];
            grid.U[fict_idx] = grid.U[real_idx];
        }
        break;
    }
    case BoundaryType::PERIODIC: {
        for (int i = 0; i < Nf; ++i) {
            int fict_idx = i;
            int real_idx = Nx + i; 

            grid.W[fict_idx] = grid.W[real_idx];
            grid.U[fict_idx] = grid.U[real_idx];
        }
        break;
    }
    }


    switch (cfg.right_boundary) {
    case BoundaryType::WALL: {

        for (int i = 0; i < Nf; ++i) {
            int fict_idx = Nf + Nx + i;     
            int real_idx = Nf + Nx - 1 - i; 


            grid.W[fict_idx].rho = grid.W[real_idx].rho;
            grid.W[fict_idx].p = grid.W[real_idx].p;
            grid.W[fict_idx].u = -grid.W[real_idx].u;


            grid.U[fict_idx].rho = grid.U[real_idx].rho;
            grid.U[fict_idx].rhou = -grid.U[real_idx].rhou;
            grid.U[fict_idx].E = grid.U[real_idx].E;
        }
        break;
    }
    case BoundaryType::FREE: {

        for (int i = 0; i < Nf; ++i) {
            int fict_idx = Nf + Nx + i;
            int real_idx = Nf + Nx - 1; 

            grid.W[fict_idx] = grid.W[real_idx];
            grid.U[fict_idx] = grid.U[real_idx];
        }
        break;
    }
    case BoundaryType::PERIODIC: {

        for (int i = 0; i < Nf; ++i) {

            int fict_idx = Nf + Nx + i;
            int real_idx = Nf + i;

            grid.W[fict_idx] = grid.W[real_idx];
            grid.U[fict_idx] = grid.U[real_idx];
        }
        break;
    }
    }
}