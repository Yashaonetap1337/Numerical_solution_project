#include "boundary_conditions.h"

void apply_boundary_conditions(Grid& grid) {
    // транслирующие √”

    // копируем состо€ние из первой реальной €чейки во все левые фиктивные.
    // (это создает нулевой градиент на границе, позвол€€ потоку вытекать).
    for (int i = 0; i < grid.num_fict; ++i) {
        int fict_idx = i;
        int real_idx = grid.num_fict;

        grid.W[fict_idx] = grid.W[real_idx];
        grid.U[fict_idx] = grid.U[real_idx];
    }

    
    // копируем состо€ние из последней реальной €чейки во все правые фиктивные.
    for (int i = 0; i < grid.num_fict; ++i) {
        int fict_idx = grid.Nx + grid.num_fict + i;
        int real_idx = grid.Nx + grid.num_fict - 1;

        grid.W[fict_idx] = grid.W[real_idx];
        grid.U[fict_idx] = grid.U[real_idx];
    }




    //// отражающие стенки
    //for (int i = 0; i < grid.num_fict; ++i) {
    //    
    //    int fict_idx = grid.num_fict - 1 - i;
    //    int real_idx = grid.num_fict + i;

    //    
    //    grid.W[fict_idx].rho = grid.W[real_idx].rho;
    //    grid.W[fict_idx].p = grid.W[real_idx].p;
    //    grid.W[fict_idx].u = -grid.W[real_idx].u;

    //    
    //    grid.U[fict_idx].rho = grid.U[real_idx].rho;
    //    grid.U[fict_idx].E = grid.U[real_idx].E;
    //    grid.U[fict_idx].rhou = -grid.U[real_idx].rhou;
    //}

    //for (int i = 0; i < grid.num_fict; ++i) {
    //    
    //    int fict_idx = grid.Nx + grid.num_fict + i;
    //    int real_idx = grid.Nx + grid.num_fict - 1 - i;
    //    
    // 
    //    grid.W[fict_idx].rho = grid.W[real_idx].rho;
    //    grid.W[fict_idx].p = grid.W[real_idx].p;
    //    grid.W[fict_idx].u = -grid.W[real_idx].u;

    //    grid.U[fict_idx].rho = grid.U[real_idx].rho;
    //    grid.U[fict_idx].E = grid.U[real_idx].E;
    //    grid.U[fict_idx].rhou = -grid.U[real_idx].rhou;
    //}
}