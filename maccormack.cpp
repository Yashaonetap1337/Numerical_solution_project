#include "maccormack.h"
#include "euler_utils.h"
#include <vector>

void maccormack_step(Grid& grid, double dt, const Config& cfg) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    std::vector<Conserved> U_n = grid.U;
    std::vector<Flux> F_n(total_cells);
    for (int i = 0; i < total_cells; ++i) F_n[i] = physToFlux(grid.W[i], gamma);

    // ÏÐÅÄÈÊÒÎÐ 
    std::vector<Conserved> U_pred(total_cells);
    for (int i = 0; i < total_cells - 1; ++i) {
        U_pred[i] = U_n[i] - (dt / dx) * (F_n[i + 1] - F_n[i]);
    }
    U_pred[total_cells - 1] = U_pred[total_cells - 2]; // ÃÓ

    // Ïîòîêè ïðåäèêòîðà
    std::vector<Flux> F_pred(total_cells);
    for (int i = 0; i < total_cells; ++i) {
        State W_pred = consToPhys(U_pred[i], gamma);
        F_pred[i] = physToFlux(W_pred, gamma);
    }

    //  ÊÎÐÐÅÊÒÎÐ 
    for (int i = 1; i < total_cells; ++i) {
        Conserved corrector_term = (dt / dx) * (F_pred[i] - F_pred[i - 1]);
        grid.U[i] = 0.5 * (U_n[i] + U_pred[i] - corrector_term);
    }
}