#include "eno.h"
#include "choice_of_riemann_solvers.h" 
#include "euler_utils.h"
#include <vector>
#include <cmath>
#include <algorithm>


static void eno3_reconstruction(const std::vector<double>& v, int i,
    double& val_L, double& val_R)
{

    // Шаблон S2 
    double smoothness_S2 = std::abs((v[i] - v[i - 1]) - (v[i - 1] - v[i - 2]));
    // Шаблон S1
    double smoothness_S1 = std::abs((v[i + 1] - v[i]) - (v[i] - v[i - 1]));
    // Шаблон S0 
    double smoothness_S0 = std::abs((v[i + 2] - v[i + 1]) - (v[i + 1] - v[i]));

    //const double v_m2 = v[i - 2], v_m1 = v[i - 1], v_0 = v[i], v_p1 = v[i + 1], v_p2 = v[i + 2];

    //// Полиномы для 3-х шаблонов
    //const double p0 = (2.0 / 6.0) * v_m2 - (7.0 / 6.0) * v_m1 + (11.0 / 6.0) * v_0;
    //const double p1 = (-1.0 / 6.0) * v_m1 + (5.0 / 6.0) * v_0 + (2.0 / 6.0) * v_p1;
    //const double p2 = (2.0 / 6.0) * v_0 + (5.0 / 6.0) * v_p1 - (1.0 / 6.0) * v_p2;

    if (smoothness_S2 < smoothness_S1 && smoothness_S2 < smoothness_S0) {
        // шаблон S2 
        val_L = (1.0 / 3.0) * v[i - 2] - (7.0 / 6.0) * v[i - 1] + (11.0 / 6.0) * v[i];
    }
    else if (smoothness_S1 < smoothness_S2 && smoothness_S1 < smoothness_S0) {
        // шаблон S1 
        val_L = (-1.0 / 6.0) * v[i - 1] + (5.0 / 6.0) * v[i] + (2.0 / 6.0) * v[i + 1];
    }
    else {
        //  шаблон S0 
        val_L = (1.0 / 3.0) * v[i] + (5.0 / 6.0) * v[i + 1] - (1.0 / 6.0) * v[i + 2];
    }



    // Шаблон S2 
    double smoothness_S2_R = std::abs((v[i+1] - v[i ]) - (v[i ] - v[i - 1]));
    // Шаблон S1
    double smoothness_S1_R = std::abs((v[i + 2] - v[i+1]) - (v[i+1] - v[i]));
    // Шаблон S0 
    double smoothness_S0_R = std::abs((v[i + 3] - v[i + 2]) - (v[i + 2] - v[i+1]));

    //const double v_m2 = v[i - 2], v_m1 = v[i - 1], v_0 = v[i], v_p1 = v[i + 1], v_p2 = v[i + 2];

    //// Полиномы для 3-х шаблонов 
    //const double p0 = (-1.0 / 6.0) * v_m2 + (5.0 / 6.0) * v_m1 + (2.0 / 6.0) * v_0;
    //const double p1 = (2.0 / 6.0) * v_m1 + (5.0 / 6.0) * v_0 - (1.0 / 6.0) * v_p1;
    //const double p2 = (11.0 / 6.0) * v_0 - (7.0 / 6.0) * v_p1 + (2.0 / 6.0) * v_p2;



    if (smoothness_S2_R < smoothness_S1_R && smoothness_S2_R < smoothness_S0_R) {
        // шаблон S2 
        val_R = (1.0 / 3.0) * v[i - 1] - (7.0 / 6.0) * v[i] + (11.0 / 6.0) * v[i+1];
    }
    else if (smoothness_S1_R < smoothness_S2_R && smoothness_S1_R < smoothness_S0_R) {
        // шаблон S1
        val_R = (-1.0 / 6.0) * v[i] + (5.0 / 6.0) * v[i+1] + (2.0 / 6.0) * v[i + 2];
    }
    else {
        // шаблон S0 
        val_R = (1.0 / 3.0) * v[i+1] + (5.0 / 6.0) * v[i + 2] - (1.0 / 6.0) * v[i + 3];
    }


}


void eno_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    std::vector<double> rho(total_cells), u(total_cells), p(total_cells);
    for (int i = 0; i < total_cells; ++i) {
        rho[i] = grid.W[i].rho; u[i] = grid.W[i].u; p[i] = grid.W[i].p;
    }

    for (int i = 0; i <= grid.Nx; ++i) {
        int cell_i = i + grid.num_fict - 1;


        if (cell_i < 2 || cell_i > total_cells - 4) {
            const State W_L = grid.W[cell_i];
            const State W_R = grid.W[cell_i + 1];
            const State state_at_interface = solve_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
            fluxes[i] = physToFlux(state_at_interface, gamma);
            continue;
        }

        State W_L_interface, W_R_interface;

        eno3_reconstruction(rho, cell_i, W_L_interface.rho, W_R_interface.rho);
        eno3_reconstruction(u, cell_i, W_L_interface.u, W_R_interface.u);
        eno3_reconstruction(p, cell_i, W_L_interface.p, W_R_interface.p);


        const double EPS = 1e-9;
        if (W_L_interface.rho < EPS) W_L_interface.rho = grid.W[cell_i].rho;
        if (W_L_interface.p < EPS)   W_L_interface.p = grid.W[cell_i].p;
        if (W_R_interface.rho < EPS) W_R_interface.rho = grid.W[cell_i + 1].rho;
        if (W_R_interface.p < EPS)   W_R_interface.p = grid.W[cell_i + 1].p;


        const State state_at_interface = solve_riemann_problem(W_L_interface, W_R_interface, 0.0, cfg, cfg.approx_type);
        fluxes[i] = physToFlux(state_at_interface, gamma);
    }
}