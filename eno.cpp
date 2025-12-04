#include "eno.h"
#include "euler_utils.h"
#include "choice_of_riemann_solvers.h" 
#include <cmath>
#include <algorithm>
#include <iostream>


static int determine_stencil_shift(const std::vector<double>& v, int i) {

    double d1_left = std::abs(v[i] - v[i - 1]);
    double d1_right = std::abs(v[i + 1] - v[i]);

    int k_min = i;
    if (d1_left < d1_right) {
        k_min = i - 1;
    }

    double d2_left = std::abs(v[k_min + 1] - 2.0 * v[k_min] + v[k_min - 1]);
    double d2_right = std::abs(v[k_min + 2] - 2.0 * v[k_min + 1] + v[k_min]);

    if (d2_left < d2_right) {
        k_min = k_min - 1;
    }
    int r = i - k_min;

    return r;
}

static double eno3_reconstruct_right_face(const std::vector<double>& v, int i) {
    int r = determine_stencil_shift(v, i);

    if (r == 0) {
        // Шаблон {i, i+1, i+2}
        return (1.0 / 3.0) * v[i] + (5.0 / 6.0) * v[i + 1] - (1.0 / 6.0) * v[i + 2];
    }
    else if (r == 1) {
        // Шаблон {i-1, i, i+1}
        return -(1.0 / 6.0) * v[i - 1] + (5.0 / 6.0) * v[i] + (1.0 / 3.0) * v[i + 1];
    }
    else { // r == 2
        // Шаблон {i-2, i-1, i}
        return (1.0 / 3.0) * v[i - 2] - (7.0 / 6.0) * v[i - 1] + (11.0 / 6.0) * v[i];
    }
}

static double eno3_reconstruct_left_face(const std::vector<double>& v, int i) {
    int r = determine_stencil_shift(v, i);

    if (r == 0) {
        // Шаблон {i, i+1, i+2} 
        return (11.0 / 6.0) * v[i] - (7.0 / 6.0) * v[i + 1] + (1.0 / 3.0) * v[i + 2];
    }
    else if (r == 1) {
        // Шаблон {i-1, i, i+1}
        return (1.0 / 3.0) * v[i - 1] + (5.0 / 6.0) * v[i] - (1.0 / 6.0) * v[i + 1];
    }
    else { // r == 2
        // Шаблон {i-2, i-1, i}
        return -(1.0 / 6.0) * v[i - 2] + (5.0 / 6.0) * v[i - 1] + (1.0 / 3.0) * v[i];
    }
}


void eno_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    std::vector<double> v1(total_cells), v2(total_cells), v3(total_cells);

    if (cfg.var_type == TypesOfVarForReconstruction::NONCONSERVATIVE) {
        for (int i = 0; i < total_cells; ++i) {
            v1[i] = grid.W[i].rho;
            v2[i] = grid.W[i].u;
            v3[i] = grid.W[i].p;
        }
    }
    else {
        for (int i = 0; i < total_cells; ++i) {
            v1[i] = grid.U[i].rho;
            v2[i] = grid.U[i].rhou;
            v3[i] = grid.U[i].E;
        }
    }

    for (int i = 0; i <= grid.Nx; ++i) {

        int idx_L = i + grid.num_fict - 1; 
        int idx_R = i + grid.num_fict;     

        State state_L, state_R;

        double val1_L = eno3_reconstruct_right_face(v1, idx_L);
        double val2_L = eno3_reconstruct_right_face(v2, idx_L);
        double val3_L = eno3_reconstruct_right_face(v3, idx_L);


        double val1_R = eno3_reconstruct_left_face(v1, idx_R);
        double val2_R = eno3_reconstruct_left_face(v2, idx_R);
        double val3_R = eno3_reconstruct_left_face(v3, idx_R);


        if (cfg.var_type == TypesOfVarForReconstruction::NONCONSERVATIVE) {
            state_L.rho = val1_L; state_L.u = val2_L; state_L.p = val3_L;
            state_R.rho = val1_R; state_R.u = val2_R; state_R.p = val3_R;
        }
        else {
            Conserved U_L = { val1_L, val2_L, val3_L };
            Conserved U_R = { val1_R, val2_R, val3_R };
            state_L = consToPhys(U_L, gamma);
            state_R = consToPhys(U_R, gamma);
        }

        fluxes[i] = compute_interface_flux(state_L, state_R, cfg);
    }
}