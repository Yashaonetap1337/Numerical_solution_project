#pragma once
#include "types.h"

// omega[i][j] — mass fraction of unreacted explosive (0..1).
// Carried alongside the grid, updated in-place each time step.
using OmegaField = std::vector<std::vector<double>>;

// Allocate omega field sized to match the grid (including fictitious cells).
OmegaField make_omega_field(const Grid& grid);

// One full Mader time step.
// W      — primitive variables at t   (input)
// W_new  — primitive variables at t+dt (output)
// omega  — reaction progress field, updated in-place
void MaderTimeStep(
    const Grid&   grid,
    const Config& cfg,
    const std::vector<std::vector<State>>& W,
          std::vector<std::vector<State>>& W_new,
    OmegaField&   omega,
    double dt
);
