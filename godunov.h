#pragma once
#include "types.h"

void godunov_step(Grid& grid, double dt, const Config& cfg);
void acoustic_solver(Grid& grid, const Config& cfg);