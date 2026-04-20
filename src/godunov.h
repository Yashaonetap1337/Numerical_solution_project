#pragma once
#include "types.h"

void godunov_flux_computation(const Grid& grid, const Config& cfg, std::vector<std::vector<Flux>>& fluxes_x, std::vector<std::vector<Flux>>& fluxes_y);