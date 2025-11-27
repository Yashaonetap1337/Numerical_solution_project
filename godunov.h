#pragma once
#include "types.h"

void godunov_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes);