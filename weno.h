#pragma once
#include "types.h"
#include <vector>

void weno_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes);
