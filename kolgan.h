#pragma once
#include "types.h"
#include <vector>

void kolgan_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes);
