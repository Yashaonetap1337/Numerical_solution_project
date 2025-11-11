#pragma once
#include "types.h"
#include <vector>

void kolgan_step(Grid& grid, double dt, const Config& cfg, std::vector<Flux>& fluxes);
