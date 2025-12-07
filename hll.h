#pragma once

#include "types.h"
#include <vector>
void hll_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes);

void hll_davis_wave_speeds(const State& W_L, const State& W_R, double gamma, double& S_L, double& S_R);

Flux hll_flux(const State& W_L, const State& W_R, double gamma);