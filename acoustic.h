#pragma once
#include "types.h"

// ‘ункци€, выполн€юща€ один шаг по времени акустическим методом
void acoustic_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes);