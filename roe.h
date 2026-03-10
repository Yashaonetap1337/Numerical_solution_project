#pragma once
#include "types.h"

Flux solve_roe_flux(const State& W_L, const State& W_R, const Config& cfg, char dir);