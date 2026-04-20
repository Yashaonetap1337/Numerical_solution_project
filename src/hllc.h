#pragma once
#include "types.h"

enum class WaveSpeedMethod { ISOENTROPIC, LINEARIZED, HYBRID };

Flux hllc_flux(const State& W_L, const State& W_R, double gamma, char dir,
    WaveSpeedMethod method = WaveSpeedMethod::HYBRID);