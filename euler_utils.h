#pragma once
#include "types.h"

Conserved physToCons(const State& W, double gamma);
State consToPhys(const Conserved& U, double gamma);
Flux fluxX(const State& W, double gamma);
Flux fluxY(const State& W, double gamma);
double soundSpeed(const State& W, double gamma);