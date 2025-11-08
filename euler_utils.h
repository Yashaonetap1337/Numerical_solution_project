#pragma once
#include "types.h"

// перевод из физических в консервативные
Conserved physToCons(const State& W, double gamma);

// перевод из консервативных в физические
State consToPhys(const Conserved& U, double gamma);

// вычисление потока по физическим переменным
Flux physToFlux(const State& W, double gamma);

// вычисление скорости звука
double soundSpeed(const State& W, double gamma);