#pragma once
#include "types.h"
#include "grid.h"

// Метод FLIC (Fluid-in-Cell), Gentry, Martin, Daly (1966)
// Двухэтапная схема: Лагранжев шаг (давление) + Эйлеров шаг (конвекция)
// В 2D используется расщепление Стрэнга: X(dt/2) -> Y(dt) -> X(dt/2)
void flic_step(Grid& grid, double dt, const Config& cfg);
