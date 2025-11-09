#pragma once
#include "types.h"

// Функция, выполняющая один шаг по времени акустическим методом
void acoustic_step(Grid& grid, double dt, const Config& cfg);
