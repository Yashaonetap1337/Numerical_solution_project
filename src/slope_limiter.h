#pragma once
#include "types.h"
#include <vector>

std::pair<double, double> reconstruct_slope_1d(const std::vector<double>& q, int i, LimiterType type);