#pragma once
#include "types.h"

State solve_general_riemann_problem(const State& W_L, const State& W_R, double xi, const Config& cfg, const ApproximationType approx_type);

