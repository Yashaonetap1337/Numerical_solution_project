#pragma once
#include "types.h"
#include <string>

void run_simulation(const Config& cfg, const std::string& outputFilename);
void generate_analytical_snapshot(const Config& cfg, double t, const std::string& outputFilename);