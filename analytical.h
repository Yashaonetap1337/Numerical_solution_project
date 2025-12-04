#pragma once
#include "types.h"
#include <string>


void generate_analytical_solution(const Config& cfg, double t, const std::string& outputFilename);
void save_analytical_solution(const Config& cfg, const std::string& outputFilename);
void generate_analytical_snapshot(const Config& cfg, double t, const std::string& filename);