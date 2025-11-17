#pragma once
#include "types.h"
#include <string>

// √енерирует аналитическое решение дл€ заданного времени
void generate_analytical_solution(const Config& cfg, double t, const std::string& outputFilename);

// —охран€ет аналитическое решение дл€ финального времени
void save_analytical_solution(const Config& cfg, const std::string& outputFilename);

// √енерирует аналитический снимок дл€ момента времени (аналог численного снимка)
void generate_analytical_snapshot(const Config& cfg, double t, const std::string& filename);