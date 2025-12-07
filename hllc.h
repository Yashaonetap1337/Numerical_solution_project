#pragma once

#include "types.h"
#include <vector>
enum class WaveSpeedMethod {
    ISOENTROPIC,    ///< Изоэнтропическая оценка (обе волны - волны разрежения)
    LINEARIZED,     ///< Линеаризованная оценка
    HYBRID          ///< Гибридная оценка (лучшая точность)
};
void hllc_wave_speeds(const State& W_L, const State& W_R, double gamma,
    WaveSpeedMethod method, double& S_L, double& S_R, double& S_M);
static Flux hllc_flux_with_method(const State& W_L, const State& W_R, double gamma,
    WaveSpeedMethod method);
void hllc_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes);
