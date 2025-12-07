#include "hll.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

void hll_davis_wave_speeds(const State& W_L, const State& W_R, double gamma, double& S_L, double& S_R) {
    // Базовая проверка на валидность состояний
    if (W_L.rho <= 0 || W_R.rho <= 0 || W_L.p <= 0 || W_R.p <= 0) {
        throw std::runtime_error("Invalid state in HLL wave speeds calculation");
    }

    double a_L = soundSpeed(W_L, gamma);
    double a_R = soundSpeed(W_R, gamma);

    // Оценки через простые волны (разрежения)
    double u_star_est = 0.5 * (W_L.u + W_R.u) + 0.5 * (a_L - a_R);
    double a_star_est = 0.5 * (a_L + a_R) + 0.25 * (gamma - 1.0) * (W_L.u - W_R.u);

    S_L = u_star_est - a_star_est;
    S_R = u_star_est + a_star_est;

    // Расширяем интервал для гарантии включения всех волн (включая ударные)
    double S_L_min = std::min(W_L.u - a_L, W_R.u - a_R);
    double S_R_max = std::max(W_L.u + a_L, W_R.u + a_R);

    S_L = std::min(S_L, S_L_min);
    S_R = std::max(S_R, S_R_max);

    // Защита от вырожденного случая
    const double eps = 1e-8;
    if (S_L >= S_R) {
        S_L = std::min(S_L_min, 0.0) - eps;
        S_R = std::max(S_R_max, 0.0) + eps;
    }
}

Flux hll_flux(const State& W_L, const State& W_R, double gamma) {
    // Вычисляем волновые скорости по Дэвису
    double S_L, S_R;
    hll_davis_wave_speeds(W_L, W_R, gamma, S_L, S_R);

    // Вычисляем консервативные переменные
    Conserved U_L = physToCons(W_L, gamma);
    Conserved U_R = physToCons(W_R, gamma);

    // Вычисляем физические потоки
    Flux F_L = physToFlux(W_L, gamma);
    Flux F_R = physToFlux(W_R, gamma);

    Flux F_hll;

    // Формула HLL потока
    if (S_L >= 0.0) {
        // Вся информация движется вправо
        F_hll = F_L;
    }
    else if (S_R <= 0.0) {
        // Вся информация движется влево
        F_hll = F_R;
    }
    else {
        // Зона взаимодействия: используем HLL приближение
        double denom = S_R - S_L;

        // Защита от деления на ноль (маловероятно после предыдущих проверок)
        if (std::abs(denom) < 1e-12) {
            // В вырожденном случае используем среднее
            F_hll.rho_f = 0.5 * (F_L.rho_f + F_R.rho_f);
            F_hll.rhou_f = 0.5 * (F_L.rhou_f + F_R.rhou_f);
            F_hll.E_f = 0.5 * (F_L.E_f + F_R.E_f);
        }
        else {
            double factor = 1.0 / denom;

            // Формула HLL потока: F_hll = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)
            F_hll.rho_f = factor * (S_R * F_L.rho_f - S_L * F_R.rho_f
                + S_L * S_R * (U_R.rho - U_L.rho));

            F_hll.rhou_f = factor * (S_R * F_L.rhou_f - S_L * F_R.rhou_f
                + S_L * S_R * (U_R.rhou - U_L.rhou));

            F_hll.E_f = factor * (S_R * F_L.E_f - S_L * F_R.E_f
                + S_L * S_R * (U_R.E - U_L.E));
        }
    }

    // Защита от некорректных значений
    if (!std::isfinite(F_hll.rho_f) || !std::isfinite(F_hll.rhou_f) || !std::isfinite(F_hll.E_f)) {
        // В случае проблем возвращаем простой upwind
        if (0.5 * (W_L.u + W_R.u) >= 0.0) {
            F_hll = F_L;
        }
        else {
            F_hll = F_R;
        }
    }

    return F_hll;
}

void hll_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const int num_faces = grid.Nx + 1;
    fluxes.resize(num_faces);

    const double gamma = cfg.phys.gamma;

    // Проходим по всем граням
    for (int i = 0; i < num_faces; ++i) {
        // Индексы ячеек слева и справа от грани
        int left_idx = i + grid.num_fict - 1;
        int right_idx = i + grid.num_fict;

        // Получаем состояния на грани (без реконструкции - метод первого порядка)
        const State& W_L = grid.W[left_idx];
        const State& W_R = grid.W[right_idx];

        // Вычисляем HLL поток
        fluxes[i] = hll_flux(W_L, W_R, gamma);
    }
}