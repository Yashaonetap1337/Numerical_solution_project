#include "types.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <iostream>

// Для 1-й волны (u - c)
State calc_sonic_state_1(const State& W_start, double gamma) {
    State W_sonic;
    double c_start = std::sqrt(gamma * W_start.p / W_start.rho);

    double k = (gamma - 1.0) / (gamma + 1.0);
    double c_sonic = k * (W_start.u + 2.0 * c_start / (gamma - 1.0));

    W_sonic.u = c_sonic;
    W_sonic.rho = W_start.rho * std::pow(c_sonic / c_start, 2.0 / (gamma - 1.0));
    W_sonic.p = W_start.p * std::pow(c_sonic / c_start, 2.0 * gamma / (gamma - 1.0));

    return W_sonic;
}

// Для 3-й волны (u + c)
State calc_sonic_state_3(const State& W_start, double gamma) {
    State W_sonic;
    double c_start = std::sqrt(gamma * W_start.p / W_start.rho);

    double k = (gamma - 1.0) / (gamma + 1.0);
    double J_minus = W_start.u - 2.0 * c_start / (gamma - 1.0);
    double c_sonic = -J_minus * k;

    W_sonic.u = -c_sonic;
    W_sonic.rho = W_start.rho * std::pow(c_sonic / c_start, 2.0 / (gamma - 1.0));
    W_sonic.p = W_start.p * std::pow(c_sonic / c_start, 2.0 * gamma / (gamma - 1.0));

    return W_sonic;
}

Flux solve_osher_flux(const State& W_L, const State& W_R, const Config& cfg) {
    const double gamma = cfg.phys.gamma;
    const double gm1 = gamma - 1.0;
    const double gp1 = gamma + 1.0;
    const double g_pow = 2.0 * gamma / gm1;

    double c_L = std::sqrt(gamma * W_L.p / W_L.rho);
    double c_R = std::sqrt(gamma * W_R.p / W_R.rho);

    // Находим промежуточные состояния W_1/3 и W_2/3
    double H_L = W_L.u + 2.0 * c_L / gm1;
    double H_R = W_R.u - 2.0 * c_R / gm1;

    double z = gm1 / (2.0 * gamma);
    double K_L = std::pow(W_L.p, z) / c_L; 
    double K_R = std::pow(W_R.p, z) / c_R;


    double u_star = (K_L * H_L + K_R * H_R) / (K_L + K_R);


    double c_L_star = c_L + 0.5 * gm1 * (W_L.u - u_star);
    double c_R_star = c_R + 0.5 * gm1 * (u_star - W_R.u);


    State W_13, W_23;
    W_13.u = u_star;
    W_23.u = u_star;


    W_13.rho = W_L.rho * std::pow(c_L_star / c_L, 2.0 / gm1);
    W_13.p = W_L.p * std::pow(c_L_star / c_L, g_pow);

    W_23.rho = W_R.rho * std::pow(c_R_star / c_R, 2.0 / gm1);
    W_23.p = W_R.p * std::pow(c_R_star / c_R, g_pow);


    Flux F_total = physToFlux(W_L, gamma);

    // ПУТЬ 1: L -> 1/3 (Волна u-c)
    double lambda_start = W_L.u - c_L;
    double lambda_end = W_13.u - c_L_star;

    if (lambda_start >= 0 && lambda_end >= 0) {
        // Все положительно, вклад 0
    }
    else if (lambda_start <= 0 && lambda_end <= 0) {
        // Все отрицательно, добавляем F(1/3) - F(L)
        F_total = F_total + (physToFlux(W_13, gamma) - physToFlux(W_L, gamma));
    }
    else {
        State W_sonic = calc_sonic_state_1(W_L, gamma);
        Flux F_sonic = physToFlux(W_sonic, gamma);

        if (lambda_start < 0) { // Разрежение 
            // Интеграл от L до Sonic
            F_total = F_total + (F_sonic - physToFlux(W_L, gamma));
        }
        else { // Сжатие
            // Интеграл от Sonic до 1/3
            F_total = F_total + (physToFlux(W_13, gamma) - F_sonic);
        }
    }

    // ПУТЬ 2: 1/3 -> 2/3 (Контактный разрыв, u)
    if (u_star < 0) {
        F_total = F_total + (physToFlux(W_23, gamma) - physToFlux(W_13, gamma));
    }

    // ПУТЬ 3: 2/3 -> R (Волна u+c)
    lambda_start = W_23.u + c_R_star;
    double lambda_end_R = W_R.u + c_R;

    if (lambda_start >= 0 && lambda_end_R >= 0) {
        // Вклад 0
    }
    else if (lambda_start <= 0 && lambda_end_R <= 0) {
        // Вклад F(R) - F(2/3)
        F_total = F_total + (physToFlux(W_R, gamma) - physToFlux(W_23, gamma));
    }
    else {
        // Звуковая точка
        State W_sonic = calc_sonic_state_3(W_R, gamma);
        Flux F_sonic = physToFlux(W_sonic, gamma);

        if (lambda_start < 0) {
            F_total = F_total + (F_sonic - physToFlux(W_23, gamma));
        }
        else {
            F_total = F_total + (physToFlux(W_R, gamma) - F_sonic);
        }
    }
    // -- - ФИНАЛЬНАЯ ЗАЩИТА ОТ NAN-- -
        // Если Ошер совсем сломался, возвращаем поток Лакса-Фридрихса (Русанова) или просто поток стены
        if (std::isnan(F_total.rho_f) || std::isnan(F_total.rhou_f) || std::isnan(F_total.E_f) ||
            std::isinf(F_total.rho_f) || std::isinf(F_total.rhou_f) || std::isinf(F_total.E_f))
        {
            // Аварийный вариант: возвращаем простой средний поток (это грубо, но не крашнет программу)
            Flux F_L_phys = physToFlux(W_L, gamma);
            Flux F_R_phys = physToFlux(W_R, gamma);

            F_total.rho_f = 0.5 * (F_L_phys.rho_f + F_R_phys.rho_f);
            F_total.rhou_f = 0.5 * (F_L_phys.rhou_f + F_R_phys.rhou_f);
            F_total.E_f = 0.5 * (F_L_phys.E_f + F_R_phys.E_f);
        }

    return F_total;
}