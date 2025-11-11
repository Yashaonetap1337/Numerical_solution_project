#include "riemann_solver.h"
#include "euler_utils.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>


// определяем константы для итерационного решателя
constexpr int MAX_ITERATIONS = 20;
constexpr double GUESSP_TOLERANCE = 1e-6;


// функция давления f(p*, W_k) и ее производная f'(p*, W_k)
// они связывают давление p* с изменением скорости через волну.
// W_k - состояние (State) в k-той области (слева L или справа R)
// p_star - искомое давление в "звездной" области
static void evaluate_pressure_functions(const double p_star, const State& W_k, const double gamma,
    double& f, double& f_prime)
{
    const double a_k = soundSpeed(W_k, gamma); // скорость звука в области k

    if (p_star > W_k.p) { // УДАРНАЯ ВОЛНА 
        const double A = 2.0 / ((gamma + 1.0) * W_k.rho);
        const double B = (gamma - 1.0) / (gamma + 1.0) * W_k.p;
        const double p_sqrt = std::sqrt(A / (p_star + B));

        f = (p_star - W_k.p) * p_sqrt;
        f_prime = (1.0 - (p_star - W_k.p) / (2.0 * (B + p_star))) * p_sqrt;
    }
    else { //  ВОЛНА РАЗРЕЖЕНИЯ 
        const double p_ratio = p_star / W_k.p;
        const double p_pow = std::pow(p_ratio, (gamma - 1.0) / (2.0 * gamma));

        f = (2.0 * a_k / (gamma - 1.0)) * (p_pow - 1.0);
        f_prime = (1.0 / (W_k.rho * a_k)) * p_pow / p_ratio;
    }
}


// решатель Римана
State solve_general_riemann_problem(const State& W_L, const State& W_R, double xi, const Config& cfg, const ApproximationType approx_type) {
    const double gamma = cfg.phys.gamma;
    

    // проверка на вакуум (если плотность или давление нулевые)
    if (W_L.rho <= 0 || W_L.p <= 0 || W_R.rho <= 0 || W_R.p <= 0) {
        throw std::runtime_error("Vacuum or invalid initial state in Riemann solver.");
    }



    double p_star;
    double u_star;
    const double a_L = soundSpeed(W_L, gamma);
    const double a_R = soundSpeed(W_R, gamma);
    const double TOL = 1e-9;
    const double du = W_R.u - W_L.u;



    // начальное приближение для давления
    switch (approx_type) {
    case ApproximationType::RAREFACTION: {
        // приближение двух волн разрежения (TRRS)
        const double p_pow_L = std::pow(W_L.p, (gamma - 1.0) / (2.0 * gamma));
        const double p_pow_R = std::pow(W_R.p, (gamma - 1.0) / (2.0 * gamma));
        const double numerator = a_L + a_R - ((gamma - 1.0) / 2.0) * (W_R.u - W_L.u);
        const double denominator = a_L / p_pow_L + a_R / p_pow_R;
        p_star = std::pow(numerator / denominator, 2.0 * gamma / (gamma - 1.0));
        break;
    }
    case ApproximationType::SHOCK: {
        // приближение двух ударных волн (TSRS) 
        const double A_L = 2.0 / ((gamma + 1.0) * W_L.rho);
        const double B_L = (gamma - 1.0) / (gamma + 1.0) * W_L.p;
        const double A_R = 2.0 / ((gamma + 1.0) * W_R.rho);
        const double B_R = (gamma - 1.0) / (gamma + 1.0) * W_R.p;

        const double g = (A_L / (A_L + A_R)) * (W_R.p - W_L.p)
            + (A_R / (A_L + A_R)) * (W_R.u - W_L.u) * (a_L + a_R);
        p_star = std::max(1e-9, 0.5 * (W_L.p + W_R.p) + 0.5 * g);
        break;
    }
    case ApproximationType::ACOUSTIC: {
        // акустическое приближение (PVRS)
        p_star = 0.5 * (W_L.p + W_R.p) - 0.125 * (W_R.u - W_L.u) * (W_L.rho + W_R.rho) * (a_L + a_R);
        break;
    }
    case ApproximationType::AVERAGE: {
        // среднее арифметическое
        p_star = 0.5 * (W_L.p + W_R.p);
        break;
    }
    default:
        throw std::runtime_error("Unknown approximation type in Riemann solver");
    }

    // защита от отрицательных давлений
    p_star = std::max(1e-9, p_star);

    // итерационный метод Ньютона
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double f_L, f_prime_L, f_R, f_prime_R;
        evaluate_pressure_functions(p_star, W_L, gamma, f_L, f_prime_L);
        evaluate_pressure_functions(p_star, W_R, gamma, f_R, f_prime_R);

        const double F = f_L + f_R + (W_R.u - W_L.u);
        const double F_prime = f_prime_L + f_prime_R;
        const double p_new = p_star - F / F_prime;

        if (std::abs(p_new - p_star) / (0.5 * (p_new + p_star)) < GUESSP_TOLERANCE) {
            break;
        }
        p_star = std::max(1e-9, p_new);
    }

    double f_L, f_prime_L, f_R, f_prime_R;
    evaluate_pressure_functions(p_star, W_L, gamma, f_L, f_prime_L);
    evaluate_pressure_functions(p_star, W_R, gamma, f_R, f_prime_R);
    u_star = 0.5 * (W_L.u + W_R.u) + 0.5 * (f_R - f_L);



    // отбор решения
    State solution_state;
    if (xi <= u_star) { // точка (0,t) находится слева от контактного разрыва 
        if (p_star <= W_L.p) { // левая волна - РАЗРЕЖЕНИЕ
            const double S_HL = W_L.u - a_L; // скорость головы волны разрежения
            if (xi <= S_HL) { // точка находится в невозмущенной области слева
                solution_state = W_L;
            }
            else {
                const double a_star_L = a_L * std::pow(p_star / W_L.p, (gamma - 1.0) / (2.0 * gamma));
                const double S_TL = u_star - a_star_L; // скорость хвоста волны разрежения
                if (xi > S_TL) { // точка находится в "звездной" области
                    solution_state.rho = W_L.rho * std::pow(p_star / W_L.p, 1.0 / gamma);
                    solution_state.u = u_star;
                    solution_state.p = p_star;
                }
                else { // точка находится ВНУТРИ веера волны разрежения
                    const double C1 = 2.0 / (gamma + 1.0);
                    const double C2 = (gamma - 1.0) / 2.0;
                    solution_state.u = C1 * (a_L + C2 * W_L.u + xi);
                    const double a_fan = C1 * (a_L + C2 * (W_L.u - xi));
                    solution_state.rho = W_L.rho * std::pow(a_fan / a_L, 2.0 / (gamma - 1.0));
                    solution_state.p = W_L.p * std::pow(a_fan / a_L, 2.0 * gamma / (gamma - 1.0));
                }
            }
        }
        else { // левая волна - УДАРНАЯ
            const double S_L = W_L.u - a_L * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / W_L.p) + (gamma - 1.0) / (2.0 * gamma));
            if (xi <= S_L) { // точка перед ударной волной
                solution_state = W_L;
            }
            else { // точка за ударной волной
                solution_state.rho = W_L.rho * (p_star / W_L.p + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * p_star / W_L.p + 1.0);
                solution_state.u = u_star;
                solution_state.p = p_star;
            }
        }
    }
    else { // точка находится справа от контактного разрыва
        if (p_star <= W_R.p) { // правая волна - РАЗРЕЖЕНИЕ
            const double S_HR = W_R.u + a_R; // скорость "головы"
            if (xi >= S_HR) { // точка в невозмущенной области справа
                solution_state = W_R;
            }
            else {
                const double a_star_R = a_R * std::pow(p_star / W_R.p, (gamma - 1.0) / (2.0 * gamma));
                const double S_TR = u_star + a_star_R; // скорость "хвоста"
                if (xi < S_TR) { // точка в "звездной" области
                    solution_state.rho = W_R.rho * std::pow(p_star / W_R.p, 1.0 / gamma);
                    solution_state.u = u_star;
                    solution_state.p = p_star;
                }
                else { // точка ВНУТРИ веера волны разрежения
                    const double C1 = 2.0 / (gamma + 1.0);
                    const double C2 = (gamma - 1.0) / 2.0;
                    solution_state.u = C1 * (-a_R + C2 * W_R.u + xi);
                    const double a_fan = C1 * (a_R - C2 * (W_R.u - xi));
                    solution_state.rho = W_R.rho * std::pow(a_fan / a_R, 2.0 / (gamma - 1.0));
                    solution_state.p = W_R.p * std::pow(a_fan / a_R, 2.0 * gamma / (gamma - 1.0));
                }
            }
        }
        else { // правая волна - УДАРНАЯ
            const double S_R = W_R.u + a_R * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / W_R.p) + (gamma - 1.0) / (2.0 * gamma));
            if (xi >= S_R) { // Точка перед ударной волной (справа)
                solution_state = W_R;
            }
            else { // Точка за ударной волной
                solution_state.rho = W_R.rho * (p_star / W_R.p + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * p_star / W_R.p + 1.0);
                solution_state.u = u_star;
                solution_state.p = p_star;
            }
        }
    }

    return solution_state;
}
