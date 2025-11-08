#include "godunov.h"
#include "euler_utils.h" 
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>



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

void acoustic_solver(Grid& grid, const Config& cfg)
{
    State W_L = cfg.phys.left;
    State W_R = cfg.phys.right;
    const double gamma = cfg.phys.gamma;

    const double c_L = soundSpeed(W_L, gamma);
    const double c_R = soundSpeed(W_R, gamma);

    std::cout << "=== ACOUSTIC SOLVER DEBUG ===" << std::endl;
    std::cout << "Left state: rho=" << W_L.rho << " u=" << W_L.u << " p=" << W_L.p << " c=" << c_L << std::endl;
    std::cout << "Right state: rho=" << W_R.rho << " u=" << W_R.u << " p=" << W_R.p << " c=" << c_R << std::endl;

    // Вычисляем p* и u*
    double numerator_p = W_R.rho * W_L.rho * c_R * c_L * ((W_L.u - W_R.u) + W_R.p / (W_R.rho * c_R) + W_L.p / (W_L.rho * c_L));
    double denominator_p = W_L.rho * c_L + W_R.rho * c_R;
    double p_star = std::max(1e-9, numerator_p / denominator_p);

    double numerator_u = (W_L.u * W_L.rho * c_L - W_R.u * W_R.rho * c_R) + (W_L.p - W_R.p);
    double denominator_u = W_L.rho * c_L + W_R.rho * c_R;
    double u_star = numerator_u / denominator_u;

    // Плотности в звездных областях
    double rho_star_L = -(W_L.p - p_star) / (c_L * c_L) + W_L.rho;
    double rho_star_R = -(W_R.p - p_star) / (c_R * c_R) + W_R.rho;

    // ИСПРАВЛЕННЫЕ ГРАНИЦЫ ВОЛН - учитываем скорость газа!
    const double x_left_wave = (W_L.u - c_L) * cfg.grid.t_final;
    const double x_contact = u_star * cfg.grid.t_final;
    const double x_right_wave = (W_R.u + c_R) * cfg.grid.t_final;

    std::cout << "CORRECTED Wave boundaries: left=" << x_left_wave
        << " contact=" << x_contact << " right=" << x_right_wave << std::endl;

    // Проверим порядок границ
    if (x_left_wave > x_right_wave) {
        std::cout << "WARNING: Wave boundaries are inverted! This indicates strong shock conditions." << std::endl;
    }

    int region1_count = 0, region2_count = 0, region3_count = 0, region4_count = 0;

    for (int i = 0; i < grid.W.size(); ++i) {
        double x = grid.x_centers[i];

        // ИСПРАВЛЕННАЯ ЛОГИКА ОБЛАСТЕЙ
        if (x <= x_left_wave) {
            grid.W[i] = W_L;
            region1_count++;
        }
        else if (x <= x_contact) {
            grid.W[i].rho = rho_star_L;
            grid.W[i].u = u_star;
            grid.W[i].p = p_star;
            region2_count++;
        }
        else if (x <= x_right_wave) {
            grid.W[i].rho = rho_star_R;
            grid.W[i].u = u_star;
            grid.W[i].p = p_star;
            region3_count++;
        }
        else {
            grid.W[i] = W_R;
            region4_count++;
        }
    }

    std::cout << "--- CORRECTED Region statistics ---" << std::endl;
    std::cout << "Region 1 (left): " << region1_count << " cells" << std::endl;
    std::cout << "Region 2 (star left): " << region2_count << " cells" << std::endl;
    std::cout << "Region 3 (star right): " << region3_count << " cells" << std::endl;
    std::cout << "Region 4 (right): " << region4_count << " cells" << std::endl;

    // Дополнительная проверка для проблемной зоны
    std::cout << "--- Checking transition areas ---" << std::endl;
    for (int i = 0; i < grid.W.size(); ++i) {
        double x = grid.x_centers[i];
        if (std::abs(x - x_contact) < 0.01 || std::abs(x - x_right_wave) < 0.01) {
            std::cout << "TRANSITION: x=" << x << " rho=" << grid.W[i].rho
                << " (contact at " << x_contact << ", right wave at " << x_right_wave << ")" << std::endl;
        }
    }
}


// решатель Римана
static State solve_exact_riemann_problem(const State& W_L, const State& W_R, const double gamma, const ApproximationType approx_type) {
    // итерационный поиск давления p* и скорости u*
    double p_star;
    double u_star;
    const double a_L = soundSpeed(W_L, gamma);
    const double a_R = soundSpeed(W_R, gamma);
    const double TOL = 1e-9;//вот тут надо будет побаловаться чтобы чекнуть чзх втф
    const double du = W_R.u - W_L.u; //нахуя? хз



   // используем приближение для двух волн разрежения 
    switch (approx_type) {
        case ApproximationType::RAREFACTION: {
            // Приближение двух волн разрежения
            const double p_pow_L = std::pow(W_L.p, (gamma - 1.0) / (2.0 * gamma));
            const double p_pow_R = std::pow(W_R.p, (gamma - 1.0) / (2.0 * gamma));

            const double numerator = a_L + a_R - ((gamma - 1.0) / 2.0) * (W_R.u - W_L.u);
            const double denominator = a_L / p_pow_L + a_R / p_pow_R;

            p_star = std::pow(numerator / denominator, 2.0 * gamma / (gamma - 1.0));
            break;
        }

        case ApproximationType::SHOCK: {
            // Приближение двух ударных волн
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
            // Акустическое приближение
            p_star = 0.5 * (W_L.p + W_R.p)
                - 0.125 * (W_R.u - W_L.u) * (W_L.rho + W_R.rho) * (a_L + a_R);
            break;
        }

        case ApproximationType::AVERAGE: {
            // Среднее арифметическое
            p_star = 0.5 * (W_L.p + W_R.p);
            break;
        }

        default:
            throw std::runtime_error("Unknown approximation type in Riemann solver");
        }

    // защита от отрицательных давлений
    p_star = std::max(1e-9, p_star);

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
    if (0.0 <= u_star) { // точка (0,t) находится слева от контактного разрыва 
        if (p_star <= W_L.p) { // левая волна - РАЗРЕЖЕНИЕ
            const double S_HL = W_L.u - a_L; // скорость головы волны разрежения
            if (0.0 <= S_HL) { // точка находится в невозмущенной области слева
                solution_state = W_L;
            }
            else {
                const double a_star_L = a_L * std::pow(p_star / W_L.p, (gamma - 1.0) / (2.0 * gamma));
                const double S_TL = u_star - a_star_L; // скорость хвоста волны разрежения
                if (0.0 > S_TL) { // точка находится в "звездной" области
                    solution_state.rho = W_L.rho * std::pow(p_star / W_L.p, 1.0 / gamma);
                    solution_state.u = u_star;
                    solution_state.p = p_star;
                }
                else { // точка находится ВНУТРИ веера волны разрежения
                    const double C = 2.0 / (gamma + 1.0);
                    const double u_fan = C * (a_L + (gamma - 1.0) / 2.0 * W_L.u);
                    const double a_fan = C * (a_L + (gamma - 1.0) / 2.0 * W_L.u);
                    solution_state.u = u_fan;
                    solution_state.rho = W_L.rho * std::pow(a_fan / a_L, 2.0 / (gamma - 1.0));
                    solution_state.p = W_L.p * std::pow(a_fan / a_L, 2.0 * gamma / (gamma - 1.0));
                }
            }
        }
        else { // левая волна - УДАРНАЯ
            const double S_L = W_L.u - a_L * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / W_L.p) + (gamma - 1.0) / (2.0 * gamma));
            if (0.0 <= S_L) { // точка перед ударной волной
                solution_state = W_L;
            }
            else { // точка за ударной волной
                solution_state.rho = W_L.rho * (p_star / W_L.p + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * p_star / W_L.p + 1.0);
                solution_state.u = u_star;
                solution_state.p = p_star;
            }
        }
    }
    else { // точка (0,t) находится справа от контактного разрыва
        if (p_star <= W_R.p) { // правая волна - РАЗРЕЖЕНИЕ
            const double S_HR = W_R.u + a_R; // скорость "головы"
            if (0.0 >= S_HR) { // точка в невозмущенной области справа
                solution_state = W_R;
            }
            else {
                const double a_star_R = a_R * std::pow(p_star / W_R.p, (gamma - 1.0) / (2.0 * gamma));
                const double S_TR = u_star + a_star_R; // скорость "хвоста"
                if (0.0 < S_TR) { // точка в "звездной" области
                    solution_state.rho = W_R.rho * std::pow(p_star / W_R.p, 1.0 / gamma);
                    solution_state.u = u_star;
                    solution_state.p = p_star;
                }
                else { // точка ВНУТРИ веера волны разрежения
                    const double C1 = 2.0 / (gamma + 1.0);
                    const double C2 = (gamma - 1.0) / (gamma + 1.0);
                    const double u_fan = C1 * (-a_R + (gamma - 1.0) / 2.0 * W_R.u);
                    const double a_fan = C1 * (a_R - (gamma - 1.0) / 2.0 * W_R.u);
                    solution_state.u = u_fan;
                    solution_state.rho = W_R.rho * std::pow(a_fan / a_R, 2.0 / (gamma - 1.0));
                    solution_state.p = W_R.p * std::pow(a_fan / a_R, 2.0 * gamma / (gamma - 1.0));
                }
            }
        }
        else { // правая волна - УДАРНАЯ
            const double S_R = W_R.u + a_R * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / W_R.p) + (gamma - 1.0) / (2.0 * gamma));
            if (0.0 >= S_R) { // точка за ударной волной
                solution_state.rho = W_R.rho * (p_star / W_R.p + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * p_star / W_R.p + 1.0);
                solution_state.u = u_star;
                solution_state.p = p_star;
            }
            else { // точка перед ударной волной
                solution_state = W_R;
            }
        }
    }

    return solution_state;
}


// функция, выполняющая один шаг по времени методом Годунова
void godunov_step(Grid& grid, double dt, const Config& cfg) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;
    const auto approx_type = cfg.approx_type;
    // вектор для хранения потоков на границах ячеек.
    std::vector<Flux> fluxes(grid.Nx + 1);

    for (int i = 0; i <= grid.Nx; ++i) {
        // определяем состояния слева и справа от границы i
        const State W_L = grid.W[i + grid.num_fict - 1];
        const State W_R = grid.W[i + grid.num_fict];

        // решаем задачу Римана для этой границы
        const State state_at_interface = solve_exact_riemann_problem(W_L, W_R, gamma, approx_type);


        // вычисляем поток по найденному состоянию
        fluxes[i] = physToFlux(state_at_interface, gamma);
    }

    // цикл по ячейкам для обновления решения
    for (int i = 0; i < grid.Nx; ++i) {
        
        const int cell_idx = i + grid.num_fict; // индекс реальной ячейки в общем массиве

        
        const Flux& F_left = fluxes[i];
        const Flux& F_right = fluxes[i + 1];

        // конечно-объемная формула обновления
        grid.U[cell_idx].rho -= (dt / dx) * (F_right.rho_f - F_left.rho_f);
        grid.U[cell_idx].rhou -= (dt / dx) * (F_right.rhou_f - F_left.rhou_f);
        grid.U[cell_idx].E -= (dt / dx) * (F_right.E_f - F_left.E_f);
    }
}