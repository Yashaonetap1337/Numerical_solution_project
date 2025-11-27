#include "eno.h"
#include "choice_of_riemann_solvers.h" 
#include "euler_utils.h"
#include <vector>
#include <cmath>
#include <algorithm>

// ¬спомогательна€ функци€ дл€ ENO-реконструкции ќƒЌќ… переменной
// на границе i+1/2, использу€ трехточечный шаблон {i-1, i, i+1}
static void eno3_reconstruction(const std::vector<double>& var, int i,
    double& val_L, double& val_R)
{
    // «начени€ в €чейках, участвующих в реконструкции
    const double v_m1 = var[i - 1]; // f_{i-1}
    const double v_0 = var[i];   // f_i
    const double v_p1 = var[i + 1]; // f_{i+1}

    // –еконструкци€ дл€ W_L (значение слева от границы i+1/2) 
    // Ёто значение строитс€ внутри €чейки i, использу€ соседей

    // ¬ычисл€ем гладкость на двух возможных шаблонах дл€ €чейки i:
    // Ўаблон {i-1, i} и {i, i+1}
    double smoothness_left = std::abs(v_0 - v_m1);
    double smoothness_right = std::abs(v_p1 - v_0);


    // Ќаклон, вычисленный по левому соседу
    double delta_L = v_0 - v_m1;
    // Ќаклон, вычисленный по правому соседу
    double delta_R = v_p1 - v_0;

    // »спользуем простой ограничитель, который выбирает более гладкий наклон
    double delta = (smoothness_left < smoothness_right) ? delta_L : delta_R;

    // Ёкстраполируем значение из центра €чейки `i` на правую границу `i+1/2`
    val_L = v_0 + 0.5 * delta;

    // –еконструкци€ дл€ W_R (значение справа от границы i+1/2) 
    // Ёто значение строитс€ внутри €чейки i+1, использу€ ее соседей

    // Ўаблоны дл€ €чейки i+1: {i, i+1} и {i+1, i+2}
    const double v_p2 = var[i + 2];
    smoothness_left = std::abs(v_p1 - v_0);
    smoothness_right = std::abs(v_p2 - v_p1);

    delta_L = v_p1 - v_0;
    delta_R = v_p2 - v_p1;

    delta = (smoothness_left < smoothness_right) ? delta_L : delta_R;

    // Ёкстраполируем значение из центра €чейки `i+1` на левую границу `i+1/2`
    val_R = v_p1 - 0.5 * delta;
}


void eno_flux_computation(const Grid& grid, const Config& cfg, std::vector<Flux>& fluxes) {
    const double dx = grid.dx;
    const double gamma = cfg.phys.gamma;
    const int total_cells = grid.Nx + 2 * grid.num_fict;

    // »звлекаем скал€рные переменные в отдельные векторы
    std::vector<double> rho(total_cells), u(total_cells), p(total_cells);
    for (int i = 0; i < total_cells; ++i) {
        rho[i] = grid.W[i].rho; u[i] = grid.W[i].u; p[i] = grid.W[i].p;
    }

    for (int i = 0; i <= grid.Nx; ++i) {
        // √лобальный индекс €чейки слева от границы
        int cell_i = i + grid.num_fict - 1;

        // ƒл€ трехточечного шаблона ENO нам нужен один сосед слева и два справа
        // от `cell_i` дл€ полной реконструкции на границе `i+1/2`.
        // ”бедимс€, что мы не выходим за границы массива.
        if (cell_i < 1 || cell_i > total_cells - 3) {
            // Ќа самых кра€х сетки используем реконструкцию 1-го пор€дка (√одунов)
            const State W_L = grid.W[cell_i];
            const State W_R = grid.W[cell_i + 1];
            const State state_at_interface = solve_riemann_problem(W_L, W_R, 0.0, cfg, cfg.approx_type);
            fluxes[i] = physToFlux(state_at_interface, gamma);
            continue;
        }

        // –еконструкци€ на границе i+1/2
        State W_L_interface, W_R_interface;

        eno3_reconstruction(rho, cell_i, W_L_interface.rho, W_R_interface.rho);
        eno3_reconstruction(u, cell_i, W_L_interface.u, W_R_interface.u);
        eno3_reconstruction(p, cell_i, W_L_interface.p, W_R_interface.p);

        // –ешаем «адачу –имана с реконструированными значени€ми
        const State state_at_interface = solve_riemann_problem(W_L_interface, W_R_interface, 0.0, cfg, cfg.approx_type);
        fluxes[i] = physToFlux(state_at_interface, gamma);
    }
}

