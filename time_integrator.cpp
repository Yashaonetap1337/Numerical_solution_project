#include "time_integrator.h"
#include "euler_utils.h"
#include "boundary_conditions.h"
#include <vector>
#include <stdexcept>


#include "godunov.h"
#include "kolgan.h"
#include "eno.h"
#include "weno.h"
#include "acoustic.h"

// Вычисление правой части L(U)
// Эта функция является "мостом" между пространственной и временной дискретизацией.
// Она принимает текущее состояние сетки `grid` и вычисляет производную по пространству,
// результат (вектор) записывает в `rhs`.
static void compute_rhs(const Grid& grid, const Config& cfg, std::vector<Conserved>& rhs) {
    const double dx = grid.dx;
    std::vector<Flux> fluxes(grid.Nx + 1);

    switch (cfg.method) {
    case NumericalMethod::GODUNOV:
        godunov_flux_computation(grid, cfg, fluxes);
        break;
    case NumericalMethod::KOLGAN:
        kolgan_flux_computation(grid, cfg, fluxes);
        break;
    case NumericalMethod::ENO:
        eno_flux_computation(grid, cfg, fluxes);
        break;
    case NumericalMethod::WENO:
        weno_flux_computation(grid, cfg, fluxes);
        break;
    case NumericalMethod::ACOUSTIC:
        acoustic_flux_computation(grid, cfg, fluxes);
        break;
    default:
        throw std::runtime_error("Unknown spatial method in compute_rhs");
    }

    //Вычисляем правую часть L(U) = - (F_{i+1/2} - F_{i-1/2}) / dx
    // Убедимся, что выходной вектор имеет правильный размер.
    rhs.resize(grid.U.size());

    // Проходим по всем РЕАЛЬНЫМ ячейкам
    for (int i = 0; i < grid.Nx; ++i) {
        // Глобальный индекс текущей реальной ячейки
        const int cell_idx = i + grid.num_fict;

        // Потоки на левой и правой границах этой ячейки
        const Flux& F_left = fluxes[i];
        const Flux& F_right = fluxes[i + 1];

        // Вычисляем и сохраняем пространственную производную
        rhs[cell_idx].rho = -(F_right.rho_f - F_left.rho_f) / dx;
        rhs[cell_idx].rhou = -(F_right.rhou_f - F_left.rhou_f) / dx;
        rhs[cell_idx].E = -(F_right.E_f - F_left.E_f) / dx;
    }
}

//  Интегратор Эйлера
static void euler_step(Grid& grid, double dt, const Config& cfg) {
    std::vector<Conserved> rhs;
    compute_rhs(grid, cfg, rhs); // Вычисляем L(U^n)

    // U^{n+1} = U^n + dt * L(U^n)
    for (size_t i = 0; i < grid.U.size(); ++i) {
        grid.U[i] = grid.U[i] + dt * rhs[i];
    }
}

//  Интегратор TVD RK3 
static void tvd_rk3_step(Grid& grid, double dt, const Config& cfg) {
    const int total_cells = grid.U.size();

    std::vector<Conserved> U_n = grid.U; // Сохраняем U^n

    // Cтадия 1 -> U^{(1)} 
    std::vector<Conserved> rhs1;
    compute_rhs(grid, cfg, rhs1);
    for (int i = 0; i < total_cells; ++i) grid.U[i] = U_n[i] + dt * rhs1[i];
    apply_boundary_conditions(grid, cfg);
    for (int i = 0; i < total_cells; ++i) grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);

    // Стадия 2 -> U^{(2)}
    std::vector<Conserved> rhs2;
    compute_rhs(grid, cfg, rhs2);
    for (int i = 0; i < total_cells; ++i) grid.U[i] = (3.0 / 4.0) * U_n[i] + (1.0 / 4.0) * grid.U[i] + (1.0 / 4.0) * dt * rhs2[i];
    apply_boundary_conditions(grid, cfg);
    for (int i = 0; i < total_cells; ++i) grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);

    // Стадия 3 -> U^{n+1}
    std::vector<Conserved> rhs3;
    compute_rhs(grid, cfg, rhs3);
    for (int i = 0; i < total_cells; ++i) grid.U[i] = (1.0 / 3.0) * U_n[i] + (2.0 / 3.0) * grid.U[i] + (2.0 / 3.0) * dt * rhs3[i];
}


void time_step(Grid& grid, double dt, const Config& cfg) {

    if (cfg.time_integrator == TimeIntegrator::EULER) {
        euler_step(grid, dt, cfg);
    }
    else if (cfg.time_integrator == TimeIntegrator::TVD_RK3) {
        tvd_rk3_step(grid, dt, cfg);
    }
    else {
        throw std::runtime_error("Unknown time integrator selected!");
    }

    // Финальное применение ГУ и обновление W после полного шага
    apply_boundary_conditions(grid, cfg);
    for (size_t i = 0; i < grid.U.size(); ++i) {
        grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);
    }
}