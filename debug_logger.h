#pragma once
#include "types.h"
#include "grid.h"
#include <fstream>
#include <iostream>
#include <cmath>

// ============================================================
// Лог 1: Проверка инициализации сетки
// Выводит значения вдоль центральной вертикальной линии (x~0.5)
// и вдоль центральной горизонтальной линии (y~0.5)
// ============================================================
inline void log_initialization(const Grid& grid, const Config& cfg) {
    std::ofstream f("log_init.csv");
    f << "axis,index,coord,rho,u,v,p\n";

    const int ng = grid.num_fict;

    // Вертикальный срез: фиксируем i ~ середина по x, меняем j
    int i_mid = ng + cfg.grid.Nx / 2;
    for (int j = ng; j < grid.Ny + ng; ++j) {
        const State& W = grid.W[i_mid][j];
        f << "vertical," << j << "," << grid.y_centers[j] << ","
            << W.rho << "," << W.u << "," << W.v << "," << W.p << "\n";
    }

    // Горизонтальный срез: фиксируем j ~ середина по y, меняем i
    int j_mid = ng + cfg.grid.Ny / 2;
    for (int i = ng; i < grid.Nx + ng; ++i) {
        const State& W = grid.W[i][j_mid];
        f << "horizontal," << i << "," << grid.x_centers[i] << ","
            << W.rho << "," << W.u << "," << W.v << "," << W.p << "\n";
    }

    f.close();
    std::cout << "[DEBUG] Initialization log -> log_init.csv\n";
}

// ============================================================
// Лог 2: Проверка флюксов по Y на разрыве
// Выводит y-флюксы вблизи y_diaphragm вдоль x=const
// Вызывать ПОСЛЕ первого вычисления флюксов
// ============================================================
inline void log_y_fluxes(
    const std::vector<std::vector<Flux>>& fluxes_y,
    const Grid& grid,
    const Config& cfg)
{
    std::ofstream f("log_fluxes_y.csv");
    f << "i,j_face,x,rho_f,rhou_f,rhov_f,E_f\n";

    const int ng = grid.num_fict;
    const int Ny = grid.Ny;

    // Найдём индекс грани ближайшей к y_diaphragm
    int j_dia = 0;
    double min_dist = 1e18;
    for (int j = 0; j <= Ny; ++j) {
        // грань j находится между ячейками j+ng-1 и j+ng
        double y_face = cfg.grid.y_min + j * grid.dy;
        double dist = std::abs(y_face - cfg.grid.y_diaphragm);
        if (dist < min_dist) { min_dist = dist; j_dia = j; }
    }

    // Выводим флюксы на 3 гранях вокруг разрыва для всех x
    for (int dj = -2; dj <= 2; ++dj) {
        int jf = j_dia + dj;
        if (jf < 0 || jf > Ny) continue;
        for (int i = ng; i < grid.Nx + ng; ++i) {
            const Flux& F = fluxes_y[i][jf];
            f << i << "," << jf << "," << grid.x_centers[i] << ","
                << F.rho_f << "," << F.rhou_f << "," << F.rhov_f << "," << F.E_f << "\n";
        }
    }

    f.close();
    std::cout << "[DEBUG] Y-flux log (near y_diaphragm, j_face=" << j_dia
        << ") -> log_fluxes_y.csv\n";
}

// ============================================================
// Лог 3: Динамика — срез вдоль Y в разные моменты времени
// Вызывать каждые N шагов из главного цикла
// ============================================================
inline void log_y_profile(const Grid& grid, const Config& cfg,
    int step, double t,
    std::ofstream& f)
{
    const int ng = grid.num_fict;
    int i_mid = ng + cfg.grid.Nx / 2;  // x ~ середина

    for (int j = ng; j < grid.Ny + ng; ++j) {
        const State& W = grid.W[i_mid][j];
        f << step << "," << t << "," << grid.y_centers[j] << ","
            << W.rho << "," << W.u << "," << W.v << "," << W.p << "\n";
    }
}

inline void open_profile_log(std::ofstream& f) {
    f.open("log_y_profile.csv");
    f << "step,t,y,rho,u,v,p\n";
    std::cout << "[DEBUG] Y-profile log -> log_y_profile.csv\n";
}