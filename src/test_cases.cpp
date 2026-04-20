#include "test_cases.h"
#include <stdexcept>

void setup_test_case(Config& cfg) {
    if (cfg.phys.test_case == 0) return;

    switch (cfg.phys.test_case) {
    case 1: // Sod shock tube
        cfg.phys.left_bottom = { 1.0, 0.75, 0.0, 1.0 };
        cfg.phys.left_top = { 1.0, 0.75, 0.0, 1.0 };
        cfg.phys.right_bottom = { 0.125, 0.0, 0.0, 0.1 };
        cfg.phys.right_top = { 0.125, 0.0, 0.0, 0.1 };
        cfg.grid.x_diaphragm = 0.5;
        cfg.grid.y_diaphragm = 0.5;
        cfg.grid.t_final = 0.2;
        break;

    case 2: // 123 problem (two rarefactions)
        cfg.phys.left_bottom = { 1.0, -2.0, 0.0, 0.4 };
        cfg.phys.left_top = { 1.0, -2.0, 0.0, 0.4 };
        cfg.phys.right_bottom = { 1.0, 2.0, 0.0, 0.4 };
        cfg.phys.right_top = { 1.0, 2.0, 0.0, 0.4 };
        cfg.grid.x_diaphragm = 0.5;
        cfg.grid.y_diaphragm = 0.5;
        cfg.grid.t_final = 0.15;
        break;

    case 3: // Very strong shock
        cfg.phys.left_bottom = { 1.0, 0.0, 0.0, 1000.0 };
        cfg.phys.left_top = { 1.0, 0.0, 0.0, 1000.0 };
        cfg.phys.right_bottom = { 1.0, 0.0, 0.0, 0.01 };
        cfg.phys.right_top = { 1.0, 0.0, 0.0, 0.01 };
        cfg.grid.x_diaphragm = 0.5;
        cfg.grid.y_diaphragm = 0.5;
        cfg.grid.t_final = 0.012;
        break;

    case 4: // Two interacting blast waves
        cfg.phys.left_bottom = { 5.99924, 19.5975, 0.0, 460.894 };
        cfg.phys.left_top = { 5.99924, 19.5975, 0.0, 460.894 };
        cfg.phys.right_bottom = { 5.99242, -6.19633, 0.0, 46.0950 };
        cfg.phys.right_top = { 5.99242, -6.19633, 0.0, 46.0950 };
        cfg.grid.x_diaphragm = 0.4;
        cfg.grid.y_diaphragm = 0.5;
        cfg.grid.t_final = 0.035;
        break;

    case 5: // Strong rarefaction
        cfg.phys.left_bottom = { 1.0, -19.59745, 0.0, 1000.0 };
        cfg.phys.left_top = { 1.0, -19.59745, 0.0, 1000.0 };
        cfg.phys.right_bottom = { 1.0, -19.59745, 0.0, 0.01 };
        cfg.phys.right_top = { 1.0, -19.59745, 0.0, 0.01 };
        cfg.grid.x_diaphragm = 0.8;
        cfg.grid.y_diaphragm = 0.5;
        cfg.grid.t_final = 0.012;
        break;

    case 6: // Mader3 – обтекание угла (параметры берутся из конфига)
        // Ничего не делаем – оставляем значения из config.txt
        // Геометрия задаётся в initialize_grid
        // t_final уже прочитан из конфига
        break;

    default:
        throw std::runtime_error("Unknown test_case number. Available: 1-6.");
    }
}