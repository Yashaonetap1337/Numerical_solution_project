#include "solver.h" 
#include "grid.h"   
#include "boundary_conditions.h"
#include "godunov.h" 
#include "acoustic.h"
#include "kolgan.h"
#include "rodionov.h"
#include "time_integrator.h"
#include "eno.h"     
#include "weno.h"     
#include "euler_utils.h"
#include "riemann_solver.h"
#include "analytical.h"
#include "maccormack.h"
#include "flux_correction.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <sys/stat.h> 


// функция, которая определяет необходимое количество фиктивных ячеек для каждого метода
static int get_required_fict_cells(const Config& cfg) {
    int method_req = 0;
    switch (cfg.method) {
    case NumericalMethod::GODUNOV:
        method_req = 1;
        break;
    case NumericalMethod::ACOUSTIC:
        method_req = 1;
        break;
    case NumericalMethod::KOLGAN:
        method_req = 2;
        break;
    case NumericalMethod::RODIONOV:
        method_req = 2;
        break;
    case NumericalMethod::ENO:
        method_req = 4;
        break;
    case NumericalMethod::WENO:
        method_req = 3;
        break;
    case NumericalMethod::MACCORMACK: 
        method_req = 1;
        break;
    default:
        throw std::runtime_error("The number of fictitious cells is not defined for the selected method!");
    }
    // Если включена коррекция FCT или Вязкость (вязкость часто требует 1 соседа, FCT - 2)
    if (cfg.flux_correction == FluxCorrectionType::FCT) {
        return std::max(method_req, 2);
    }
    return method_req;
}


// функция для расчета шага по времени. 
static double calculate_timestep(const Grid& grid, const Config& cfg) {
    double max_lambda = 0.0;
    // ищем максимальную скорость волны по всем реальным ячейкам
    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i) {
        double a = soundSpeed(grid.W[i], cfg.phys.gamma);
        max_lambda = std::max(max_lambda, std::abs(grid.W[i].u) + a);
    }
    return cfg.grid.CFL * grid.dx / max_lambda;
}



void save_snapshot(const Grid& grid, const Config& cfg, int step, double t, const std::string& filename) {
    // Создаем директорию, если она не существует
    // std::filesystem::create_directories(cfg.output.snapshots_directory);

    std::string full_filename = cfg.output.snapshots_directory + "/" + filename;
    std::ofstream file(full_filename);

    if (!file.is_open()) {
        std::cerr << "Warning: Could not open snapshot file " << full_filename << std::endl;
        return;
    }

    file << "x,rho,u,p,e,step,time\n";
    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i) {
        State W = consToPhys(grid.U[i], cfg.phys.gamma);
        double e = (W.rho > 1e-9) ? W.p / (W.rho * (cfg.phys.gamma - 1.0)) : 0.0;
        file << grid.x_centers[i] << "," << W.rho << "," << W.u << "," << W.p << "," << e
            << "," << step << "," << t << "\n";
    }
    file.close();

    std::cout << "Numerical snapshot saved: " << full_filename << " (step=" << step << ", t=" << t << ")" << std::endl;
}

void run_simulation(const Config& cfg, const std::string& outputFilename) {

    int num_fict = get_required_fict_cells(cfg);
    Grid grid(cfg.grid.Nx, num_fict);


    initialize_grid(grid, cfg);

    double t = 0.0;
    int step = 0;

    // Переменные для управления снимками
    int next_snapshot_step = 0;
    double next_snapshot_time = 0.0;

    // Сохраняем начальное состояние
    if (cfg.output.snapshot_output != SnapshotOutputType::NONE) {
        std::string init_filename = "initial_state.csv";
        save_snapshot(grid, cfg, 0, 0.0, init_filename);

        // Сохраняем аналитическое решение для начального момента
        std::string analytical_init_filename = "analytical_initial_state.csv";
        generate_analytical_snapshot(cfg, 0.0, cfg.output.snapshots_directory + "/" + analytical_init_filename);

        // Инициализируем следующие моменты для снимков
        if (cfg.output.snapshot_output == SnapshotOutputType::BY_STEPS) {
            next_snapshot_step = cfg.output.snapshot_interval_steps;
        }
        else if (cfg.output.snapshot_output == SnapshotOutputType::BY_TIME) {
            next_snapshot_time = cfg.output.snapshot_interval_time;
        }
    }

    while (t < cfg.grid.t_final) {

        apply_boundary_conditions(grid, cfg);
        const double gamma = cfg.phys.gamma;
        // пересчитываем физические переменные из консервативных 
        for (int i = 0; i < grid.Nx + 2 * num_fict; ++i) {
            grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);
        }

        // расчет шага по времени 
        double dt = calculate_timestep(grid, cfg);
        if (t + dt > cfg.grid.t_final) {
            dt = cfg.grid.t_final - t;
        }

        if (cfg.method == NumericalMethod::RODIONOV || cfg.method == NumericalMethod::MACCORMACK) {

            for (int i = 0; i < grid.U.size(); ++i) grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);
            apply_boundary_conditions(grid, cfg);

            if (cfg.method == NumericalMethod::RODIONOV) {
                rodionov_step(grid, dt, cfg);
            }
            else {
                maccormack_step(grid, dt, cfg);
            }
        }
        else {
            time_step(grid, dt, cfg);
        }

        apply_flux_correction(grid, cfg);

        t += dt;
        step++;

        // Проверяем, нужно ли сохранять снимок
        bool should_save_snapshot = false;
        std::string snapshot_filename;

        if (cfg.output.snapshot_output == SnapshotOutputType::BY_STEPS && step >= next_snapshot_step) {
            should_save_snapshot = true;
            snapshot_filename = "snapshot_step_" + std::to_string(step) + ".csv";
            next_snapshot_step += cfg.output.snapshot_interval_steps;
        }
        else if (cfg.output.snapshot_output == SnapshotOutputType::BY_TIME && t >= next_snapshot_time) {
            should_save_snapshot = true;
            // Форматируем время для имени файла
            std::string time_str = std::to_string(t);
            size_t dot_pos = time_str.find('.');
            if (dot_pos != std::string::npos) {
                time_str = time_str.substr(0, dot_pos + 3); // берем 2 знака после запятой
            }
            std::replace(time_str.begin(), time_str.end(), '.', '_');
            snapshot_filename = "snapshot_time_" + time_str + ".csv";
            next_snapshot_time += cfg.output.snapshot_interval_time;
        }

        if (should_save_snapshot) {
            // Сохраняем численный снимок
            save_snapshot(grid, cfg, step, t, snapshot_filename);

            // Сохраняем аналитический снимок для того же времени
            std::string analytical_filename = "analytical_" + snapshot_filename;
            generate_analytical_snapshot(cfg, t, cfg.output.snapshots_directory + "/" + analytical_filename);
        }

        if (step % 100 == 0) {
            std::cout << "Step: " << step << ", Time: " << t << "/" << cfg.grid.t_final << std::endl;
        }
    }

    // Сохраняем финальное состояние
    if (cfg.output.snapshot_output != SnapshotOutputType::NONE) {
        std::string final_filename = "final_state.csv";
        save_snapshot(grid, cfg, step, t, final_filename);

        // Сохраняем аналитическое решение для финального времени
        std::string analytical_final_filename = "analytical_final_state.csv";
        generate_analytical_snapshot(cfg, t, cfg.output.snapshots_directory + "/" + analytical_final_filename);
    }

    // сохранение основных результатов
    std::ofstream file(outputFilename);
    if (!file.is_open()) throw std::runtime_error("Cannot open output file: " + outputFilename);


    file << "x,rho,u,p,e\n";
    for (int i = num_fict; i < grid.Nx + num_fict; ++i) {
        State final_W = consToPhys(grid.U[i], cfg.phys.gamma);


        double internal_energy = 0.0;
        if (final_W.rho > 1e-9) {
            internal_energy = final_W.p / (final_W.rho * (cfg.phys.gamma - 1.0));
        }


        file << grid.x_centers[i] << "," << final_W.rho << "," << final_W.u << "," << final_W.p << "," << internal_energy << "\n";
    }

    file.close();
    std::cout << "Simulation finished. Result saved to: " << outputFilename << std::endl;
}