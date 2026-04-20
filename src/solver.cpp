#include "solver.h"
#include "grid.h"
#include "boundary_conditions.h"
#include "euler_utils.h"
#include "time_integrator.h"
#include "rodionov.h"
#include "maccormack.h"
#include "analytical.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "debug_logger.h"
#include "MaderSolver.h"
#include <map>
#include "mesh_io.h"
#include "unstructured_solver.h"

static int get_required_fict_cells(const Config& cfg) {
    int req = 0;
    switch (cfg.method) {
    case NumericalMethod::GODUNOV:
    case NumericalMethod::ACOUSTIC:
    case NumericalMethod::MACCORMACK:
        req = 1; break;
    case NumericalMethod::KOLGAN:
    case NumericalMethod::RODIONOV:
        req = 2; break;
    case NumericalMethod::ENO:
        req = 4; break;
    case NumericalMethod::WENO:
        req = 3; break;
    case NumericalMethod::MADER:
        req = 1; break;
    default:
        throw std::runtime_error("Number of fictitious cells not defined");
    }
    if (cfg.flux_correction == FluxCorrectionType::FCT)
        req = std::max(req, 2);
    return req;
}

static double calculate_timestep(const Grid& grid, const Config& cfg) {
    double max_lambda_x = 0.0, max_lambda_y = 0.0;
    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i)
        for (int j = grid.num_fict; j < grid.Ny + grid.num_fict; ++j) {
            double a = soundSpeed(grid.W[i][j], cfg.phys.gamma);
            max_lambda_x = std::max(max_lambda_x, std::abs(grid.W[i][j].u) + a);
            max_lambda_y = std::max(max_lambda_y, std::abs(grid.W[i][j].v) + a);
        }
    double dt_x = grid.dx / (max_lambda_x + 1e-12);
    double dt_y = grid.dy / (max_lambda_y + 1e-12);
    return cfg.grid.CFL * std::min(dt_x, dt_y);
}

void save_snapshot(const Grid& grid, const Config& cfg, int step, double t, const std::string& filename) {
    std::string full = cfg.output.snapshots_directory + "/" + filename;
    std::ofstream f(full);
    if (!f.is_open()) {
        std::cerr << "Cannot open " << full << std::endl;
        return;
    }
    f << "x,y,rho,u,v,p,E,step,time\n";
    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i)
        for (int j = grid.num_fict; j < grid.Ny + grid.num_fict; ++j)
            f << grid.x_centers[i] << "," << grid.y_centers[j] << ","
              << grid.W[i][j].rho << "," << grid.W[i][j].u << ","
              << grid.W[i][j].v << "," << grid.W[i][j].p << ","
              << grid.U[i][j].E << "," << step << "," << t << "\n";
    f.close();
}

void run_simulation(const Config& cfg, const std::string& outputFilename) {

    if (!cfg.mesh_file.empty()) {
        std::map<int, BoundaryType> tag_to_bc = {
            {10, BoundaryType::INFLOW},
            {11, BoundaryType::OUTFLOW},
            {12, BoundaryType::WALL}
        };
        Mesh mesh = load_msh(cfg.mesh_file, tag_to_bc);
        run_unstructured(mesh, cfg, outputFilename);
        return;
    }

    int num_fict = get_required_fict_cells(cfg);
    Grid grid(cfg.grid.Nx, cfg.grid.Ny, num_fict);
    initialize_grid(grid, cfg);

    // --- Инициализация поля массовой доли (только для метода Mader) ---
    OmegaField omega;
    if (cfg.method == NumericalMethod::MADER) {
        omega = make_omega_field(grid);
        // По умолчанию всё — непрореагировавшее ВВ (omega = 1)
        for (int i = 0; i < grid.Nx + 2 * num_fict; ++i)
            for (int j = 0; j < grid.Ny + 2 * num_fict; ++j)
                omega[i][j] = 1.0;

        if (cfg.phys.test_case == 6) {
            // Специальная инициализация для Mader3 (обтекание угла)
            const double x_corner = 2.0, y_corner = 3.0, x_det = 0.15;
            for (int i = 0; i < grid.Nx + 2 * num_fict; ++i) {
                double x = grid.x_centers[i];
                for (int j = 0; j < grid.Ny + 2 * num_fict; ++j) {
                    double y = grid.y_centers[j];
                    if ((x < x_corner && y < y_corner) ||   // стенка
                        (x < x_det && y >= y_corner)) {     // детонатор (уже сгорел)
                        omega[i][j] = 0.0;
                    }
                }
            }
        } else {
            // Стандартная инициализация для других тестов:
            // все ячейки левее x_diaphragm считаются детонатором (omega = 0)
            for (int i = 0; i < grid.Nx + 2 * num_fict; ++i)
                for (int j = 0; j < grid.Ny + 2 * num_fict; ++j)
                    if (grid.x_centers[i] < cfg.grid.x_diaphragm)
                        omega[i][j] = 0.0;
        }
    }

    log_initialization(grid, cfg);

    std::ofstream profile_log;
    open_profile_log(profile_log);

    double t = 0.0;
    int step = 0;

    int next_snapshot_step = 0;
    double next_snapshot_time = 0.0;
    if (cfg.output.snapshot_output != SnapshotOutputType::NONE) {
        save_snapshot(grid, cfg, 0, 0.0, "initial_state.csv");
        if (cfg.output.snapshot_output == SnapshotOutputType::BY_STEPS)
            next_snapshot_step = cfg.output.snapshot_interval_steps;
        else if (cfg.output.snapshot_output == SnapshotOutputType::BY_TIME)
            next_snapshot_time = cfg.output.snapshot_interval_time;
    }

    while (t < cfg.grid.t_final) {
        for (int i = 0; i < grid.Nx + 2 * num_fict; ++i)
            for (int j = 0; j < grid.Ny + 2 * num_fict; ++j)
                grid.W[i][j] = consToPhys(grid.U[i][j], cfg.phys.gamma);
        apply_boundary_conditions(grid, cfg);

        double dt = calculate_timestep(grid, cfg);
        if (t + dt > cfg.grid.t_final) dt = cfg.grid.t_final - t;

        // Отслеживание фронта для MADER (опционально)
        if (cfg.method == NumericalMethod::MADER && step % 10 == 0) {
            double x_front = grid.x_centers[num_fict];
            int i_front = num_fict;
            for (int i = grid.Nx + num_fict - 1; i >= num_fict; --i) {
                int j = num_fict;
                if (omega[i][j] < 0.99) {
                    x_front = grid.x_centers[i];
                    i_front = i;
                    break;
                }
            }
            double P_max = 0.0;
            for (int i = std::max(num_fict, i_front - 3);
                 i <= std::min(grid.Nx + num_fict - 1, i_front + 5); ++i) {
                int j = num_fict;
                if (grid.W[i][j].p > P_max) P_max = grid.W[i][j].p;
            }
            double u_plateau = 0.0; int cnt = 0;
            for (int i = num_fict; i <= std::max(num_fict, i_front - 10); ++i) {
                int j = num_fict;
                u_plateau += grid.W[i][j].u; cnt++;
            }
            if (cnt > 0) u_plateau /= cnt;

            static std::ofstream front_log;
            if (!front_log.is_open()) {
                front_log.open("results/front_position.csv");
                front_log << "t,x_front,P_max,u_plateau\n";
            }
            front_log << t << "," << x_front << "," << P_max << "," << u_plateau << "\n";
        }

        // Вызов метода
        if (cfg.method == NumericalMethod::MADER) {
            std::vector<std::vector<State>> W_new;
            MaderTimeStep(grid, cfg, grid.W, W_new, omega, dt, t);
            for (int i = 0; i < grid.Nx + 2 * num_fict; ++i)
                for (int j = 0; j < grid.Ny + 2 * num_fict; ++j) {
                    grid.W[i][j] = W_new[i][j];
                    grid.U[i][j] = physToCons(W_new[i][j], cfg.phys.gamma);
                }
        }
        else if (cfg.method == NumericalMethod::RODIONOV) {
            rodionov_step(grid, dt, cfg);
        }
        else if (cfg.method == NumericalMethod::MACCORMACK) {
            maccormack_step(grid, dt, cfg);
        }
        else {
            time_step(grid, dt, cfg);
        }

        t += dt;
        step++;

        bool save = false;
        std::string snap_name;
        if (cfg.output.snapshot_output == SnapshotOutputType::BY_STEPS && step >= next_snapshot_step) {
            save = true;
            snap_name = "snapshot_step_" + std::to_string(step) + ".csv";
            next_snapshot_step += cfg.output.snapshot_interval_steps;
        }
        else if (cfg.output.snapshot_output == SnapshotOutputType::BY_TIME && t >= next_snapshot_time) {
            save = true;
            std::string ts = std::to_string(t);
            std::replace(ts.begin(), ts.end(), '.', '_');
            snap_name = "snapshot_time_" + ts + ".csv";
            next_snapshot_time += cfg.output.snapshot_interval_time;
        }
        if (save) save_snapshot(grid, cfg, step, t, snap_name);

        if (step % 100 == 0)
            std::cout << "Step: " << step << ", t: " << t << "/" << cfg.grid.t_final << std::endl;
        if (step % 50 == 0 || step == 1)
            log_y_profile(grid, cfg, step, t, profile_log);
    }

    profile_log.close();

    std::ofstream fout(outputFilename);
    if (!fout.is_open()) throw std::runtime_error("Cannot open output file");
    bool save_omega = (cfg.method == NumericalMethod::MADER);
    fout << "x,y,rho,u,v,p,E" << (save_omega ? ",omega" : "") << "\n";
    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i)
        for (int j = grid.num_fict; j < grid.Ny + grid.num_fict; ++j) {
            fout << grid.x_centers[i] << "," << grid.y_centers[j] << ","
                 << grid.W[i][j].rho << "," << grid.W[i][j].u << ","
                 << grid.W[i][j].v << "," << grid.W[i][j].p << ","
                 << grid.U[i][j].E;
            if (save_omega) fout << "," << omega[i][j];
            fout << "\n";
        }
    fout.close();

    if (cfg.output.snapshot_output != SnapshotOutputType::NONE)
        save_snapshot(grid, cfg, step, t, "final_state.csv");

    std::cout << "Simulation finished. Results in " << outputFilename << std::endl;
}