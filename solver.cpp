#include "solver.h" 
#include "grid.h"   
#include "boundary_conditions.h"
#include "godunov.h" 
#include "acoustic.h"
#include "kolgan.h"
#include "euler_utils.h"
#include "riemann_solver.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>





// функция для генерации аналитического снэпшота для конкретного времени t
void generate_analytical_snapshot(const Config& cfg, double t, const std::string& outputFilename) {
    std::ofstream file(outputFilename);
    if (!file.is_open()) return; // просто пропускаем, если не удалось
    file << "x,rho,u,p,e\n";

    const int points = 2000;
    const double dx = (cfg.grid.x_max - cfg.grid.x_min) / points;

    for (int i = 0; i <= points; ++i) {
        const double x = cfg.grid.x_min + i * dx;
        const double xi = (t > 1e-9) ? (x - cfg.grid.x_diaphragm) / t : 0.0;
        State W = solve_general_riemann_problem(cfg.phys.left, cfg.phys.right, xi, cfg, cfg.approx_type);
        double e = (W.rho > 1e-9) ? W.p / (W.rho * (cfg.phys.gamma - 1.0)) : 0.0;
        file << x << "," << W.rho << "," << W.u << "," << W.p << "," << e << "\n";
    }
    file.close();
}














//void save_snapshot(const Grid& grid, const std::vector<Flux>& fluxes, int step, const Config& cfg) {
//    std::string filename = "snapshot_step" + std::to_string(step) + ".csv";
//    std::ofstream file(filename);
//
//    if (!file.is_open()) {
//        std::cerr << "Warning: Could not open snapshot file " << filename << std::endl;
//        return;
//    }
//
//    file << "cell_idx,x," // Геометрия
//        << "rho,u,p,"       // Физические переменные (W)
//        << "cons_rho,cons_rhou,cons_E," // Консервативные переменные (U)
//        << "flux_rho,flux_rhou,flux_E\n"; // Потоки НА ЛЕВОЙ ГРАНИЦЕ ячейки
//
//    for (int i = 0; i < grid.Nx; ++i) {
//        const int cell_idx = i + grid.num_fict; // Глобальный индекс ячейки
//
//        const State& W = grid.W[cell_idx];
//        const Conserved& U = grid.U[cell_idx];
//        const Flux& F_left = fluxes[i]; // Поток на левой границе ячейки i
//
//        file << i << "," << grid.x_centers[cell_idx] << ","
//            << W.rho << "," << W.u << "," << W.p << ","
//            << U.rho << "," << U.rhou << "," << U.E << ","
//            << F_left.rho_f << "," << F_left.rhou_f << "," << F_left.E_f << "\n";
//    }
//
//    
//    const Flux& F_last = fluxes[grid.Nx];
//    file << grid.Nx << ",,," // Нет ячейки, только граница
//        << ",,,,"
//        << F_last.rho_f << "," << F_last.rhou_f << "," << F_last.E_f << "\n";
//
//    file.close();
//}




// функция, которая определяет необходимое количество фиктивных ячеек для каждого метода
static int get_required_fict_cells(NumericalMethod method) {
    switch (method) {
    case NumericalMethod::GODUNOV:
        return 1; 
    case NumericalMethod::ACOUSTIC:
        return 1;
    case NumericalMethod::KOLGAN:
        return 2;
    default:
        throw std::runtime_error("The number of fictitious cells is not defined for the selected method!");
    }
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



void run_simulation(const Config& cfg, const std::string& outputFilename) {
    
    int num_fict = get_required_fict_cells(cfg.method);
    Grid grid(cfg.grid.Nx, num_fict);

    
    initialize_grid(grid, cfg);

    double t = 0.0;
    int step = 0;

    const int snapshot_interval = 200; // сохранять каждые 200 шагов
    int snapshot_counter = 0;

    //// Сохраняем начальное состояние (t=0)
    //std::string num_fname_t0 = "results/" + outputFilename_base + "_t0.csv";
    //// ... (здесь можно добавить код для сохранения начального численного состояния)
    //std::string an_fname_t0 = "results/analytical_" + outputFilename_base + "_t0.csv";
    //generate_analytical_snapshot(cfg, t, an_fname_t0);
    
    while (t < cfg.grid.t_final) {
        
        apply_boundary_conditions(grid);

        // пересчитываем физические переменные из консервативных 
        for (int i = 0; i < grid.Nx + 2 * num_fict; ++i) {
            grid.W[i] = consToPhys(grid.U[i], cfg.phys.gamma);
        }

        // расчет шага по времени 
        double dt = calculate_timestep(grid, cfg);
        if (t + dt > cfg.grid.t_final) {
            dt = cfg.grid.t_final - t;
        }

        std::vector<Flux> fluxes(grid.Nx + 1);

        if (cfg.method == NumericalMethod::GODUNOV) {
            godunov_step(grid, dt, cfg, fluxes); 
        }
        else if (cfg.method == NumericalMethod::ACOUSTIC) {
            acoustic_step(grid, dt, cfg);
        }
        else if (cfg.method == NumericalMethod::KOLGAN) { 
            kolgan_step(grid, dt, cfg, fluxes);
        }
        else {
            throw std::runtime_error("Unknown or not implemented numerical method selected!");
        }

        
        t += dt;
        step++;


        //if (step % snapshot_interval == 0) {
        //    snapshot_counter++;
        //    std::cout << "Step: " << step << ", Time: " << t << "/" << cfg.grid.t_final
        //        << " (Saving snapshot #" << snapshot_counter << ")" << std::endl;

        //    // Сохраняем численный снэпшот
        //    std::string num_fname = "results/" + outputFilename_base + "_snap" + std::to_string(snapshot_counter) + ".csv";
        //    std::ofstream num_file(num_fname);
        //    num_file << "x,rho,u,p,e\n";
        //    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i) {
        //        State W = consToPhys(grid.U[i], cfg.phys.gamma);
        //        double e = (W.rho > 1e-9) ? W.p / (W.rho * (cfg.phys.gamma - 1.0)) : 0.0;
        //        num_file << grid.x_centers[i] << "," << W.rho << "," << W.u << "," << W.p << "," << e << "\n";
        //    }
        //    num_file.close();

        //    // ГЕНЕРИРУЕМ АНАЛИТИЧЕСКИЙ СНЭПШОТ ДЛЯ ЭТОГО ЖЕ МОМЕНТА ВРЕМЕНИ
        //    std::string an_fname = "results/analytical_" + outputFilename_base + "_snap" + std::to_string(snapshot_counter) + ".csv";
        //    generate_analytical_snapshot(cfg, t, an_fname);
        //}




        //if (step % snapshot_interval == 0) {
        //    std::cout << "Step: " << step << ", Time: " << t << "/" << cfg.grid.t_final
        //        << " (Saving snapshot...)" << std::endl;
        //    save_snapshot(grid, fluxes, step, cfg);
        //}
        //else if (step % 1 == 0) { // Для остальных шагов - старый вывод
        //    std::cout << "Step: " << step << ", Time: " << t << "/" << cfg.grid.t_final << std::endl;
        //}




        
        //// отладочная информация 
        //if (step % 100 == 0) {

        //    double min_pressure = 1e30;
        //    double min_density = 1e30;

        //    double x_at_min_p = 0.0;
        //    double x_at_min_rho = 0.0;

        //    // ищем минимальные p и rho
        //    for (int i = grid.num_fict; i < grid.Nx + grid.num_fict; ++i) {
        //        State current_W = consToPhys(grid.U[i], cfg.phys.gamma);

        //        
        //        if (current_W.p < min_pressure) {
        //            min_pressure = current_W.p;
        //            x_at_min_p = grid.x_centers[i]; 
        //        }

        //        
        //        if (current_W.rho < min_density) {
        //            min_density = current_W.rho;
        //            x_at_min_rho = grid.x_centers[i];
        //        }
        //    }

        //    
        //    std::cout << "Step: " << step
        //        << ", Time: " << t << "/" << cfg.grid.t_final
        //        << ", Min P: " << min_pressure << " at x=" << x_at_min_p
        //        << ", Min Rho: " << min_density << " at x=" << x_at_min_rho
        //        << std::endl;
        //}
        



        if (step % 100 == 0) {
            std::cout << "Step: " << step << ", Time: " << t << "/" << cfg.grid.t_final << std::endl;
        }
    }

    // сохранение результатов
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
