#include "solver.h" 
#include "grid.h"   
#include "boundary_conditions.h"
#include "godunov.h" 
#include "euler_utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>


// функция, которая определяет необходимое количество фиктивных ячеек для каждого метода
static int get_required_fict_cells(NumericalMethod method) {
    switch (method) {
    case NumericalMethod::GODUNOV:
        return 1;
    case NumericalMethod::ACOUSTIC:
        return 1;
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

    if (cfg.method == NumericalMethod::ACOUSTIC) {
        acoustic_solver(grid, cfg);

        // Сохранение результата напрямую из grid.W
        std::ofstream file(outputFilename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outputFilename);
        }

        file << "x,rho,u,p,e\n";
        for (int i = num_fict; i < grid.Nx + num_fict; ++i) {
            const State& W = grid.W[i];
            double internal_energy = (W.rho > 1e-9)
                ? W.p / (W.rho * (cfg.phys.gamma - 1.0))
                : 0.0;
            file << grid.x_centers[i] << ","
                << W.rho << ","
                << W.u << ","
                << W.p << ","
                << internal_energy << "\n";
        }
        file.close();
        std::cout << "Acoustic solution finished. Result saved to: " << outputFilename << std::endl;
        return;
    }

    double t = 0.0;
    int step = 0;
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

        // численный метод
        if (cfg.method == NumericalMethod::GODUNOV) {
            godunov_step(grid, dt, cfg); 
        }
        
        else {
            //throw std::runtime_error("Unknown or not implemented numerical method selected!");
            break;
        }

        
        t += dt;
        step++;

        
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