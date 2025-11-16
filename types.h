#pragma once

#include <vector>
#include <string>


// состояние в физических переменных (ρ, u, p)
struct State {
    double rho = 0.0; 
    double u = 0.0;   
    double p = 0.0;   
};

// состояние в консервативных переменных (U1, U2, U3)
struct Conserved {
    double rho = 0.0;  // U1: Плотность массы
    double rhou = 0.0; // U2: Плотность импульса
    double E = 0.0;    // U3: Полная энергия
};

// поток консервативных переменных (F1, F2, F3)
struct Flux {
    double rho_f = 0.0;
    double rhou_f = 0.0;
    double E_f = 0.0;
};


struct PhysicalParameters {
    int test_case = 0;
    double gamma;
    State left;  
    State right; 
};

struct GridParameters {
    double x_min, x_max;
    double x_diaphragm;
    int Nx; 
    double t_final;
    double CFL;
};


struct Grid {
    int Nx; // количество реальных ячеек
    int num_fict; // количество фиктивных ячеек 
    double dx;

    std::vector<double> x_centers; 
    std::vector<State> W;          // вектор физических переменных
    std::vector<Conserved> U;      // вектор консервативных переменных

    
    Grid(int nx, int fict) : Nx(nx), num_fict(fict) {
        const int total_cells = nx + 2 * fict;
        x_centers.resize(total_cells);
        W.resize(total_cells);
        U.resize(total_cells);
    }
};



enum class NumericalMethod { 
    GODUNOV, 
    ACOUSTIC,
    KOLGAN,
    RODIONOV
};


enum class ApproximationType {
    RAREFACTION,
    SHOCK,
    ACOUSTIC,
    AVERAGE
};

enum class TypesOfVarForReconstruction {
    CONSERVATIVE,
    NONCONSERVATIVE
};

struct Config {
    PhysicalParameters phys;
    GridParameters grid;
    NumericalMethod method;
    ApproximationType approx_type;
    TypesOfVarForReconstruction var_type;
};


Config readConfig(const std::string& filename);
// объявление функции для работы с аналитическим решением
void save_analytical_solution(int test_num, const std::string& outputFilename);
