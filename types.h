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

    // --- НАЧАЛО НОВЫХ ОПЕРАТОРОВ ---

    // 1. "Научим" складывать два объекта Conserved (U1 + U2)
    Conserved operator+(const Conserved& other) const {
        return { rho + other.rho, rhou + other.rhou, E + other.E };
    }

    // 2. "Научим" умножать объект Conserved на число (U * a)
    Conserved operator*(double scalar) const {
        return { rho * scalar, rhou * scalar, E * scalar };
    }
};

// 3. "Научим" умножать число на объект Conserved (a * U)
//    Это нужно для выражений вида `dt * rhs1[i]`
inline Conserved operator*(double scalar, const Conserved& cons) {
    return cons * scalar; // Используем уже определенный оператор выше
}


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

enum class SnapshotOutputType {
    NONE,
    BY_STEPS,
    BY_TIME
};

enum class BoundaryType {
    WALL,
    FREE,
    PERIODIC
};

struct OutputParams {
    SnapshotOutputType snapshot_output = SnapshotOutputType::NONE;
    int snapshot_interval_steps = 100;
    double snapshot_interval_time = 0.02;
    std::string snapshots_directory = "snapshots";
};


// Реализуемые методы
enum class NumericalMethod {
    GODUNOV,
    ACOUSTIC,
    KOLGAN,
    RODIONOV,
    ENO,  
    WENO   
};

// Тип аппроксимации начального давления
enum class ApproximationType {
    RAREFACTION,
    SHOCK,
    ACOUSTIC,
    AVERAGE
};


enum class RiemannSolverType {
    EXACT,    
    ACOUSTIC  
};

// Тип используемых переменных
enum class TypesOfVarForReconstruction {
    CONSERVATIVE,
    NONCONSERVATIVE
};

enum class TimeIntegrator {
    EULER,
    TVD_RK3
};


// Лол
struct Config {
    PhysicalParameters phys;
    GridParameters grid;
    NumericalMethod method;
    ApproximationType approx_type;
    TimeIntegrator time_integrator;
    RiemannSolverType riemann_solver_type;
    TypesOfVarForReconstruction var_type;
    OutputParams output;
    BoundaryType left_boundary;
    BoundaryType right_boundary;
};



Config readConfig(const std::string& filename);
// объявление функции для работы с аналитическим решением
void save_analytical_solution(int test_num, const std::string& outputFilename);