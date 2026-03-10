#pragma once

#include <vector>
#include <string>

struct State {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double p = 0.0;
};

struct Conserved {
    double rho = 0.0;
    double rhou = 0.0;
    double rhov = 0.0;
    double E = 0.0;

    Conserved operator+(const Conserved& other) const {
        return { rho + other.rho, rhou + other.rhou, rhov + other.rhov, E + other.E };
    }
    Conserved operator-(const Conserved& other) const {
        return { rho - other.rho, rhou - other.rhou, rhov - other.rhov, E - other.E };
    }
    Conserved operator*(double scalar) const {
        return { rho * scalar, rhou * scalar, rhov * scalar, E * scalar };
    }
};

inline Conserved operator*(double scalar, const Conserved& cons) {
    return cons * scalar;
}

struct Flux {
    double rho_f = 0.0;
    double rhou_f = 0.0;
    double rhov_f = 0.0;
    double E_f = 0.0;

    Flux operator+(const Flux& other) const {
        return { rho_f + other.rho_f, rhou_f + other.rhou_f, rhov_f + other.rhov_f, E_f + other.E_f };
    }
    Flux operator-(const Flux& other) const {
        return { rho_f - other.rho_f, rhou_f - other.rhou_f, rhov_f - other.rhov_f, E_f - other.E_f };
    }
    Flux operator*(double scalar) const {
        return { rho_f * scalar, rhou_f * scalar, rhov_f * scalar, E_f * scalar };
    }
};

inline Conserved operator*(double scalar, const Flux& flux) {
    return { scalar * flux.rho_f, scalar * flux.rhou_f, scalar * flux.rhov_f, scalar * flux.E_f };
}

struct PhysicalParameters {
    int test_case = 0;
    double gamma;
    State left_bottom;
    State left_top;
    State right_bottom;
    State right_top;
};

struct GridParameters {
    double x_min, x_max;
    double y_min, y_max;
    double x_diaphragm, y_diaphragm;
    int Nx, Ny;
    double t_final;
    double CFL;
};

struct Grid {
    int Nx, Ny;
    int num_fict;
    double dx, dy;
    std::vector<double> x_centers;
    std::vector<double> y_centers;
    std::vector<std::vector<State>> W;
    std::vector<std::vector<Conserved>> U;

    Grid(int nx, int ny, int fict)
        : Nx(nx), Ny(ny), num_fict(fict),
        x_centers(nx + 2 * fict),
        y_centers(ny + 2 * fict),
        W(nx + 2 * fict, std::vector<State>(ny + 2 * fict)),
        U(nx + 2 * fict, std::vector<Conserved>(ny + 2 * fict)) {}
};

enum class SnapshotOutputType { NONE, BY_STEPS, BY_TIME };
enum class BoundaryType { WALL, FREE, PERIODIC };
enum class NumericalMethod { GODUNOV, ACOUSTIC, KOLGAN, RODIONOV, ENO, WENO, MACCORMACK };
enum class ApproximationType { RAREFACTION, SHOCK, ACOUSTIC, AVERAGE };
enum class RiemannSolverType { EXACT, ACOUSTIC, HLL, HLLC, ROE, RUSANOV, OSHER };
enum class TypesOfVarForReconstruction { CONSERVATIVE, NONCONSERVATIVE };
enum class TimeIntegrator { EULER, TVD_RK3 };
enum class FluxCorrectionType { NONE, VISCOSITY, FCT };
enum class LimiterType { NONE, MINMOD, SUPERBEE, VAN_LEER, KOLGAN_1972, KOLGAN_1975, OSHER_1984 };

struct OutputParams {
    SnapshotOutputType snapshot_output = SnapshotOutputType::NONE;
    int snapshot_interval_steps = 100;
    double snapshot_interval_time = 0.02;
    std::string snapshots_directory = "snapshots";
};

struct Config {
    PhysicalParameters phys;
    GridParameters grid;
    NumericalMethod method;
    ApproximationType approx_type;
    TimeIntegrator time_integrator;
    RiemannSolverType riemann_solver_type;
    TypesOfVarForReconstruction var_type;
    OutputParams output;
    BoundaryType left_boundary = BoundaryType::FREE;
    BoundaryType right_boundary = BoundaryType::FREE;
    BoundaryType bottom_boundary = BoundaryType::FREE;
    BoundaryType top_boundary = BoundaryType::FREE;
    FluxCorrectionType flux_correction = FluxCorrectionType::NONE;
    double viscosity_coeff = 0.1;
    bool run_riemann_comparison = false;
    LimiterType limiter = LimiterType::MINMOD;
    LimiterType slope_limiter = LimiterType::NONE;
    LimiterType flux_limiter = LimiterType::NONE;
};

Config readConfig(const std::string& filename);