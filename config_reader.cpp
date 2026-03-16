#include "types.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

static double stringToDouble(std::string str) {
    std::replace(str.begin(), str.end(), ',', '.');
    try {
        return std::stod(str);
    }
    catch (...) {
        throw std::invalid_argument("Invalid number: " + str);
    }
}

static NumericalMethod parseMethod(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "godunov") return NumericalMethod::GODUNOV;
    if (lower == "acoustic") return NumericalMethod::ACOUSTIC;
    if (lower == "kolgan") return NumericalMethod::KOLGAN;
    if (lower == "rodionov") return NumericalMethod::RODIONOV;
    if (lower == "eno") return NumericalMethod::ENO;
    if (lower == "weno") return NumericalMethod::WENO;
    if (lower == "maccormack") return NumericalMethod::MACCORMACK;
    if (lower == "mader") return NumericalMethod::MADER;
    throw std::runtime_error("Unknown method: " + s);
}

static ApproximationType parseApproxType(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "rarefaction") return ApproximationType::RAREFACTION;
    if (lower == "shock") return ApproximationType::SHOCK;
    if (lower == "acoustic") return ApproximationType::ACOUSTIC;
    if (lower == "average") return ApproximationType::AVERAGE;
    throw std::runtime_error("Unknown approx type: " + s);
}

static TimeIntegrator parseTimeIntegrator(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "euler") return TimeIntegrator::EULER;
    if (lower == "tvd_rk3") return TimeIntegrator::TVD_RK3;
    throw std::runtime_error("Unknown time integrator: " + s);
}

static RiemannSolverType parseRiemannSolver(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "exact") return RiemannSolverType::EXACT;
    if (lower == "acoustic") return RiemannSolverType::ACOUSTIC;
    if (lower == "roe") return RiemannSolverType::ROE;
    if (lower == "hll") return RiemannSolverType::HLL;
    if (lower == "hllc") return RiemannSolverType::HLLC;
    if (lower == "rusanov") return RiemannSolverType::RUSANOV;
    if (lower == "osher") return RiemannSolverType::OSHER;
    throw std::runtime_error("Unknown Riemann solver: " + s);
}

static TypesOfVarForReconstruction parseVarType(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "conservative") return TypesOfVarForReconstruction::CONSERVATIVE;
    if (lower == "nonconservative") return TypesOfVarForReconstruction::NONCONSERVATIVE;
    throw std::runtime_error("Unknown var type: " + s);
}

static SnapshotOutputType parseSnapshotOutput(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "none") return SnapshotOutputType::NONE;
    if (lower == "by_steps") return SnapshotOutputType::BY_STEPS;
    if (lower == "by_time") return SnapshotOutputType::BY_TIME;
    throw std::runtime_error("Unknown snapshot output: " + s);
}

static FluxCorrectionType parseFluxCorrection(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "none") return FluxCorrectionType::NONE;
    if (lower == "viscosity") return FluxCorrectionType::VISCOSITY;
    if (lower == "fct") return FluxCorrectionType::FCT;
    return FluxCorrectionType::NONE;
}

static LimiterType parseLimiter(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "none") return LimiterType::NONE;
    if (lower == "minmod") return LimiterType::MINMOD;
    if (lower == "superbee") return LimiterType::SUPERBEE;
    if (lower == "van_leer" || lower == "vanleer") return LimiterType::VAN_LEER;
    if (lower == "kolgan_1972") return LimiterType::KOLGAN_1972;
    if (lower == "kolgan_1975") return LimiterType::KOLGAN_1975;
    if (lower == "osher_1984") return LimiterType::OSHER_1984;
    throw std::runtime_error("Unknown limiter: " + s);
}

Config readConfig(const std::string& filename) {
    Config cfg;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open config file");

    std::string line, section;
    int line_num = 0;
    while (std::getline(file, line)) {
        line_num++;
        size_t comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);
        line = trim(line);
        if (line.empty()) continue;

        if (line.front() == '[' && line.back() == ']') {
            section = trim(line.substr(1, line.size() - 2));
            continue;
        }

        size_t eq = line.find('=');
        if (eq == std::string::npos)
            throw std::runtime_error("Invalid line " + std::to_string(line_num));

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        if (key.empty() || val.empty())
            throw std::runtime_error("Empty key/value line " + std::to_string(line_num));

        try {
            if (section == "PhysicalParameters") {
                double v = stringToDouble(val);
                if (key == "test_case") cfg.phys.test_case = std::stoi(val);
                else if (key == "gamma") cfg.phys.gamma = v;
                else if (key == "rho_left_bottom") cfg.phys.left_bottom.rho = v;
                else if (key == "u_left_bottom") cfg.phys.left_bottom.u = v;
                else if (key == "v_left_bottom") cfg.phys.left_bottom.v = v;
                else if (key == "p_left_bottom") cfg.phys.left_bottom.p = v;
                else if (key == "rho_left_top") cfg.phys.left_top.rho = v;
                else if (key == "u_left_top") cfg.phys.left_top.u = v;
                else if (key == "v_left_top") cfg.phys.left_top.v = v;
                else if (key == "p_left_top") cfg.phys.left_top.p = v;
                else if (key == "rho_right_bottom") cfg.phys.right_bottom.rho = v;
                else if (key == "u_right_bottom") cfg.phys.right_bottom.u = v;
                else if (key == "v_right_bottom") cfg.phys.right_bottom.v = v;
                else if (key == "p_right_bottom") cfg.phys.right_bottom.p = v;
                else if (key == "rho_right_top") cfg.phys.right_top.rho = v;
                else if (key == "u_right_top") cfg.phys.right_top.u = v;
                else if (key == "v_right_top") cfg.phys.right_top.v = v;
                else if (key == "p_right_top") cfg.phys.right_top.p = v;
                else if (key == "Rgas") cfg.phys.Rgas = v;
                else if (key == "M_molar") cfg.phys.M_molar = v;
                else if (key == "E_act") cfg.phys.E_act = v;
                else if (key == "Z_freq") cfg.phys.Z_freq = v;
                else if (key == "VISC") cfg.phys.VISC = v;
                else if (key == "MINWT") cfg.phys.MINWT = v;
                else if (key == "MINGRHO") cfg.phys.MINGRHO = v;
                else if (key == "GASW") cfg.phys.GASW = v;
         
                else throw std::runtime_error("Unknown physical key: " + key);
            }
            else if (section == "GridParameters") {
                if (key == "x_min") cfg.grid.x_min = stringToDouble(val);
                else if (key == "x_max") cfg.grid.x_max = stringToDouble(val);
                else if (key == "y_min") cfg.grid.y_min = stringToDouble(val);
                else if (key == "y_max") cfg.grid.y_max = stringToDouble(val);
                else if (key == "x_diaphragm") cfg.grid.x_diaphragm = stringToDouble(val);
                else if (key == "y_diaphragm") cfg.grid.y_diaphragm = stringToDouble(val);
                else if (key == "Nx") cfg.grid.Nx = std::stoi(val);
                else if (key == "Ny") cfg.grid.Ny = std::stoi(val);
                else if (key == "t_final") cfg.grid.t_final = stringToDouble(val);
                else if (key == "CFL") cfg.grid.CFL = stringToDouble(val);
                else throw std::runtime_error("Unknown grid key: " + key);
            }
            else if (section == "NumericalMethod") {
                if (key == "method") cfg.method = parseMethod(val);
                else if (key == "approximation") cfg.approx_type = parseApproxType(val);
                else if (key == "reconstruction") cfg.var_type = parseVarType(val);
                else if (key == "time_integrator") cfg.time_integrator = parseTimeIntegrator(val);
                else if (key == "riemann_solver") cfg.riemann_solver_type = parseRiemannSolver(val);
                else if (key == "flux_correction") cfg.flux_correction = parseFluxCorrection(val);
                else if (key == "visc_coeff") cfg.viscosity_coeff = stringToDouble(val);
                else if (key == "limiter") cfg.limiter = parseLimiter(val);
                else if (key == "slope_limiter") cfg.slope_limiter = parseLimiter(val);
                else if (key == "flux_limiter") cfg.flux_limiter = parseLimiter(val);
                else throw std::runtime_error("Unknown method key: " + key);
            }
            else if (section == "Output") {
                if (key == "snapshot_output") cfg.output.snapshot_output = parseSnapshotOutput(val);
                else if (key == "snapshot_interval_steps") cfg.output.snapshot_interval_steps = std::stoi(val);
                else if (key == "snapshot_interval_time") cfg.output.snapshot_interval_time = stringToDouble(val);
                else if (key == "snapshots_directory") cfg.output.snapshots_directory = val;
                else if (key == "run_comparison") {
                    std::string low = val;
                    std::transform(low.begin(), low.end(), low.begin(), ::tolower);
                    cfg.run_riemann_comparison = (low == "true" || low == "1" || low == "yes");
                }
                else throw std::runtime_error("Unknown output key: " + key);
            }
            else if (section == "BoundaryConditions") {
                if (key == "left_boundary") {
                    if (val == "WALL") cfg.left_boundary = BoundaryType::WALL;
                    else if (val == "FREE") cfg.left_boundary = BoundaryType::FREE;
                    else if (val == "PERIODIC") cfg.left_boundary = BoundaryType::PERIODIC;
                    else throw std::runtime_error("Invalid left boundary");
                }
                else if (key == "right_boundary") {
                    if (val == "WALL") cfg.right_boundary = BoundaryType::WALL;
                    else if (val == "FREE") cfg.right_boundary = BoundaryType::FREE;
                    else if (val == "PERIODIC") cfg.right_boundary = BoundaryType::PERIODIC;
                    else throw std::runtime_error("Invalid right boundary");
                }
                else if (key == "bottom_boundary") {
                    if (val == "WALL") cfg.bottom_boundary = BoundaryType::WALL;
                    else if (val == "FREE") cfg.bottom_boundary = BoundaryType::FREE;
                    else if (val == "PERIODIC") cfg.bottom_boundary = BoundaryType::PERIODIC;
                    else throw std::runtime_error("Invalid bottom boundary");
                }
                else if (key == "top_boundary") {
                    if (val == "WALL") cfg.top_boundary = BoundaryType::WALL;
                    else if (val == "FREE") cfg.top_boundary = BoundaryType::FREE;
                    else if (val == "PERIODIC") cfg.top_boundary = BoundaryType::PERIODIC;
                    else throw std::runtime_error("Invalid top boundary");
                }
                else throw std::runtime_error("Unknown boundary key: " + key);
            }
            else throw std::runtime_error("Unknown section: " + section);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Error line " + std::to_string(line_num) + ": " + e.what());
        }
    }

    if (cfg.grid.x_min >= cfg.grid.x_max) throw std::runtime_error("x_min >= x_max");
    if (cfg.grid.y_min >= cfg.grid.y_max) throw std::runtime_error("y_min >= y_max");
    if (cfg.grid.Nx <= 0 || cfg.grid.Ny <= 0) throw std::runtime_error("Nx/Ny must be positive");
    if (cfg.grid.t_final <= 0) throw std::runtime_error("t_final must be positive");
    if (cfg.grid.CFL <= 0 || cfg.grid.CFL > 1) throw std::runtime_error("CFL must be in (0,1]");
    if (cfg.phys.gamma <= 1) throw std::runtime_error("gamma must be >1");
    if (cfg.phys.left_bottom.rho <= 0 || cfg.phys.left_top.rho <= 0 ||
        cfg.phys.right_bottom.rho <= 0 || cfg.phys.right_top.rho <= 0)
        throw std::runtime_error("Density must be positive");
    if (cfg.phys.left_bottom.p <= 0 || cfg.phys.left_top.p <= 0 ||
        cfg.phys.right_bottom.p <= 0 || cfg.phys.right_top.p <= 0)
        throw std::runtime_error("Pressure must be positive");

    return cfg;
}