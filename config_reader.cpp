#include "types.h" 

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cctype>
#include <algorithm>


// функция для удаления пробелов в начале и конце строки
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

// функция для преобразования строки в double,
// корректно обрабатывает как точку, так и запятую в качестве разделителя.
static double stringToDouble(std::string str) {
    std::replace(str.begin(), str.end(), ',', '.');
    try {
        return std::stod(str);
    }
    catch (const std::invalid_argument& e) {
        throw std::invalid_argument("Invalid number format: " + str);
    }
    catch (const std::out_of_range& e) {
        throw std::out_of_range("Number out of range: " + str);
    }
}

// функция для парсинга численного метода
static NumericalMethod parseMethod(const std::string& s) {
    std::string lower = s;
    // преобразуем строку в нижний регистр для регистронезависимого сравнения
    std::transform(lower.begin(), lower.end(), lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    if (lower == "godunov") return NumericalMethod::GODUNOV;
    if (lower == "acoustic") return NumericalMethod::ACOUSTIC;
    if (lower == "kolgan") return NumericalMethod::KOLGAN;
    if (lower == "rodionov") return NumericalMethod::RODIONOV;
    if (lower == "eno") return NumericalMethod::ENO;
    if (lower == "weno") return NumericalMethod::WENO;
    if (lower == "maccormack") return NumericalMethod::MACCORMACK;

    throw std::runtime_error("Unknown numerical method: " + s);
}

static ApproximationType parseApproxType(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    if (lower == "rarefaction") return ApproximationType::RAREFACTION;
    if (lower == "shock")       return ApproximationType::SHOCK;
    if (lower == "acoustic")    return ApproximationType::ACOUSTIC;
    if (lower == "average")     return ApproximationType::AVERAGE;
    throw std::runtime_error("Unknown approximation type: " + s);
}

// функция для парсинга типа интегратора
static TimeIntegrator parseTimeIntegrator(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "euler") return TimeIntegrator::EULER;
    if (lower == "tvd_rk3") return TimeIntegrator::TVD_RK3;
    throw std::runtime_error("Unknown time integrator type: " + s);
}

// функция для парсинга типа решателя Римана
static RiemannSolverType parseRiemannSolverType(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "exact") return RiemannSolverType::EXACT;
    if (lower == "acoustic") return RiemannSolverType::ACOUSTIC;
    if (lower == "roe") return RiemannSolverType::ROE;
    throw std::runtime_error("Unknown Riemann solver type: " + s);
}

static TypesOfVarForReconstruction parseVarType(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    if (lower == "nonconservative") return TypesOfVarForReconstruction::NONCONSERVATIVE;
    if (lower == "conservative")       return TypesOfVarForReconstruction::CONSERVATIVE;
    throw std::runtime_error("Unknown type for reconstruction: " + s);
}

static SnapshotOutputType parseSnapshotOutputType(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    if (lower == "none") return SnapshotOutputType::NONE;
    if (lower == "by_steps") return SnapshotOutputType::BY_STEPS;
    if (lower == "by_time") return SnapshotOutputType::BY_TIME;
    throw std::runtime_error("Unknown snapshot output type: " + s);
}

static FluxCorrectionType parseFluxCorrection(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower == "none") return FluxCorrectionType::NONE;
    if (lower == "viscosity") return FluxCorrectionType::VISCOSITY;
    if (lower == "fct") return FluxCorrectionType::FCT;
    return FluxCorrectionType::NONE;
}

// функция для чтения конфигурационного файла
Config readConfig(const std::string& filename) {
    Config cfg;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + filename);
    }

    std::string line;
    std::string currentSection;

    int line_num = 0;
    while (std::getline(file, line)) {
        line_num++;

        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        line = trim(line);
        if (line.empty()) {
            continue;
        }

        if (line.front() == '[' && line.back() == ']') {
            currentSection = trim(line.substr(1, line.size() - 2));
            continue;
        }

        size_t eqPos = line.find('=');
        if (eqPos == std::string::npos) {
            throw std::runtime_error("Invalid format on line " + std::to_string(line_num) + ": " + line);
        }

        std::string key = trim(line.substr(0, eqPos));
        std::string valueStr = trim(line.substr(eqPos + 1));
        if (key.empty() || valueStr.empty()) {
            throw std::runtime_error("Empty key or value on line " + std::to_string(line_num) + ": " + line);
        }

        try {
            if (currentSection == "PhysicalParameters") {
                double value = stringToDouble(valueStr);
                if (key == "test_case") cfg.phys.test_case = std::stoi(valueStr);
                else if (key == "gamma")     cfg.phys.gamma = value;
                else if (key == "rho_left")  cfg.phys.left.rho = value;
                else if (key == "u_left")    cfg.phys.left.u = value;
                else if (key == "p_left")    cfg.phys.left.p = value;
                else if (key == "rho_right") cfg.phys.right.rho = value;
                else if (key == "u_right")   cfg.phys.right.u = value;
                else if (key == "p_right")   cfg.phys.right.p = value;
                else throw std::runtime_error("Unknown physical parameter: " + key);
            }
            else if (currentSection == "GridParameters") {
                if (key == "x_min")   cfg.grid.x_min = stringToDouble(valueStr);
                else if (key == "x_max")   cfg.grid.x_max = stringToDouble(valueStr);
                else if (key == "x_diaphragm") cfg.grid.x_diaphragm = stringToDouble(valueStr);
                else if (key == "Nx")      cfg.grid.Nx = std::stoi(valueStr);
                else if (key == "t_final") cfg.grid.t_final = stringToDouble(valueStr);
                else if (key == "CFL")     cfg.grid.CFL = stringToDouble(valueStr);
                else throw std::runtime_error("Unknown grid parameter: " + key);
            }
            else if (currentSection == "NumericalMethod") {
                if (key == "method") {
                    cfg.method = parseMethod(valueStr);
                }
                else if (key == "approximation") {
                    cfg.approx_type = parseApproxType(valueStr);
                }
                else if (key == "reconstruction") {
                    cfg.var_type = parseVarType(valueStr);
                }
                else if (key == "time_integrator") {
                    cfg.time_integrator = parseTimeIntegrator(valueStr);
                }
                else if (key == "riemann_solver") {
                    cfg.riemann_solver_type = parseRiemannSolverType(valueStr);
                }
                else if (key == "flux_correction") {
                    cfg.flux_correction = parseFluxCorrection(valueStr);
                }
                else if (key == "visc_coeff") {
                    cfg.viscosity_coeff = stringToDouble(valueStr);
                }
                else {
                    throw std::runtime_error("Unknown method parameter: " + key);
                }
            }
            // секция чтения параметров для хранения данных для скринов
            else if (currentSection == "Output") {
                if (key == "snapshot_output") {
                    cfg.output.snapshot_output = parseSnapshotOutputType(valueStr);
                }
                else if (key == "snapshot_interval_steps") {
                    cfg.output.snapshot_interval_steps = std::stoi(valueStr);
                }
                else if (key == "snapshot_interval_time") {
                    cfg.output.snapshot_interval_time = stringToDouble(valueStr);
                }
                else if (key == "snapshots_directory") {
                    cfg.output.snapshots_directory = valueStr;
                }
                else throw std::runtime_error("Unknown output parameter: " + key);
            }
            else if (currentSection == "BoundaryConditions") {
                if (key == "left_boundary") {
                    if (valueStr == "WALL") cfg.left_boundary = BoundaryType::WALL;
                    else if (valueStr == "FREE") cfg.left_boundary = BoundaryType::FREE;
                    else if (valueStr == "PERIODIC") cfg.left_boundary = BoundaryType::PERIODIC;
                    else throw std::runtime_error("Unknown left boundary type: " + valueStr);
                }
                else if (key == "right_boundary") {
                    if (valueStr == "WALL") cfg.right_boundary = BoundaryType::WALL;
                    else if (valueStr == "FREE") cfg.right_boundary = BoundaryType::FREE;
                    else if (valueStr == "PERIODIC") cfg.right_boundary = BoundaryType::PERIODIC;
                    else throw std::runtime_error("Unknown right boundary type: " + valueStr);
                }
            }
            else if (!currentSection.empty()) {
                throw std::runtime_error("Unknown section: " + currentSection);
            }
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Error parsing line " + std::to_string(line_num) + " (" + line + ") - " + e.what());
        }
    }

    // валидация прочитанных параметров
    if (cfg.grid.x_min >= cfg.grid.x_max)
        throw std::runtime_error("Validation error: x_min must be less than x_max.");
    if (cfg.grid.Nx <= 0)
        throw std::runtime_error("Validation error: Nx must be a positive integer.");
    if (cfg.grid.t_final <= 0)
        throw std::runtime_error("Validation error: t_final must be positive.");
    if (cfg.grid.CFL <= 0 || cfg.grid.CFL > 1.0)
        throw std::runtime_error("Validation error: CFL must be in the interval (0, 1].");
    if (cfg.phys.gamma <= 1.0)
        throw std::runtime_error("Validation error: gamma (ratio of specific heats) must be greater than 1.");
    if (cfg.phys.left.rho <= 0 || cfg.phys.right.rho <= 0)
        throw std::runtime_error("Validation error: density (rho) must be positive.");
    if (cfg.phys.left.p <= 0 || cfg.phys.right.p <= 0)
        throw std::runtime_error("Validation error: pressure (p) must be positive.");


    return cfg;
}