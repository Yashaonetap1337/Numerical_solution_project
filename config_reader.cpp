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
                else {
                    throw std::runtime_error("Unknown method parameter: " + key);
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