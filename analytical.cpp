#include "types.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>


struct AnalyticalSolution {
    std::vector<double> x, rho, p, u, e;
};


AnalyticalSolution get_analytical_solution(int test_num) {
    std::string filename = "toro_" + std::to_string(test_num) + "_exact.txt";
    AnalyticalSolution solution;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        double val_x, val_rho, val_p, val_u, val_e, unused;
        if (ss >> val_x >> val_rho >> val_p >> val_u >> val_e >> unused) {
            solution.x.push_back(val_x);
            solution.rho.push_back(val_rho);
            solution.p.push_back(val_p);
            solution.u.push_back(val_u);
            solution.e.push_back(val_e);
        }
    }
    file.close();
    return solution;
}

void save_analytical_solution(int test_num, const std::string& outputFilename) {
    AnalyticalSolution sol = get_analytical_solution(test_num);

    std::ofstream out(outputFilename);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to create file: " + outputFilename);
    }

    // запись CSV
    out << "x,rho,p,u,e\n";
    for (size_t i = 0; i < sol.x.size(); ++i) {
        out << sol.x[i] << ","
            << sol.rho[i] << ","
            << sol.p[i] << ","
            << sol.u[i] << ","
            << sol.e[i] << "\n";
    }
    out.close();
}