#include "test_cases.h"

#include <stdexcept>

void setup_test_case(Config& cfg) {
    if (cfg.phys.test_case == 0) {
        return;
    }
    switch (cfg.phys.test_case) {
    case 1:
        
        cfg.phys.left.rho = 1.0;    cfg.phys.left.u = 0.75;    cfg.phys.left.p = 1.0;
        cfg.phys.right.rho = 0.125;  cfg.phys.right.u = 0.0;   cfg.phys.right.p = 0.1;
        
        
        cfg.grid.x_diaphragm = 0.3;
        cfg.grid.t_final = 0.2;
        break;

    case 2: 
        
        cfg.phys.left.rho = 1.0;    cfg.phys.left.u = -2.0;   cfg.phys.left.p = 0.4;
        cfg.phys.right.rho = 1.0;    cfg.phys.right.u = 2.0;   cfg.phys.right.p = 0.4;
        
        
        cfg.grid.x_diaphragm = 0.5;
        cfg.grid.t_final = 0.15;
        break;

    case 3: 
        
        cfg.phys.left.rho = 1.0;    cfg.phys.left.u = 0.0;    cfg.phys.left.p = 1000.0;
        cfg.phys.right.rho = 1.0;    cfg.phys.right.u = 0.0;   cfg.phys.right.p = 0.01;
        
        
        cfg.grid.x_diaphragm = 0.5;
        cfg.grid.t_final = 0.012;
        break;

    case 4: 
        
        cfg.phys.left.rho = 5.99924; cfg.phys.left.u = 19.5975;  cfg.phys.left.p = 460.894;
        cfg.phys.right.rho = 5.99242; cfg.phys.right.u = -6.19633; cfg.phys.right.p = 46.0950;
        
        
        cfg.grid.x_diaphragm = 0.4;
        cfg.grid.t_final = 0.035;
        break;

    case 5: 
        
        cfg.phys.left.rho = 1.0;    cfg.phys.left.u = -19.59745; cfg.phys.left.p = 1000.0;
        cfg.phys.right.rho = 1.0;    cfg.phys.right.u = -19.59745; cfg.phys.right.p = 0.01;
        
        
        cfg.grid.x_diaphragm = 0.8;
        cfg.grid.t_final = 0.012;
        break;

    default:
        throw std::runtime_error("Unknown test_case number provided in config.txt. Available: 1-5.");
    }
}