#include "MaderSolver.h"
#include <cmath>
#include <algorithm>
#include <vector>

OmegaField make_omega_field(const Grid& grid) {
    size_t NX = (size_t)(grid.Nx + 2 * grid.num_fict);
    size_t NY = (size_t)(grid.Ny + 2 * grid.num_fict);
    return OmegaField(NX, std::vector<double>(NY, 0.0));
}

void MaderTimeStep(const Grid& grid, const Config& cfg, const std::vector<std::vector<State>>& W,
    std::vector<std::vector<State>>& W_new, OmegaField& omega, double dt) {

    const int ng = grid.num_fict;
    const int Nx = grid.Nx; const int Ny = grid.Ny;
    const size_t NX = (size_t)(Nx + 2 * ng);
    const size_t NY = (size_t)(Ny + 2 * ng);
    const double dx = grid.dx; const double dy = grid.dy;
    const double gm1 = cfg.phys.gamma - 1.0;

    auto A = [&](double v = 0.0) { return std::vector<std::vector<double>>(NX, std::vector<double>(NY, v)); };
    auto P = A(), qx = A(), qy = A(), u_t = A(), v_t = A(), r_t = A(), rE_t = A();
    auto DM = A(), DE = A(), DPU = A(), DPV = A();

    // 1. Давление и вязкость
    for (size_t i = 0; i < NX; i++)
        for (size_t j = 0; j < NY; j++) P[i][j] = std::max(W[i][j].p, 1e-6);

    for (size_t i = 1; i < NX - 1; i++) {
        for (size_t j = 0; j < NY; j++) {
            double r = std::max(W[i][j].rho, 1e-6);
            double du = W[i + 1][j].u - W[i - 1][j].u;
            if (du < 0) qx[i][j] = cfg.viscosity_coeff * r * du * du;
        }
    }

    // 2. Лагранжев шаг
    for (int i = ng; i < Nx + ng; i++) {
        for (int j = ng; j < Ny + ng; j++) {
            double r = std::max(W[i][j].rho, 1e-6);
            // Ускорение
            double du = -(dt / (r * dx * 2.0)) * (P[i + 1][j] - P[i - 1][j] + qx[i + 1][j] - qx[i - 1][j]);
            u_t[i][j] = W[i][j].u + du;
            v_t[i][j] = W[i][j].v; // В 1D v не меняется

            // Плотность
            double div = (W[i + 1][j].u - W[i - 1][j].u) / (2.0 * dx);
            r_t[i][j] = r / (1.0 + dt * div + 1e-12);

            const double Q_HEAT = 5.0;

            // ЭНЕРГИЯ (Полная удельная энергия): E = I + 0.5*u^2
            double I_old = P[i][j] / (gm1 * r);
            double I_new = I_old * std::pow(r_t[i][j] / r, gm1); // Адиабатическое приближение
            rE_t[i][j] = r_t[i][j] * (I_new + 0.5 * (u_t[i][j] * u_t[i][j]));
        }
    }

    // 3. Перенос (Remap)
    for (int i = ng; i <= Nx + ng; i++) {
        for (int j = ng; j < Ny + ng; j++) {
            double alpha = 0.5 * (u_t[i - 1][j] + u_t[i][j]) * dt / dx;
            alpha = std::max(-0.2, std::min(0.2, alpha));
            int d = (alpha >= 0) ? i - 1 : i;
            int a = (alpha >= 0) ? i : i - 1;
            double flow = r_t[d][j] * std::abs(alpha);
            DM[a][j] += flow;  DM[d][j] -= flow;
            DE[a][j] += (rE_t[d][j] / r_t[d][j]) * flow; DE[d][j] -= (rE_t[d][j] / r_t[d][j]) * flow;
            DPU[a][j] += u_t[d][j] * flow; DPU[d][j] -= u_t[d][j] * flow;
        }
    }

    // 4. Финализация
    W_new = W;
    for (int i = ng; i < Nx + ng; i++) {
        for (int j = ng; j < Ny + ng; j++) {
            double rn = std::max(1e-6, r_t[i][j] + DM[i][j]);
            double un = (r_t[i][j] * u_t[i][j] + DPU[i][j]) / rn;
            double En = (rE_t[i][j] + DE[i][j]) / rn;
            double In = std::max(1e-6, En - 0.5 * un * un);

            W_new[i][j].rho = rn;
            W_new[i][j].u = un;
            W_new[i][j].v = 0.0;
            W_new[i][j].p = gm1 * rn * In;
        }
    }
}