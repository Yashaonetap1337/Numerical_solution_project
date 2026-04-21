#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <fstream>

using namespace std;

// Параметры расчетной области и среды
const int NX = 41;          // Количество ячеек по X
const int NY = 41;          // Количество ячеек по Y
const double LX = 1.0;      // Длина области
const double LY = 1.0;      // Высота области
const double RHO = 1.0;     // Плотность
const double RE = 100.0;    // Число Рейнольдса
const double U_LID = 1.0;   // Скорость крышки

// Параметры SIMPLE
const double ALPHA_U = 0.7; // Релаксация для скорости
const double ALPHA_P = 0.3; // Релаксация для давления
const double EPS = 1e-6;    // Критерий сходимости
const int MAX_ITER = 5000;  // Макс. итераций

struct Solver {
    double dx, dy, nu;
    // CHECK: staggered_grid: u(NX+1, NY), v(NX, NY+1), p(NX, NY)
    vector<vector<double>> u, v, p, p_corr;
    vector<vector<double>> u_old, v_old;
    vector<vector<double>> auP, avP; // Эффективные диагонали
    
    Solver() :
        dx(LX / NX), dy(LY / NY), nu(U_LID* LX / RE),
        u(NX + 1, vector<double>(NY, 0.0)),
        v(NX, vector<double>(NY + 1, 0.0)),
        p(NX, vector<double>(NY, 0.0)),
        p_corr(NX, vector<double>(NY, 0.0)),
        u_old(NX + 1, vector<double>(NY, 0.0)),
        v_old(NX, vector<double>(NY + 1, 0.0)),
        auP(NX + 1, vector<double>(NY, 0.0)),
        avP(NX, vector<double>(NY + 1, 0.0)) {}
    //CHECK: simple_loop
    void solve() {
        for (int iter = 1; iter <= MAX_ITER; ++iter) {
            u_old = u;
            v_old = v;

            // 1. Решение уравнений импульса (Предиктор)
            solveMomentum();

            // 2. Решение уравнения поправки давления (Пуассон)
            solvePressureCorrection();

            // 3. Коррекция скоростей и давления
            double max_res = updateFlow();

            if (iter % 100 == 0 || iter == 1) {
                cout << "Iter: " << iter << " | Max Residual (dU): " << scientific << max_res << endl;
            }

            if (max_res < EPS) {
                cout << "Converged at iteration " << iter << endl;
                break;
            }
        }
    }

    void solveMomentum() {
        // CHECK: predictor_u
        for (int i = 1; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                double ue = 0.5 * (u[i][j] + u[i + 1][j]);
                double uw = 0.5 * (u[i][j] + u[i - 1][j]);
                double vn = 0.5 * (v[i - 1][j + 1] + v[i][j + 1]);
                double vs = 0.5 * (v[i - 1][j] + v[i][j]);

                double Fe = RHO * ue * dy;
                double Fw = RHO * uw * dy;
                double Fn = RHO * vn * dx;
                double Fs = RHO * vs * dx;

                double aE = max(-Fe, 0.0) + nu * dy / dx;
                double aW = max(Fw, 0.0) + nu * dy / dx;
                double aN, aS, uN_ghost, uS_ghost;

                // Граничные условия для U по Y (стенки)
                if (j == NY - 1) { // Крышка
                    aN = max(-Fn, 0.0) + nu * dx / (dy * 0.5);
                    uN_ghost = 2.0 * U_LID - u[i][j];
                }
                else {
                    aN = max(-Fn, 0.0) + nu * dx / dy;
                    uN_ghost = u[i][j + 1];
                }

                if (j == 0) { // Дно
                    aS = max(Fs, 0.0) + nu * dx / (dy * 0.5);
                    uS_ghost = -u[i][j];
                }
                else {
                    aS = max(Fs, 0.0) + nu * dx / dy;
                    uS_ghost = u[i][j - 1];
                }

                double aP = aE + aW + aN + aS;
                double H = aE * u[i + 1][j] + aW * u[i - 1][j] + aN * uN_ghost + aS * uS_ghost;

                double dP = (p[i][j] - p[i - 1][j]) * dy;

                // Под-релаксация (формула 21)
                u[i][j] = (1.0 - ALPHA_U) * u_old[i][j] + (ALPHA_U / aP) * (H - dP);
                auP[i][j] = aP / ALPHA_U; // Эффективная диагональ
            }
        }

        // // CHECK: predictor_v
        for (int i = 0; i < NX; ++i) {
            for (int j = 1; j < NY; ++j) {
                double ue = 0.5 * (u[i + 1][j] + u[i + 1][j - 1]);
                double uw = 0.5 * (u[i][j] + u[i][j - 1]);
                double vn = 0.5 * (v[i][j] + v[i][j + 1]);
                double vs = 0.5 * (v[i][j] + v[i][j - 1]);

                double Fe = RHO * ue * dy;
                double Fw = RHO * uw * dy;
                double Fn = RHO * vn * dx;
                double Fs = RHO * vs * dx;

                double aN = max(-Fn, 0.0) + nu * dx / dy;
                double aS = max(Fs, 0.0) + nu * dx / dy;
                double aE, aW, vE_ghost, vW_ghost;

                if (i == NX - 1) { // Правая стенка
                    aE = max(-Fe, 0.0) + nu * dy / (dx * 0.5);
                    vE_ghost = -v[i][j];
                }
                else {
                    aE = max(-Fe, 0.0) + nu * dy / dx;
                    vE_ghost = v[i + 1][j];
                }

                if (i == 0) { // Левая стенка
                    aW = max(Fw, 0.0) + nu * dy / (dx * 0.5);
                    vW_ghost = -v[i][j];
                }
                else {
                    aW = max(Fw, 0.0) + nu * dy / dx;
                    vW_ghost = v[i - 1][j];
                }

                double aP = aE + aW + aN + aS;
                double H = aE * vE_ghost + aW * vW_ghost + aN * v[i][j + 1] + aS * v[i][j - 1];

                double dP = (p[i][j] - p[i][j - 1]) * dx;

                v[i][j] = (1.0 - ALPHA_U) * v_old[i][j] + (ALPHA_U / aP) * (H - dP);
                avP[i][j] = aP / ALPHA_U;
            }
        }
    }

    void solvePressureCorrection() {
        fill(p_corr.begin(), p_corr.end(), vector<double>(NY, 0.0));

        for (int solve_it = 0; solve_it < 100; ++solve_it) { // Итерации решателя (Гаусс-Зейдель)
            for (int i = 0; i < NX; ++i) {
                for (int j = 0; j < NY; ++j) {
                    if (i == 0 && j == 0) continue; // Фиксация давления p'[0,0] = 0

                    double aE = (i < NX - 1) ? (dy * dy / auP[i + 1][j]) : 0;
                    double aW = (i > 0) ? (dy * dy / auP[i][j]) : 0;
                    double aN = (j < NY - 1) ? (dx * dx / avP[i][j + 1]) : 0;
                    double aS = (j > 0) ? (dx * dx / avP[i][j]) : 0;
                    double aP = aE + aW + aN + aS;
                    //CHECK: poisson
                    double b = -((u[i + 1][j] - u[i][j]) * dy + (v[i][j + 1] - v[i][j]) * dx);

                    double pE = (i < NX - 1) ? p_corr[i + 1][j] : 0;
                    double pW = (i > 0) ? p_corr[i - 1][j] : 0;
                    double pN = (j < NY - 1) ? p_corr[i][j + 1] : 0;
                    double pS = (j > 0) ? p_corr[i][j - 1] : 0;
                    
                    p_corr[i][j] = (aE * pE + aW * pW + aN * pN + aS * pS + b) / aP;
                }
            }
        }
    }

    double updateFlow() {
        double max_du = 0;

        //CHECK: correction
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                p[i][j] += ALPHA_P * p_corr[i][j];
            }
        }

        // Коррекция U
        for (int i = 1; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                double correction = (dy / auP[i][j]) * (p_corr[i][j] - p_corr[i - 1][j]);
                u[i][j] -= correction;
                max_du = max(max_du, abs(u[i][j] - u_old[i][j]));
            }
        }

        // Коррекция V
        for (int i = 0; i < NX; ++i) {
            for (int j = 1; j < NY; ++j) {
                double correction = (dx / avP[i][j]) * (p_corr[i][j] - p_corr[i][j - 1]);
                v[i][j] -= correction;
                max_du = max(max_du, abs(v[i][j] - v_old[i][j]));
            }
        }

        return max_du;
    }
};

int main() {
    cout << "Starting SIMPLE Solver for Lid-Driven Cavity (Re=" << RE << ")..." << endl;
    Solver solver;
    solver.solve();

    ofstream out("cavity_flow.csv");
    out << "x,y,u,v" << endl;
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double x = (i + 0.5) * solver.dx;
            double y = (j + 0.5) * solver.dy;

            // Интерполяция со смещенных граней в центр ячейки
            double uc = 0.5 * (solver.u[i][j] + solver.u[i + 1][j]);
            double vc = 0.5 * (solver.v[i][j] + solver.v[i][j + 1]);

            out << x << "," << y << "," << uc << "," << vc << endl;
        }
    }
    out.close();
    cout << "Data saved to cavity_flow.csv" << endl;

    // Вывод центральных профилей для проверки
    cout << "\nU-velocity profile at x = 0.5:" << endl;
    int mid_x = NX / 2;
    for (int j = NY - 1; j >= 0; --j) {
        double y = (j + 0.5) * (LY / NY);
        cout << fixed << setprecision(4) << "y: " << y << " | u: " << solver.u[mid_x][j] << endl;
    }


    return 0;
}