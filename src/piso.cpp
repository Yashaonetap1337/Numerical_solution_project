#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

// Параметры сетки и физики
const int NX = 31;
const int NY = 31;
const double LX = 1.0, LY = 1.0;
const double RHO = 1.0, RE = 100.0, U_LID = 1.0;

// Параметры PISO
const double DT = 0.005;      // Шаг по времени (должен быть малым!)
const double T_END = 5.0;     // Общее время расчета
const int N_CORR = 2;         // Количество корректоров (стандарт для PISO - 2)

struct PISO_Solver {
    double dx, dy, nu;
    vector<double> u, v, p, u_n, v_n;
    vector<double> auP, avP;

    PISO_Solver() :
        dx(LX / NX), dy(LY / NY), nu(U_LID* LX / RE),
        u((NX + 1)* NY, 0), v(NX* (NY + 1), 0), p(NX* NY, 0),
        u_n((NX + 1)* NY, 0), v_n(NX* (NY + 1), 0),
        auP((NX + 1)* NY, 0), avP(NX* (NY + 1), 0) {}

    bool isBlocked(int i, int j) { return false; } // Для каверны нет препятствий

    // --- Предиктор Импульса ---
    void momentumPredictor() {
        u_n = u; v_n = v; // Сохраняем значения с предыдущего шага времени

        // Расчет U*
        for (int i = 1; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double uc = u[i * NY + j], uE = u[(i + 1) * NY + j], uW = u[(i - 1) * NY + j];
                double ue = 0.5 * (uc + uE), uw = 0.5 * (uc + uW);
                double vn = 0.5 * (v[(i - 1) * (NY + 1) + j + 1] + v[i * (NY + 1) + j + 1]);
                double vs = 0.5 * (v[(i - 1) * (NY + 1) + j] + v[i * (NY + 1) + j]);

                double Fe = RHO * ue * dy, Fw = RHO * uw * dy, Fn = RHO * vn * dx, Fs = RHO * vs * dx;
                double aE = max(-Fe, 0.0) + nu * dy / dx;
                double aW = max(Fw, 0.0) + nu * dy / dx;
                double aN = (j == NY - 1) ? max(-Fn, 0.0) + nu * dx / (0.5 * dy) : max(-Fn, 0.0) + nu * dx / dy;
                double aS = (j == 0) ? max(Fs, 0.0) + nu * dx / (0.5 * dy) : max(Fs, 0.0) + nu * dx / dy;

                double uNv = (j == NY - 1) ? 2 * U_LID - uc : u[i * NY + j + 1];
                double uSv = (j == 0) ? -uc : u[i * NY + j - 1];

                double aP_c = aE + aW + aN + aS + RHO * dx * dy / DT; // Добавлен временной член
                double H = aE * uE + aW * uW + aN * uNv + aS * uSv + (RHO * dx * dy / DT) * u_n[i * NY + j];

                double dP = (p[i * NY + j] - p[(i - 1) * NY + j]) * dy;
                u[i * NY + j] = (H - dP) / aP_c;
                auP[i * NY + j] = aP_c;
            }
        }

        // Расчет V* (аналогично)
        for (int i = 0; i < NX; i++) {
            for (int j = 1; j < NY; j++) {
                double vc = v[i * (NY + 1) + j], vN = v[i * (NY + 1) + j + 1], vS = v[i * (NY + 1) + j - 1];
                double vn = 0.5 * (vc + vN), vsf = 0.5 * (vc + vS);
                double ue = 0.5 * (u[(i + 1) * NY + j - 1] + u[(i + 1) * NY + j]);
                double uw = 0.5 * (u[i * NY + j - 1] + u[i * NY + j]);

                double Fe = RHO * ue * dy, Fw = RHO * uw * dy, Fn = RHO * vn * dx, Fs = RHO * vsf * dx;
                double aN = max(-Fn, 0.0) + nu * dx / dy;
                double aS = max(Fs, 0.0) + nu * dx / dy;
                double aE = (i == NX - 1) ? max(-Fe, 0.0) + nu * dy / (0.5 * dx) : max(-Fe, 0.0) + nu * dy / dx;
                double aW = (i == 0) ? max(Fw, 0.0) + nu * dy / (0.5 * dx) : max(Fw, 0.0) + nu * dy / dx;

                double vEv = (i == NX - 1) ? -vc : v[(i + 1) * (NY + 1) + j];
                double vWv = (i == 0) ? -vc : v[(i - 1) * (NY + 1) + j];

                double aP_c = aE + aW + aN + aS + RHO * dx * dy / DT;
                double H = aE * vEv + aW * vWv + aN * vN + aS * vS + (RHO * dx * dy / DT) * v_n[i * (NY + 1) + j];

                double dP = (p[i * NY + j] - p[i * NY + j - 1]) * dx;
                v[i * (NY + 1) + j] = (H - dP) / aP_c;
                avP[i * (NY + 1) + j] = aP_c;
            }
        }
    }

    // --- Решатель уравнения Пуассона ---
    void solvePressure() {
        vector<double> p_corr(NX * NY, 0);
        for (int step = 0; step < 100; step++) { // Итерации Гаусса-Зейделя
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    if (i == 0 && j == 0) continue; // Точка фиксации давления

                    double aE = (i < NX - 1) ? (dy * dy / auP[(i + 1) * NY + j]) : 0;
                    double aW = (i > 0) ? (dy * dy / auP[i * NY + j]) : 0;
                    double aN = (j < NY - 1) ? (dx * dx / avP[i * (NY + 1) + j + 1]) : 0;
                    double aS = (j > 0) ? (dx * dx / avP[i * (NY + 1) + j]) : 0;
                    double aP = aE + aW + aN + aS;

                    double b = -((u[(i + 1) * NY + j] - u[i * NY + j]) * dy + (v[i * (NY + 1) + j + 1] - v[i * (NY + 1) + j]) * dx);

                    double pE = (i < NX - 1) ? p_corr[(i + 1) * NY + j] : 0;
                    double pW = (i > 0) ? p_corr[(i - 1) * NY + j] : 0;
                    double pN = (j < NY - 1) ? p_corr[i * NY + j + 1] : 0;
                    double pS = (j > 0) ? p_corr[i * NY + j - 1] : 0;

                    p_corr[i * NY + j] = (aE * pE + aW * pW + aN * pN + aS * pS + b) / aP;
                }
            }
        }

        // Коррекция полей
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) p[i * NY + j] += p_corr[i * NY + j];
        }
        for (int i = 1; i < NX; i++) {
            for (int j = 0; j < NY; j++) u[i * NY + j] -= (dy / auP[i * NY + j]) * (p_corr[i * NY + j] - p_corr[(i - 1) * NY + j]);
        }
        for (int i = 0; i < NX; i++) {
            for (int j = 1; j < NY; j++) v[i * (NY + 1) + j] -= (dx / avP[i * (NY + 1) + j]) * (p_corr[i * NY + j] - p_corr[i * NY + j - 1]);
        }
    }

    void run() {
        double t = 0;
        int step = 0;
        //CHECK: PISO_loop
        while (t < T_END) {
            momentumPredictor();
            for (int c = 0; c < N_CORR; c++) {
                solvePressure();
            }
            t += DT;
            step++;
            if (step % 20 == 0) cout << "Time: " << t << " s" << endl;
        }
        saveToCSV();
    }

    void saveToCSV() {
        ofstream out("cavity_piso.csv");
        out << "x,y,u,v\n";
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double uc = 0.5 * (u[i * NY + j] + u[(i + 1) * NY + j]);
                double vc = 0.5 * (v[i * (NY + 1) + j] + v[i * (NY + 1) + j + 1]);
                out << (i + 0.5) * dx << "," << (j + 0.5) * dy << "," << uc << "," << vc << "\n";
            }
        }
        out.close();
    }
};

int main() {
    PISO_Solver solver;
    solver.run();
    return 0;
}