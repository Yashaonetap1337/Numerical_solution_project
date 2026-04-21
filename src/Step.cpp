#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

// Параметры PISO
const int N_CORR = 2;         // Количество корректоров

struct PISO_Solver {
    int NX, NY;
    double LX, LY, dx, dy, nu, RHO, RE, U_LID;
    double DT, T_END;
    string test_name;

    vector<double> u, v, p, u_n, v_n;
    vector<double> auP, avP;
    vector<bool> blocked_cells;

    PISO_Solver(string test) : test_name(test), RHO(1.0), RE(100.0), U_LID(1.0) {
        if (test == "cavity") {
            NX = 31; NY = 31;
            LX = 1.0; LY = 1.0;
            DT = 0.005;
            T_END = 5.0;
        }
        else if (test == "taylor_green") {
            NX = 32; NY = 32;
            LX = 1.0; LY = 1.0;
            RE = 10.0;
            DT = 0.005;
            T_END = 0.5;
        }
        else if (test == "step") {
            NX = 300; NY = 40;
            LX = 15.0; LY = 1.0;
            DT = 0.01;
            T_END = 4.0;
        }

        dx = LX / NX; dy = LY / NY;
        nu = U_LID * (test == "step" ? 0.5 : 1.0) / RE;

        u = vector<double>((NX + 1) * NY, 0);
        v = vector<double>(NX * (NY + 1), 0);
        p = vector<double>(NX * NY, 0);
        u_n = vector<double>((NX + 1) * NY, 0);
        v_n = vector<double>(NX * (NY + 1), 0);
        auP = vector<double>((NX + 1) * NY, 0);
        avP = vector<double>(NX * (NY + 1), 0);
        blocked_cells = vector<bool>(NX * NY, false);

        initializeTest();
    }

    void initializeTest() {
        if (test_name == "taylor_green") {
            const double pi = M_PI;
            for (int i = 0; i <= NX; i++) {
                for (int j = 0; j < NY; j++) {
                    double x = i * dx;
                    double y = (j + 0.5) * dy;
                    u[i * NY + j] = -cos(pi * x) * sin(pi * y);
                }
            }
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j <= NY; j++) {
                    double x = (i + 0.5) * dx;
                    double y = j * dy;
                    v[i * (NY + 1) + j] = sin(pi * x) * cos(pi * y);
                }
            }
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    double x = (i + 0.5) * dx;
                    double y = (j + 0.5) * dy;
                    p[i * NY + j] = -0.25 * (cos(2 * pi * x) + cos(2 * pi * y));
                }
            }
        }
        else if (test_name == "step") {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    double x = (i + 0.5) * dx;
                    double y = (j + 0.5) * dy;
                    if (x < 2.5 && y < 0.5) {
                        blocked_cells[i * NY + j] = true;
                    }
                }
            }
        }
    }

    bool isBlocked(int i, int j) {
        if (test_name != "step") return false;
        if (i < 0 || i >= NX || j < 0 || j >= NY) return false;
        return blocked_cells[i * NY + j];
    }

    void applyBoundaryConditions(double t) {
        if (test_name == "cavity") {
            for (int i = 0; i <= NX; i++) {
                for (int j = 0; j < NY; j++) {
                    if (i == 0 || i == NX) u[i * NY + j] = 0;
                    if (j == 0) u[i * NY + j] = 0;
                    if (j == NY - 1) u[i * NY + j] = U_LID;
                }
            }
            for (int i = 0; i < NX; i++) {
                v[i * (NY + 1) + 0] = 0;
                v[i * (NY + 1) + NY] = 0;
            }
        }
        else if (test_name == "taylor_green") {
            // Граничные условия из точного решения
            const double pi = M_PI;
            double decay = exp(-2 * pi * pi * nu * t);

            for (int j = 0; j < NY; j++) {
                double y = (j + 0.5) * dy;
                u[0 * NY + j] = -cos(0) * sin(pi * y) * decay;
                u[NX * NY + j] = -cos(pi * LX) * sin(pi * y) * decay;
            }
            for (int i = 0; i < NX; i++) {
                double x = (i + 0.5) * dx;
                v[i * (NY + 1) + 0] = sin(pi * x) * cos(0) * decay;
                v[i * (NY + 1) + NY] = sin(pi * x) * cos(pi * LY) * decay;
            }
        }
        else if (test_name == "step") {
            double U_avg = 1.0;
            double h = 0.5;
            double H = 1.0;
            //CHECK: STEP_INLET
            for (int j = 0; j < NY; j++) {
                double y = (j + 0.5) * dy;
                if (y > h) {
                    u[0 * NY + j] = 6.0 * U_avg * (y - h) * (H - y) / pow(H - h, 2);
                }
                else {
                    u[0 * NY + j] = 0;
                }
            }

            for (int j = 0; j < NY; j++) {
                u[NX * NY + j] = u[(NX - 1) * NY + j];
            }
            for (int j = 0; j <= NY; j++) {
                v[(NX - 1) * (NY + 1) + j] = v[(NX - 2) * (NY + 1) + j];
            }

            for (int i = 0; i < NX; i++) {
                v[i * (NY + 1) + 0] = 0;
                v[i * (NY + 1) + NY] = 0;
            }
            //CHECK: STEP_BLOCK
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    if (isBlocked(i, j)) {
                        u[i * NY + j] = 0;
                        u[(i + 1) * NY + j] = 0;
                        v[i * (NY + 1) + j] = 0;
                        v[i * (NY + 1) + j + 1] = 0;
                    }
                }
            }
        }
    }

    void momentumPredictor() {
        u_n = u; v_n = v;

        // U*
        for (int i = 1; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (isBlocked(i, j) || isBlocked(i - 1, j)) continue;

                double uc = u[i * NY + j], uE = u[(i + 1) * NY + j], uW = u[(i - 1) * NY + j];
                double ue = 0.5 * (uc + uE), uw = 0.5 * (uc + uW);
                double vn = 0.5 * (v[(i - 1) * (NY + 1) + j + 1] + v[i * (NY + 1) + j + 1]);
                double vs = 0.5 * (v[(i - 1) * (NY + 1) + j] + v[i * (NY + 1) + j]);

                double Fe = RHO * ue * dy, Fw = RHO * uw * dy, Fn = RHO * vn * dx, Fs = RHO * vs * dx;
                double aE = max(-Fe, 0.0) + nu * dy / dx;
                double aW = max(Fw, 0.0) + nu * dy / dx;
                double aN = (j == NY - 1) ? max(-Fn, 0.0) + nu * dx / (0.5 * dy) : max(-Fn, 0.0) + nu * dx / dy;
                double aS = (j == 0) ? max(Fs, 0.0) + nu * dx / (0.5 * dy) : max(Fs, 0.0) + nu * dx / dy;

                double uNv = (j == NY - 1) ? (test_name == "cavity" ? 2 * U_LID - uc : -uc) : u[i * NY + j + 1];
                double uSv = (j == 0) ? -uc : u[i * NY + j - 1];

                double aP_c = aE + aW + aN + aS + RHO * dx * dy / DT;
                double H = aE * uE + aW * uW + aN * uNv + aS * uSv + (RHO * dx * dy / DT) * u_n[i * NY + j];

                double dP = (p[i * NY + j] - p[(i - 1) * NY + j]) * dy;
                u[i * NY + j] = (H - dP) / aP_c;
                auP[i * NY + j] = aP_c;
            }
        }

        // V*
        for (int i = 0; i < NX; i++) {
            for (int j = 1; j < NY; j++) {
                if (isBlocked(i, j) || isBlocked(i, j - 1)) continue;

                double vc = v[i * (NY + 1) + j], vN = v[i * (NY + 1) + j + 1], vS = v[i * (NY + 1) + j - 1];
                double vn = 0.5 * (vc + vN), vsf = 0.5 * (vc + vS);
                double ue = 0.5 * (u[(i + 1) * NY + j - 1] + u[(i + 1) * NY + j]);
                double uw = 0.5 * (u[i * NY + j - 1] + u[i * NY + j]);

                double Fe = RHO * ue * dy, Fw = RHO * uw * dy, Fn = RHO * vn * dx, Fs = RHO * vsf * dx;
                double aN = max(-Fn, 0.0) + nu * dx / dy;
                double aS = max(Fs, 0.0) + nu * dx / dy;
                double aE = (i == NX - 1) ? max(-Fe, 0.0) + nu * dy / (test_name == "step" ? dx : 0.5 * dx) : max(-Fe, 0.0) + nu * dy / dx;
                double aW = (i == 0) ? max(Fw, 0.0) + nu * dy / (0.5 * dx) : max(Fw, 0.0) + nu * dy / dx;

                double vEv = (i == NX - 1) ? (test_name == "step" ? vc : -vc) : v[(i + 1) * (NY + 1) + j];
                double vWv = (i == 0) ? -vc : v[(i - 1) * (NY + 1) + j];

                double aP_c = aE + aW + aN + aS + RHO * dx * dy / DT;
                double H = aE * vEv + aW * vWv + aN * vN + aS * vS + (RHO * dx * dy / DT) * v_n[i * (NY + 1) + j];

                double dP = (p[i * NY + j] - p[i * NY + j - 1]) * dx;
                v[i * (NY + 1) + j] = (H - dP) / aP_c;
                avP[i * (NY + 1) + j] = aP_c;
            }
        }
    }

    void solvePressure() {
        vector<double> p_corr(NX * NY, 0);
        for (int step = 0; step < 100; step++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    if (isBlocked(i, j)) continue;
                    if (i == 0 && j == 0) continue;

                    double aE = (i < NX - 1 && !isBlocked(i + 1, j) && auP[(i + 1) * NY + j] > 1e-12) ? (dy * dy / auP[(i + 1) * NY + j]) : 0;
                    double aW = (i > 0 && !isBlocked(i - 1, j) && auP[i * NY + j] > 1e-12) ? (dy * dy / auP[i * NY + j]) : 0;
                    double aN = (j < NY - 1 && !isBlocked(i, j + 1) && avP[i * (NY + 1) + j + 1] > 1e-12) ? (dx * dx / avP[i * (NY + 1) + j + 1]) : 0;
                    double aS = (j > 0 && !isBlocked(i, j - 1) && avP[i * (NY + 1) + j] > 1e-12) ? (dx * dx / avP[i * (NY + 1) + j]) : 0;
                    double aP = aE + aW + aN + aS;

                    if (aP < 1e-12) continue;

                    double b = -((u[(i + 1) * NY + j] - u[i * NY + j]) * dy + (v[i * (NY + 1) + j + 1] - v[i * (NY + 1) + j]) * dx);

                    double pE = (i < NX - 1) ? p_corr[(i + 1) * NY + j] : 0;
                    double pW = (i > 0) ? p_corr[(i - 1) * NY + j] : 0;
                    double pN = (j < NY - 1) ? p_corr[i * NY + j + 1] : 0;
                    double pS = (j > 0) ? p_corr[i * NY + j - 1] : 0;

                    p_corr[i * NY + j] = (aE * pE + aW * pW + aN * pN + aS * pS + b) / aP;
                }
            }
        }

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (!isBlocked(i, j)) p[i * NY + j] += p_corr[i * NY + j];
            }
        }
        for (int i = 1; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (!isBlocked(i, j) && !isBlocked(i - 1, j) && auP[i * NY + j] > 1e-12) {
                    u[i * NY + j] -= (dy / auP[i * NY + j]) * (p_corr[i * NY + j] - p_corr[(i - 1) * NY + j]);
                }
            }
        }
        for (int i = 0; i < NX; i++) {
            for (int j = 1; j < NY; j++) {
                if (!isBlocked(i, j) && !isBlocked(i, j - 1) && avP[i * (NY + 1) + j] > 1e-12) {
                    v[i * (NY + 1) + j] -= (dx / avP[i * (NY + 1) + j]) * (p_corr[i * NY + j] - p_corr[i * NY + j - 1]);
                }
            }
        }
    }

    void run() {
        cout << "Starting PISO solver for test: " << test_name << endl;
        cout << "Grid: " << NX << "x" << NY << ", Re = " << RE << ", DT = " << DT << endl;

        double t = 0;
        int step = 0;
        while (t < T_END) {
            applyBoundaryConditions(t);
            momentumPredictor();
            for (int c = 0; c < N_CORR; c++) {
                solvePressure();
            }
            applyBoundaryConditions(t + DT);

            t += DT;
            step++;
            if (step % 20 == 0) cout << "Time: " << fixed << setprecision(3) << t << " s" << endl;
        }
        saveToCSV();
    }

    void saveToCSV() {
        string filename = "piso_" + test_name + ".csv";
        ofstream out(filename);
        out << "x,y,u,v,p\n";
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (isBlocked(i, j)) continue;

                double uc = 0.5 * (u[i * NY + j] + u[(i + 1) * NY + j]);
                double vc = 0.5 * (v[i * (NY + 1) + j] + v[i * (NY + 1) + j + 1]);
                double pc = p[i * NY + j];
                out << (i + 0.5) * dx << "," << (j + 0.5) * dy << "," << uc << "," << vc << "," << pc << "\n";
            }
        }
        out.close();
        cout << "Data saved to " << filename << endl;

        if (test_name == "step") {
            cout << "\nRecirculation length analysis:" << endl;
            double h = 0.5;
            double step_end = 2.5;
            for (int i = 0; i < NX; i++) {
                double x = (i + 0.5) * dx;
                if (x > step_end) {
                    double u_bottom = 0.5 * (u[i * NY + 0] + u[(i + 1) * NY + 0]);
                    if (u_bottom > 0 && i > 0) {
                        double x_prev = (i - 0.5) * dx;
                        double u_prev = 0.5 * (u[(i - 1) * NY + 0] + u[i * NY + 0]);
                        if (u_prev <= 0) {
                            double x_r = x - step_end;
                            cout << "Reattachment point: x_r = " << x_r << " (" << x_r / h << " * h)" << endl;
                            break;
                        }
                    }
                }
            }
        }
    }
};

int main() {
    vector<string> tests = { "step" };

    for (const auto& test : tests) {
        cout << "\n" << string(60, '=') << endl;
        PISO_Solver solver(test);
        solver.run();
        cout << string(60, '=') << "\n" << endl;
    }

    return 0;
}
