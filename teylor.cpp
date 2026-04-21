#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const int NX = 41, NY = 41;
const double LX = 1.0, LY = 1.0;
const double NU = 0.01;
const double RHO = 1.0;
const double DT = 0.001;
const double T_END = 10.0;

struct SolverTG_CSV {
    double dx, dy;
    vector<double> u, v, p;
    ofstream decay_file;

    SolverTG_CSV() : dx(LX / NX), dy(LY / NY),
        u((NX + 1)* NY, 0), v(NX* (NY + 1), 0), p(NX* NY, 0) {
        decay_file.open("decay.csv");
        decay_file << "t,max_u,analytical,error_l2\n"; // Заголовок CSV
    }

    // CHECK: TG_INIT
    double ex_u(double x, double y, double t) { return -cos(M_PI * x) * sin(M_PI * y) * exp(-2.0 * M_PI * M_PI * NU * t); }
    double ex_v(double x, double y, double t) { return sin(M_PI * x) * cos(M_PI * y) * exp(-2.0 * M_PI * M_PI * NU * t); }

    void save_field_csv(double t) {
        ofstream f("field.csv");
        f << "x,y,u,v,mag\n";
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double xc = (i + 0.5) * dx;
                double yc = (j + 0.5) * dy;
                // Усредняем скорость в центр ячейки
                double uc = 0.5 * (u[i * NY + j] + u[(i + 1) * NY + j]);
                double vc = 0.5 * (v[i * (NY + 1) + j] + v[i * (NY + 1) + j + 1]);
                f << xc << "," << yc << "," << uc << "," << vc << "," << sqrt(uc * uc + vc * vc) << "\n";
            }
        }
        f.close();
    }

    void run() {
        for (int i = 0; i <= NX; i++)
            for (int j = 0; j < NY; j++) u[i * NY + j] = ex_u(i * dx, (j + 0.5) * dy, 0);
        for (int i = 0; i < NX; i++)
            for (int j = 0; j <= NY; j++) v[i * (NY + 1) + j] = ex_v((i + 0.5) * dx, j * dy, 0);

        double t = 0;
        int step = 0;
        while (t < T_END) {
            t += DT; step++;


            // CHECK: TG_BC
            // Границы
            for (int j = 0; j < NY; j++) {
                u[0 * NY + j] = ex_u(0, (j + 0.5) * dy, t);
                u[NX * NY + j] = ex_u(LX, (j + 0.5) * dy, t);
            }
            for (int i = 0; i < NX; i++) {
                v[i * (NY + 1) + 0] = ex_v((i + 0.5) * dx, 0, t);
                v[i * (NY + 1) + NY] = ex_v((i + 0.5) * dx, LY, t);
            }

            // PISO (Предиктор + Коррекция)
            vector<double> u_s = u, v_s = v;
            for (int i = 1; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    double ut = (j == NY - 1) ? ex_u(i * dx, LY, t) : u[i * NY + j + 1];
                    double ub = (j == 0) ? ex_u(i * dx, 0, t) : u[i * NY + j - 1];
                    double lap = (u[(i + 1) * NY + j] - 2 * u[i * NY + j] + u[(i - 1) * NY + j]) / (dx * dx) + (ut - 2 * u[i * NY + j] + ub) / (dy * dy);
                    u_s[i * NY + j] = u[i * NY + j] + DT * (NU * lap - (p[i * NY + j] - p[(i - 1) * NY + j]) / dx);
                }
            }
            for (int i = 0; i < NX; i++) {
                for (int j = 1; j < NY; j++) {
                    double vr = (i == NX - 1) ? ex_v(LX, j * dy, t) : v[(i + 1) * (NY + 1) + j];
                    double vl = (i == 0) ? ex_v(0, j * dy, t) : v[(i - 1) * (NY + 1) + j];
                    double lap = (vr - 2 * v[i * (NY + 1) + j] + vl) / (dx * dx) + (v[i * (NY + 1) + j + 1] - 2 * v[i * (NY + 1) + j] + v[i * (NY + 1) + j - 1]) / (dy * dy);
                    v_s[i * (NY + 1) + j] = v[i * (NY + 1) + j] + DT * (NU * lap - (p[i * NY + j] - p[i * NY + j - 1]) / dy);
                }
            }
            // Давление
            for (int iter = 0; iter < 40; iter++) {
                for (int i = 0; i < NX; i++) {
                    for (int j = 0; j < NY; j++) {
                        double div = (u_s[(i + 1) * NY + j] - u_s[i * NY + j]) / dx + (v_s[i * (NY + 1) + j + 1] - v_s[i * (NY + 1) + j]) / dy;
                        double pE = (i == NX - 1) ? p[i * NY + j] : p[(i + 1) * NY + j];
                        double pW = (i == 0) ? p[i * NY + j] : p[(i - 1) * NY + j];
                        double pN = (j == NY - 1) ? p[i * NY + j] : p[i * NY + j + 1];
                        double pS = (j == 0) ? p[i * NY + j] : p[i * NY + j - 1];
                        p[i * NY + j] = 0.25 * (pE + pW + pN + pS - dx * dy * div * RHO / DT);
                    }
                }
            }
            // Обновление
            for (int i = 1; i < NX; i++)
                for (int j = 0; j < NY; j++) u[i * NY + j] = u_s[i * NY + j] - DT / RHO * (p[i * NY + j] - p[(i - 1) * NY + j]) / dx;
            for (int i = 0; i < NX; i++)
                for (int j = 1; j < NY; j++) v[i * (NY + 1) + j] = v_s[i * (NY + 1) + j] - DT / RHO * (p[i * NY + j] - p[i * NY + j - 1]) / dy;

            // Логирование затухания
            if (step % 100 == 0) {
                double max_u = 0, l2 = 0;
                for (int i = 0; i <= NX; i++) {
                    for (int j = 0; j < NY; j++) {
                        max_u = max(max_u, abs(u[i * NY + j]));
                        l2 += pow(u[i * NY + j] - ex_u(i * dx, (j + 0.5) * dy, t), 2);
                    }
                }
                decay_file << t << "," << max_u << "," << exp(-2 * M_PI * M_PI * NU * t) << "," << sqrt(l2 / (NX * NY)) << "\n";
                cout << "t = " << t << " saved." << endl;
            }
        }
        save_field_csv(t);
        decay_file.close();
    }
};

int main() {
    SolverTG_CSV().run();
    return 0;
}