#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

// Параметры сетки и физики
const int NX = 31;
const int NY = 31;
const double LX = 1.0, LY = 1.0;
const double RHO = 1.0, RE = 100.0, U_LID = 1.0;

// Параметры PIMPLE
const double DT = 0.02;           // Шаг по времени (может быть больше, чем для PISO!)
const double T_END = 10.0;        // Общее время расчета
const int N_OUTER = 3;            // Количество внешних итераций (SIMPLE)
const int N_CORR = 1;             // Количество коррекций давления (PISO)
const double ALPHA_U = 0.7;       // Релаксация скорости (для промежуточных итераций)
const double ALPHA_P = 0.3;       // Релаксация давления (для промежуточных итераций)

struct PIMPLE_Solver {
    double dx, dy, nu;
    vector<double> u, v, p;           // Текущие поля (результат шага)
    vector<double> u_n, v_n;          // Поля с начала временного шага
    vector<double> u_work, v_work, p_work; // Рабочие массивы для внешних итераций
    vector<double> auP, avP;          // Эффективные диагонали

    PIMPLE_Solver() :
        dx(LX / NX), dy(LY / NY), nu(U_LID* LX / RE),
        u((NX + 1)* NY, 0), v(NX* (NY + 1), 0), p(NX* NY, 0),
        u_n((NX + 1)* NY, 0), v_n(NX* (NY + 1), 0),
        u_work((NX + 1)* NY, 0), v_work(NX* (NY + 1), 0), p_work(NX* NY, 0),
        auP((NX + 1)* NY, 0), avP(NX* (NY + 1), 0) {}

    // --- Предиктор импульса для U ---
    void momentumPredictorU(double alpha_u) {
        for (int i = 1; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                // Интерполяция скоростей на грани контрольного объема
                double uc = u_work[i * NY + j];
                double uE = u_work[(i + 1) * NY + j];
                double uW = u_work[(i - 1) * NY + j];

                double ue = 0.5 * (uc + uE);
                double uw = 0.5 * (uc + uW);

                double vn = 0.5 * (v_work[(i - 1) * (NY + 1) + j + 1] + v_work[i * (NY + 1) + j + 1]);
                double vs = 0.5 * (v_work[(i - 1) * (NY + 1) + j] + v_work[i * (NY + 1) + j]);

                // Массовые потоки
                double Fe = RHO * ue * dy;
                double Fw = RHO * uw * dy;
                double Fn = RHO * vn * dx;
                double Fs = RHO * vs * dx;

                // Коэффициенты (upwind + диффузия)
                double aE = max(-Fe, 0.0) + nu * dy / dx;
                double aW = max(Fw, 0.0) + nu * dy / dx;

                double aN, aS, uN_ghost, uS_ghost;

                // Граничные условия
                if (j == NY - 1) { // Верхняя крышка (движущаяся стенка)
                    aN = max(-Fn, 0.0) + nu * dx / (0.5 * dy);
                    uN_ghost = 2 * U_LID - uc;
                }
                else {
                    aN = max(-Fn, 0.0) + nu * dx / dy;
                    uN_ghost = u_work[i * NY + j + 1];
                }

                if (j == 0) { // Нижняя стенка
                    aS = max(Fs, 0.0) + nu * dx / (0.5 * dy);
                    uS_ghost = -uc;
                }
                else {
                    aS = max(Fs, 0.0) + nu * dx / dy;
                    uS_ghost = u_work[i * NY + j - 1];
                }

                // Диагональ с временным членом
                double aP_base = aE + aW + aN + aS + RHO * dx * dy / DT;

                // Оператор H (вклады соседей + временной член)
                double H = aE * uE + aW * uW + aN * uN_ghost + aS * uS_ghost
                    + (RHO * dx * dy / DT) * u_n[i * NY + j];

                // Градиент давления
                double dP = (p_work[i * NY + j] - p_work[(i - 1) * NY + j]) * dy;

                // Предикторная скорость с под-релаксацией
                u_work[i * NY + j] = (1.0 - alpha_u) * u_work[i * NY + j]
                    + (alpha_u / aP_base) * (H - dP);

                // Эффективная диагональ
                auP[i * NY + j] = aP_base / alpha_u;
            }
        }

        // Граничные условия для u
        for (int j = 0; j < NY; j++) {
            u_work[0 * NY + j] = 0;      // Левая стенка
            u_work[NX * NY + j] = 0;     // Правая стенка
        }
    }

    // --- Предиктор импульса для V ---
    void momentumPredictorV(double alpha_u) {
        for (int i = 0; i < NX; i++) {
            for (int j = 1; j < NY; j++) {
                // Интерполяция скоростей на грани контрольного объема
                double vc = v_work[i * (NY + 1) + j];
                double vN = v_work[i * (NY + 1) + j + 1];
                double vS = v_work[i * (NY + 1) + j - 1];

                double vn = 0.5 * (vc + vN);
                double vsf = 0.5 * (vc + vS);

                double ue = 0.5 * (u_work[(i + 1) * NY + j - 1] + u_work[(i + 1) * NY + j]);
                double uw = 0.5 * (u_work[i * NY + j - 1] + u_work[i * NY + j]);

                // Массовые потоки
                double Fe = RHO * ue * dy;
                double Fw = RHO * uw * dy;
                double Fn = RHO * vn * dx;
                double Fs = RHO * vsf * dx;

                // Коэффициенты
                double aN = max(-Fn, 0.0) + nu * dx / dy;
                double aS = max(Fs, 0.0) + nu * dx / dy;

                double aE, aW, vE_ghost, vW_ghost;

                // Граничные условия
                if (i == NX - 1) { // Правая стенка
                    aE = max(-Fe, 0.0) + nu * dy / (0.5 * dx);
                    vE_ghost = -vc;
                }
                else {
                    aE = max(-Fe, 0.0) + nu * dy / dx;
                    vE_ghost = v_work[(i + 1) * (NY + 1) + j];
                }

                if (i == 0) { // Левая стенка
                    aW = max(Fw, 0.0) + nu * dy / (0.5 * dx);
                    vW_ghost = -vc;
                }
                else {
                    aW = max(Fw, 0.0) + nu * dy / dx;
                    vW_ghost = v_work[(i - 1) * (NY + 1) + j];
                }

                // Диагональ с временным членом
                double aP_base = aE + aW + aN + aS + RHO * dx * dy / DT;

                // Оператор H
                double H = aE * vE_ghost + aW * vW_ghost + aN * vN + aS * vS
                    + (RHO * dx * dy / DT) * v_n[i * (NY + 1) + j];

                // Градиент давления
                double dP = (p_work[i * NY + j] - p_work[i * NY + j - 1]) * dx;

                // Предикторная скорость с под-релаксацией
                v_work[i * (NY + 1) + j] = (1.0 - alpha_u) * v_work[i * (NY + 1) + j]
                    + (alpha_u / aP_base) * (H - dP);

                // Эффективная диагональ
                avP[i * (NY + 1) + j] = aP_base / alpha_u;
            }
        }

        // Граничные условия для v
        for (int i = 0; i < NX; i++) {
            v_work[i * (NY + 1) + 0] = 0;      // Нижняя стенка
            v_work[i * (NY + 1) + NY] = 0;     // Верхняя стенка
        }
    }

    // --- Решатель уравнения Пуассона ---
    void solvePressure(double alpha_p) {
        vector<double> p_corr(NX * NY, 0);

        // Итерации Гаусса-Зейделя
        for (int step = 0; step < 100; step++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // Фиксация давления в точке (0,0)
                    if (i == 0 && j == 0) {
                        p_corr[i * NY + j] = 0;
                        continue;
                    }

                    // Коэффициенты уравнения Пуассона
                    double aE = (i < NX - 1) ? (dy * dy / auP[(i + 1) * NY + j]) : 0;
                    double aW = (i > 0) ? (dy * dy / auP[i * NY + j]) : 0;
                    double aN = (j < NY - 1) ? (dx * dx / avP[i * (NY + 1) + j + 1]) : 0;
                    double aS = (j > 0) ? (dx * dx / avP[i * (NY + 1) + j]) : 0;
                    double aP = aE + aW + aN + aS;

                    // Невязка массы (правая часть со знаком минус!)
                    double b = -((u_work[(i + 1) * NY + j] - u_work[i * NY + j]) * dy
                        + (v_work[i * (NY + 1) + j + 1] - v_work[i * (NY + 1) + j]) * dx);

                    double pE = (i < NX - 1) ? p_corr[(i + 1) * NY + j] : 0;
                    double pW = (i > 0) ? p_corr[(i - 1) * NY + j] : 0;
                    double pN = (j < NY - 1) ? p_corr[i * NY + j + 1] : 0;
                    double pS = (j > 0) ? p_corr[i * NY + j - 1] : 0;

                    if (aP > 1e-12) {
                        p_corr[i * NY + j] = (aE * pE + aW * pW + aN * pN + aS * pS + b) / aP;
                    }
                }
            }
        }

        // Коррекция давления с релаксацией
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                p_work[i * NY + j] += alpha_p * p_corr[i * NY + j];
            }
        }

        // Коррекция скоростей (без релаксации!)
        for (int i = 1; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                u_work[i * NY + j] -= (dy / auP[i * NY + j]) * (p_corr[i * NY + j] - p_corr[(i - 1) * NY + j]);
            }
        }

        for (int i = 0; i < NX; i++) {
            for (int j = 1; j < NY; j++) {
                v_work[i * (NY + 1) + j] -= (dx / avP[i * (NY + 1) + j]) * (p_corr[i * NY + j] - p_corr[i * NY + j - 1]);
            }
        }
    }
    //CHECK: PIMPLE_loop
    void run() {
        double t = 0;
        int time_step = 0;

        cout << "Starting PIMPLE Solver for Lid-Driven Cavity (Re=" << RE << ")..." << endl;
        cout << "DT=" << DT << ", T_END=" << T_END << ", N_OUTER=" << N_OUTER
            << ", N_CORR=" << N_CORR << endl;

        // ЦИКЛ ПО ВРЕМЕНИ
        while (t < T_END) {
            // Сохраняем поля с начала временного шага
            u_n = u;
            v_n = v;

            // Инициализируем рабочие массивы
            u_work = u;
            v_work = v;
            p_work = p;

            // ВНЕШНИЙ ЦИКЛ (SIMPLE итерации)
            for (int outer = 0; outer < N_OUTER; outer++) {
                // Выбор режима релаксации
                double alpha_u_current, alpha_p_current;

                if (outer < N_OUTER - 1) {
                    // Промежуточные итерации - режим SIMPLE (с релаксацией)
                    alpha_u_current = ALPHA_U;
                    alpha_p_current = ALPHA_P;
                }
                else {
                    // Последняя итерация - режим PISO (без релаксации)
                    alpha_u_current = 1.0;
                    alpha_p_current = 1.0;
                }

                // Предиктор импульса
                momentumPredictorU(alpha_u_current);
                momentumPredictorV(alpha_u_current);

                // ВНУТРЕННИЙ ЦИКЛ (PISO коррекции)
                for (int corr = 0; corr < N_CORR; corr++) {
                    solvePressure(alpha_p_current);
                }
            }

            // Обновляем результаты временного шага
            u = u_work;
            v = v_work;
            p = p_work;

            t += DT;
            time_step++;

            if (time_step % 20 == 0) {
                cout << "Time step: " << time_step << ", t = " << fixed << setprecision(3) << t << " s" << endl;

                // Проверка невязки массы
                double max_div = 0;
                for (int i = 0; i < NX; i++) {
                    for (int j = 0; j < NY; j++) {
                        double div = (u[(i + 1) * NY + j] - u[i * NY + j]) * dy
                            + (v[i * (NY + 1) + j + 1] - v[i * (NY + 1) + j]) * dx;
                        max_div = max(max_div, abs(div));
                    }
                }
                cout << "  Max divergence: " << scientific << max_div << endl;
            }
        }

        cout << "\nSimulation completed!" << endl;
        saveToCSV();
    }

    void saveToCSV() {
        ofstream out("cavity_pimple.csv");
        out << "x,y,u,v,p\n";
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double x = (i + 0.5) * dx;
                double y = (j + 0.5) * dy;

                // Интерполяция скоростей в центр ячейки
                double uc = 0.5 * (u[i * NY + j] + u[(i + 1) * NY + j]);
                double vc = 0.5 * (v[i * (NY + 1) + j] + v[i * (NY + 1) + j + 1]);
                double pc = p[i * NY + j];

                out << x << "," << y << "," << uc << "," << vc << "," << pc << "\n";
            }
        }
        out.close();
        cout << "Results saved to cavity_pimple.csv" << endl;

        // Вывод профиля U на центральной линии
        cout << "\nU-velocity profile at x = 0.5 (for verification with Ghia et al.):" << endl;
        int mid_x = NX / 2;
        cout << "y\tu\n";
        for (int j = NY - 1; j >= 0; j--) {
            double y = (j + 0.5) * dy;
            cout << fixed << setprecision(4) << y << "\t" << u[mid_x * NY + j] << endl;
        }
    }
};

int main() {
    PIMPLE_Solver solver;
    solver.run();
    return 0;
}