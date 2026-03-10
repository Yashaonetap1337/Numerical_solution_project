#include "acoustic_riemann.h"
#include "euler_utils.h"
#include <algorithm>

State solve_acoustic_riemann_problem(const State& W_L, const State& W_R, double gamma) {
    double aL = soundSpeed(W_L, gamma);
    double aR = soundSpeed(W_R, gamma);
    double rho_avg = 0.5 * (W_L.rho + W_R.rho);
    double a_avg = 0.5 * (aL + aR);
    double p_star = 0.5 * (W_L.p + W_R.p) - 0.5 * (W_R.u - W_L.u) * rho_avg * a_avg;
    double u_star = 0.5 * (W_L.u + W_R.u) - 0.5 * (W_R.p - W_L.p) / (rho_avg * a_avg);

    State S;
    double SL = W_L.u - aL;
    double SR = W_R.u + aR;
    if (SL > 0.0) {
        S = W_L;
    } else if (SR < 0.0) {
        S = W_R;
    } else {
        S.rho = rho_avg;
        S.u = u_star;
        S.p = std::max(1e-9, p_star);
        S.v = (u_star >= 0.0) ? W_L.v : W_R.v;
    }
    return S;
}