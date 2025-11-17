#include "euler_utils.h"

#include <cmath>
#include <iostream>

Conserved physToCons(const State& W, double gamma) {
    Conserved U;
    U.rho = W.rho;
    U.rhou = W.rho * W.u;
    U.E = W.p / (gamma - 1.0) + 0.5 * W.rho * W.u * W.u;
    return U;
}

State consToPhys(const Conserved& U, double gamma) {
    State W;
    W.rho = U.rho;
    W.u = U.rhou / U.rho;

    double pressure = (gamma - 1.0) * (U.E - 0.5 * W.u * U.rhou);

    W.p = std::max(1e-9, pressure);
    return W;
}

Flux physToFlux(const State& W, double gamma) {
    Flux F;
    F.rho_f = W.rho * W.u;
    F.rhou_f = W.rho * W.u * W.u + W.p;
    F.E_f = (W.p / (gamma - 1.0) + 0.5 * W.rho * W.u * W.u + W.p) * W.u;
    return F;
}

double soundSpeed(const State& W, double gamma) {
    // скорость звука не может быть вычислена для отрицательного давления/плотности
    if (W.p <= 0 || W.rho <= 0) return 0.0;
    return std::sqrt(gamma * W.p / W.rho);
}