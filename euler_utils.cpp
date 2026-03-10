#include "euler_utils.h"
#include <cmath>
#include <algorithm>

Conserved physToCons(const State& W, double gamma) {
    Conserved U;
    U.rho = W.rho;
    U.rhou = W.rho * W.u;
    U.rhov = W.rho * W.v;
    U.E = W.p / (gamma - 1.0) + 0.5 * W.rho * (W.u * W.u + W.v * W.v);
    return U;
}

State consToPhys(const Conserved& U, double gamma) {
    State W;
    W.rho = std::max(1e-9, U.rho);
    W.u = U.rhou / W.rho;
    W.v = U.rhov / W.rho;
    double kinetic = 0.5 * (U.rhou * W.u + U.rhov * W.v);
    double p = (gamma - 1.0) * (U.E - kinetic);
    W.p = std::max(1e-9, p);
    return W;
}

Flux fluxX(const State& W, double gamma) {
    double E_cons = W.p / (gamma - 1.0) + 0.5 * W.rho * (W.u * W.u + W.v * W.v);
    Flux F;
    F.rho_f = W.rho * W.u;
    F.rhou_f = W.rho * W.u * W.u + W.p;
    F.rhov_f = W.rho * W.u * W.v;
    F.E_f = W.u * (E_cons + W.p);
    return F;
}

Flux fluxY(const State& W, double gamma) {
    double E_cons = W.p / (gamma - 1.0) + 0.5 * W.rho * (W.u * W.u + W.v * W.v);
    Flux G;
    G.rho_f = W.rho * W.v;
    G.rhou_f = W.rho * W.v * W.u;
    G.rhov_f = W.rho * W.v * W.v + W.p;
    G.E_f = W.v * (E_cons + W.p);
    return G;
}

double soundSpeed(const State& W, double gamma) {
    if (W.p <= 0 || W.rho <= 0) return 0.0;
    return std::sqrt(gamma * W.p / W.rho);
}