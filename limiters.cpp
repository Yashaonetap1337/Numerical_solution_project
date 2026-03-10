#include "limiters.h"
#include <cmath>
#include <algorithm>


static double sign(double x) {
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
}


static double lim_minmod(double a, double b) {
    if (a * b <= 0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}


static double lim_superbee(double a, double b) {
    if (a * b <= 0) return 0.0;
    double abs_a = std::abs(a);
    double abs_b = std::abs(b);


    double term1 = (2.0 * abs_a < abs_b) ? 2.0 * abs_a : abs_b;
    double term2 = (abs_a < 2.0 * abs_b) ? abs_a : 2.0 * abs_b;

    return sign(a) * std::max(term1, term2);
}


static double lim_van_leer(double a, double b) {
    if (a * b <= 0) return 0.0;
    return (2.0 * a * b) / (a + b);
}


static double lim_kolgan_1972(double a, double b) {
    if (a * b <= 0) return 0.0;
    if (std::abs(a) < std::abs(b)) return a;
    else return b;
}

static double lim_kolgan_1975(double a, double b) {

    double c = 0.5 * (a + b);
    if (a * b <= 0) return 0.0;
    double abs_a = std::abs(a);
    double abs_b = std::abs(b);
    double abs_c = std::abs(c);

    double min_val = std::min({ abs_a, abs_b, abs_c });

    return sign(a) * min_val;
}

static double lim_osher_1984(double a, double b) {
    double ab = a * b;

    if (ab <= 0) return 0.0;
    if (a * a < ab) return a;
    if (b * b < ab) return b;

    return a;
}
double compute_limited_slope(double a, double b, LimiterType type) {
    switch (type) {
    case LimiterType::NONE:
        return 0.0;

    case LimiterType::MINMOD:
        return lim_minmod(a, b);

    case LimiterType::SUPERBEE:
        return lim_superbee(a, b);

    case LimiterType::VAN_LEER:
        return lim_van_leer(a, b);

    case LimiterType::KOLGAN_1972:
        return lim_kolgan_1972(a, b);

    case LimiterType::KOLGAN_1975:
        return lim_kolgan_1975(a, b);

    case LimiterType::OSHER_1984: 
        return lim_osher_1984(a, b);

    default:
        return lim_minmod(a, b); 
    }
}