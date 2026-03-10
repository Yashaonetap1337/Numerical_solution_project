#include "limiters_math.h"
#include <cmath>
#include <algorithm>
#include <vector>

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
    double term1 = std::min(2.0 * abs_a, abs_b);
    double term2 = std::min(abs_a, 2.0 * abs_b);
    return sign(a) * std::max(term1, term2);
}

static double lim_van_leer(double a, double b) {
    if (a * b <= 0) return 0.0;
    return (2.0 * a * b) / (a + b);
}

static double lim_kolgan_75(double a, double b) {
    double c = 0.5 * (a + b);
    if (a * b <= 0) return 0.0;
    return sign(a) * std::min({ std::abs(a), std::abs(b), std::abs(c) });
}

static double lim_osher_84(double a, double b) {
    if (a * b <= 0) return 0.0;
    if (a * a < a * b) return a;
    if (b * b < a * b) return b;
    return a;
}

double math_limiter(double a, double b, LimiterType type) {
    switch (type) {
    case LimiterType::MINMOD:      return lim_minmod(a, b);
    case LimiterType::SUPERBEE:    return lim_superbee(a, b);
    case LimiterType::VAN_LEER:    return lim_van_leer(a, b);
    case LimiterType::KOLGAN_1972: return lim_minmod(a, b);
    case LimiterType::KOLGAN_1975: return lim_kolgan_75(a, b);
    case LimiterType::OSHER_1984:  return lim_osher_84(a, b);
    case LimiterType::NONE:        return 0.0;
    default: return lim_minmod(a, b);
    }
}