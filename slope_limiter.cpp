#include "slope_limiter.h"
#include "limiters_math.h"

std::pair<double, double> reconstruct_slope_1d(const std::vector<double>& q, int i, LimiterType type) {
    double dq_L = q[i] - q[i - 1];
    double dq_R = q[i + 1] - q[i];

    double slope = math_limiter(dq_L, dq_R, type);

    return { q[i] - 0.5 * slope, q[i] + 0.5 * slope };
}