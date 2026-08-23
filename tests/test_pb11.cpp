#include "pb11.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

double relative_difference(double lhs, double rhs) {
    return std::abs(lhs - rhs) / std::abs(rhs);
}

void check_close(double actual, double expected, double relative_tolerance,
                 const std::string& message) {
    const double scale = std::max(std::abs(actual), std::abs(expected));
    check(std::abs(actual - expected) <= relative_tolerance * scale, message);
}

}  // namespace

int main() {
    check_close(pb11_sfactor(0.148), 3537.9844888, 1.0e-8,
                "148 keV S-factor regression value");
    check(pb11_sfactor(0.148) > 5.0 * pb11_sfactor(0.140),
          "148 keV narrow resonance is resolved");

    const double s_left_0400 = pb11_sfactor(std::nextafter(0.400, 0.0));
    const double s_right_0400 = pb11_sfactor(std::nextafter(0.400, 1.0));
    const double s_left_0668 = pb11_sfactor(std::nextafter(0.668, 0.0));
    const double s_right_0668 = pb11_sfactor(std::nextafter(0.668, 1.0));
    check(relative_difference(s_left_0400, s_right_0400) < 2.0e-3,
          "S-factor joins at 0.400 MeV within rounded-parameter precision");
    check(relative_difference(s_left_0668, s_right_0668) < 1.0e-2,
          "S-factor joins at 0.668 MeV within rounded-parameter precision");

    check(pb11_cross_section(0.0) == 0.0, "cross section has zero-energy limit");
    check(std::isfinite(pb11_cross_section(std::numeric_limits<double>::denorm_min())),
          "cross section is finite at a subnormal positive energy");
    check(pb11_cross_section(0.6) > 1.0 && pb11_cross_section(0.6) < 1.5,
          "main cross-section resonance has Figure 1 magnitude");
    check(pb11_cross_section(2.3) > 0.4 && pb11_cross_section(2.3) < 0.8,
          "second broad cross-section resonance has Figure 1 magnitude");

    check(std::isnan(pb11_sfactor(-1.0)), "negative S-factor energy is rejected");
    check(std::isnan(pb11_sfactor(9.760001)), "energy above fit range is rejected");
    check(std::isnan(pb11_cross_section(-1.0)),
          "negative cross-section energy is rejected");
    check(std::isnan(pb11_reactivity_integral(0.0)),
          "non-positive integral temperature is rejected");
    check(std::isnan(pb11_reactivity_fast(9.999)),
          "fast approximation below published range is rejected");
    check(std::isnan(pb11_reactivity_fast(500.001)),
          "fast approximation above published range is rejected");

    constexpr std::array<double, 14> temperatures{
        10.0, 20.0, 30.0, 50.0, 60.0, 69.0, 69.9,
        70.0, 70.1, 100.0, 200.0, 300.0, 400.0, 500.0};
    for (double temperature : temperatures) {
        const double integral = pb11_reactivity_integral(temperature);
        const double fast = pb11_reactivity_fast(temperature);
        check(std::isfinite(integral) && integral > 0.0,
              "integral reactivity is finite and positive");
        check(std::isfinite(fast) && fast > 0.0,
              "fast reactivity is finite and positive");
        const double paper_bound = temperature < 70.0 ? 0.0205 : 0.0105;
        check(relative_difference(fast, integral) < paper_bound,
              "analytic reactivity satisfies the paper error bound at T=" +
                  std::to_string(temperature) + " keV");
    }

    const double fast_below = pb11_reactivity_fast(69.999);
    const double fast_at = pb11_reactivity_fast(70.0);
    check(relative_difference(fast_below, fast_at) < 0.025,
          "LT and HT fits have no large jump at 70 keV");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All p-11B tests passed\n";
    return 0;
}
