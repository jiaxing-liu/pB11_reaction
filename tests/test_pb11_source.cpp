#include "pb11_source.h"

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr double kJoulesPerKeV = 1.602176634e-16;
constexpr double kFusionEnergyPerReactionJ =
    8.68e6 * 1.602176634e-19;

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool close_relative(double actual, double expected, double tolerance) {
    if (actual == expected) {
        return true;
    }
    if (!std::isfinite(actual) || !std::isfinite(expected)) {
        return false;
    }
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return scale > 0.0 && std::abs(actual - expected) <= tolerance * scale;
}

bool all_zero(const pb11_instant_source_v1& source) {
    return source.reaction_rate_m3_s == 0.0 &&
           source.proton_source_m3_s == 0.0 &&
           source.boron_source_m3_s == 0.0 &&
           source.helium4_source_m3_s == 0.0 &&
           source.fusion_power_W_m3 == 0.0 &&
           source.electron_power_W_m3 == 0.0 &&
           source.ion_power_W_m3 == 0.0;
}

int call_source(double temperature_keV, double np_m3, double nB_m3,
                double electron_fraction, int method,
                pb11_instant_source_v1* source) {
    return pb11_c_instant_thermal_source(
        temperature_keV * kJoulesPerKeV, np_m3, nB_m3, electron_fraction,
        method, source);
}

void check_failure(double kT_J, double np_m3, double nB_m3,
                   double electron_fraction, int method, int expected_status,
                   const std::string& message) {
    pb11_instant_source_v1 source{
        -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    const int status = pb11_c_instant_thermal_source(
        kT_J, np_m3, nB_m3, electron_fraction, method, &source);
    check(status == expected_status, message + " status");
    check(all_zero(source), message + " clears every output field");
}

}  // namespace

int main() {
    const double np_m3 = 2.0e19;
    const double nB_m3 = 3.0e19;
    const double electron_fraction = 0.30;
    const int method = PB11_REACTIVITY_INTEGRAL;

    double reactivity = 0.0;
    const int reactivity_status =
        pb11_c_reactivity(100.0, method, &reactivity);
    check(reactivity_status == PB11_STATUS_OK,
          "reference reactivity call succeeds");
    check(std::isfinite(reactivity) && reactivity > 0.0,
          "reference reactivity is finite and positive");

    pb11_instant_source_v1 source{};
    int status = call_source(100.0, np_m3, nB_m3, electron_fraction, method,
                             &source);
    check(status == PB11_STATUS_OK,
          "SI instant source succeeds with integral method");
    check(std::isfinite(source.reaction_rate_m3_s) &&
              source.reaction_rate_m3_s > 0.0,
          "reaction source is finite and positive");

    /* Compare independently with the existing nuclear reactivity ABI. */
    const double expected_rate = np_m3 * nB_m3 * reactivity;
    const double expected_power =
        expected_rate * kFusionEnergyPerReactionJ;
    check(close_relative(source.reaction_rate_m3_s, expected_rate, 1.0e-12),
          "joule-to-keV conversion agrees with existing reactivity ABI");
    check(close_relative(source.fusion_power_W_m3, expected_power, 1.0e-12),
          "fusion power uses the specified Q value");
    check(close_relative(source.proton_source_m3_s, -expected_rate, 1.0e-12),
          "proton source consumes one fuel nucleus per reaction");
    check(close_relative(source.boron_source_m3_s, -expected_rate, 1.0e-12),
          "boron source consumes one fuel nucleus per reaction");
    check(close_relative(source.helium4_source_m3_s, 3.0 * expected_rate,
                         1.0e-12),
          "helium source has three-alpha stoichiometry");
    check(close_relative(source.electron_power_W_m3,
                         electron_fraction * expected_power, 1.0e-12),
          "electron power follows the requested fraction");
    check(close_relative(source.ion_power_W_m3,
                         (1.0 - electron_fraction) * expected_power, 1.0e-12),
          "ion power follows the complementary fraction");
    check(close_relative(source.electron_power_W_m3 + source.ion_power_W_m3,
                         source.fusion_power_W_m3, 1.0e-14),
          "electron and ion power sum to fusion power");

    /* Each fuel density scales the reaction and power sources linearly. */
    pb11_instant_source_v1 doubled_protons{};
    pb11_instant_source_v1 doubled_boron{};
    status = call_source(100.0, 2.0 * np_m3, nB_m3, electron_fraction,
                         PB11_REACTIVITY_FAST, &doubled_protons);
    check(status == PB11_STATUS_OK, "doubling proton density succeeds");
    status = call_source(100.0, np_m3, 2.0 * nB_m3, electron_fraction,
                         PB11_REACTIVITY_FAST, &doubled_boron);
    check(status == PB11_STATUS_OK, "doubling boron density succeeds");
    status = call_source(100.0, np_m3, nB_m3, electron_fraction,
                         PB11_REACTIVITY_FAST, &source);
    check(status == PB11_STATUS_OK, "fast source succeeds");
    check(close_relative(doubled_protons.reaction_rate_m3_s,
                         2.0 * source.reaction_rate_m3_s, 1.0e-14),
          "reaction rate scales with proton density");
    check(close_relative(doubled_boron.fusion_power_W_m3,
                         2.0 * source.fusion_power_W_m3, 1.0e-14),
          "fusion power scales with boron density");

    status = call_source(100.0, 0.0, nB_m3, electron_fraction, method,
                         &source);
    check(status == PB11_STATUS_OK && all_zero(source),
          "zero proton density gives an all-zero source");
    status = call_source(100.0, np_m3, 0.0, electron_fraction, method,
                         &source);
    check(status == PB11_STATUS_OK && all_zero(source),
          "zero boron density gives an all-zero source");

    pb11_instant_source_v1 electron_only{};
    pb11_instant_source_v1 ion_only{};
    status = call_source(100.0, np_m3, nB_m3, 0.0, PB11_REACTIVITY_FAST,
                         &ion_only);
    check(status == PB11_STATUS_OK, "zero electron fraction succeeds");
    status = call_source(100.0, np_m3, nB_m3, 1.0, PB11_REACTIVITY_FAST,
                         &electron_only);
    check(status == PB11_STATUS_OK, "unit electron fraction succeeds");
    check(ion_only.electron_power_W_m3 == 0.0 &&
              close_relative(ion_only.ion_power_W_m3,
                             ion_only.fusion_power_W_m3, 1.0e-14),
          "zero electron fraction deposits all power to ions");
    check(electron_only.ion_power_W_m3 == 0.0 &&
              close_relative(electron_only.electron_power_W_m3,
                             electron_only.fusion_power_W_m3, 1.0e-14),
          "unit electron fraction deposits all power to electrons");

    /* Decimal SI endpoint literals can be one ULP above 10 or 500 keV. */
    status = pb11_c_instant_thermal_source(
        1.602176634e-15, np_m3, nB_m3, electron_fraction,
        PB11_REACTIVITY_FAST, &source);
    check(status == PB11_STATUS_OK &&
              std::isfinite(source.reaction_rate_m3_s),
          "10 keV SI endpoint is accepted by the fast method");
    status = pb11_c_instant_thermal_source(
        8.01088317e-14, np_m3, nB_m3, electron_fraction,
        PB11_REACTIVITY_FAST, &source);
    check(status == PB11_STATUS_OK &&
              std::isfinite(source.reaction_rate_m3_s),
          "500 keV SI endpoint is accepted by the fast method");

    check_failure(std::numeric_limits<double>::quiet_NaN(), np_m3, nB_m3,
                  electron_fraction, method, PB11_STATUS_INVALID_ARGUMENT,
                  "NaN thermal energy");
    check_failure(100.0 * kJoulesPerKeV, std::numeric_limits<double>::infinity(),
                  nB_m3, electron_fraction, method, PB11_STATUS_INVALID_ARGUMENT,
                  "infinite proton density");
    check_failure(100.0 * kJoulesPerKeV, np_m3, nB_m3,
                  std::numeric_limits<double>::quiet_NaN(), method,
                  PB11_STATUS_INVALID_ARGUMENT, "NaN electron fraction");
    check_failure(0.0, np_m3, nB_m3, electron_fraction, method,
                  PB11_STATUS_OUT_OF_RANGE, "zero thermal energy");
    check_failure(-kJoulesPerKeV, np_m3, nB_m3, electron_fraction, method,
                  PB11_STATUS_OUT_OF_RANGE, "negative thermal energy");
    check_failure(100.0 * kJoulesPerKeV, -1.0, nB_m3, electron_fraction,
                  method, PB11_STATUS_OUT_OF_RANGE, "negative proton density");
    check_failure(100.0 * kJoulesPerKeV, np_m3, -1.0, electron_fraction,
                  method, PB11_STATUS_OUT_OF_RANGE, "negative boron density");
    check_failure(100.0 * kJoulesPerKeV, np_m3, nB_m3, -0.01, method,
                  PB11_STATUS_OUT_OF_RANGE, "negative electron fraction");
    check_failure(100.0 * kJoulesPerKeV, np_m3, nB_m3, 1.01, method,
                  PB11_STATUS_OUT_OF_RANGE, "electron fraction above one");
    check_failure(100.0 * kJoulesPerKeV, np_m3, nB_m3, electron_fraction, 99,
                  PB11_STATUS_UNKNOWN_METHOD, "unknown reactivity method");
    check_failure(std::numeric_limits<double>::max(), np_m3, nB_m3,
                  electron_fraction, method, PB11_STATUS_NUMERICAL_FAILURE,
                  "joule-to-keV conversion overflow");
    check_failure(100.0 * kJoulesPerKeV,
                  std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), electron_fraction,
                  PB11_REACTIVITY_FAST, PB11_STATUS_NUMERICAL_FAILURE,
                  "unrepresentable reaction-rate output");

    check(pb11_c_instant_thermal_source(100.0 * kJoulesPerKeV, np_m3, nB_m3,
                                        electron_fraction, method, nullptr) ==
              PB11_STATUS_NULL_OUTPUT,
          "null output pointer is reported");

    /* A density pair whose product exceeds double remains valid when R fits. */
    double low_temperature_reactivity = 0.0;
    status = pb11_c_reactivity(0.01, PB11_REACTIVITY_INTEGRAL,
                                &low_temperature_reactivity);
    check(status == PB11_STATUS_OK && low_temperature_reactivity > 0.0 &&
              std::isfinite(low_temperature_reactivity),
          "low-temperature reference reactivity is finite and positive");
    status = pb11_c_instant_thermal_source(
        0.01 * kJoulesPerKeV, 1.0e200, 1.0e200, 0.25,
        PB11_REACTIVITY_INTEGRAL, &source);
    check(status == PB11_STATUS_OK,
          "large intermediate density product with finite final outputs");
    const long double large_rate =
        static_cast<long double>(low_temperature_reactivity) * 1.0e200L *
        1.0e200L;
    const long double large_power =
        large_rate * static_cast<long double>(kFusionEnergyPerReactionJ);
    check(close_relative(source.reaction_rate_m3_s,
                         static_cast<double>(large_rate), 1.0e-12),
          "large density product preserves reaction rate");
    check(close_relative(source.fusion_power_W_m3,
                         static_cast<double>(large_power), 1.0e-12),
          "large density product preserves fusion power");

    /* Nonzero values below the double range are rejected instead of becoming 0. */
    check_failure(100.0 * kJoulesPerKeV,
                  std::numeric_limits<double>::denorm_min(),
                  std::numeric_limits<double>::denorm_min(), 0.5,
                  PB11_REACTIVITY_FAST, PB11_STATUS_NUMERICAL_FAILURE,
                  "nonrepresentable subnormal reaction rate");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All p-11B instant-source tests passed\n";
    return 0;
}
