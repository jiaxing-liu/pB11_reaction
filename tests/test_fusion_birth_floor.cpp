#include "fusion_birth_floor.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr std::size_t kSpecies = 6;
constexpr double kKiloelectronVoltJ = 1.602176634e-16;
constexpr double kD102IonTemperatureKeV = 0.10195170940710378;
constexpr double kD102FirstCenterKeV = 0.5e-10;

int failures = 0;
int checks = 0;

void check(bool condition, const std::string& message) {
    ++checks;
    if (!condition) {
        ++failures;
        std::cerr << "ERROR: " << message << '\n';
    }
}

bool same_bits(double left, double right) {
    return std::memcmp(&left, &right, sizeof(double)) == 0;
}

bool ledger_is_zero(const fusion_birth_floor_ledger_v1& ledger) {
    const double zero = 0.0;
    for (std::size_t i = 0; i < kSpecies; ++i) {
        if (!same_bits(ledger.born_number_m3[i], zero) ||
            !same_bits(ledger.born_energy_J_m3[i], zero) ||
            !same_bits(ledger.mapped_number_m3[i], zero) ||
            !same_bits(ledger.mapped_energy_J_m3[i], zero) ||
            !same_bits(ledger.ion_energy_correction_J_m3[i], zero) ||
            !same_bits(ledger.energy_residual_J_m3[i], zero)) {
            return false;
        }
    }
    return same_bits(ledger.remaining_ion_energy_J_m3, zero);
}

fusion_birth_floor_ledger_v1 sentinel_ledger() {
    fusion_birth_floor_ledger_v1 ledger{};
    for (std::size_t i = 0; i < kSpecies; ++i) {
        ledger.born_number_m3[i] = 11.0 + static_cast<double>(i);
        ledger.born_energy_J_m3[i] = 21.0 + static_cast<double>(i);
        ledger.mapped_number_m3[i] = 31.0 + static_cast<double>(i);
        ledger.mapped_energy_J_m3[i] = 41.0 + static_cast<double>(i);
        ledger.ion_energy_correction_J_m3[i] = 51.0 + static_cast<double>(i);
        ledger.energy_residual_J_m3[i] = 61.0 + static_cast<double>(i);
    }
    ledger.remaining_ion_energy_J_m3 = 71.0;
    return ledger;
}

fusion_birth_floor_options_v1 options(double first_center_J,
                                      double ion_kT_J,
                                      double ion_energy_J_m3,
                                      double max_center_over_ion_kT,
                                      double max_ion_energy_fraction) {
    fusion_birth_floor_options_v1 result{};
    result.first_center_J = first_center_J;
    result.ion_kT_J = ion_kT_J;
    result.ion_energy_J_m3 = ion_energy_J_m3;
    result.max_center_over_ion_kT = max_center_over_ion_kT;
    result.max_ion_energy_fraction = max_ion_energy_fraction;
    return result;
}

fusion_birth_floor_options_v1 ordinary_options() {
    return options(0.5, 1.0, 1000.0, 1.0, 0.5);
}

void expect_status(int actual, int expected, const std::string& case_name) {
    check(actual == expected,
          case_name + ": status=" + std::to_string(actual) +
              ", expected=" + std::to_string(expected));
}

bool close_enough(double actual, long double expected) {
    const long double scale = std::abs(expected) + std::abs(static_cast<long double>(actual));
    const long double tolerance =
        8.0L * static_cast<long double>(std::numeric_limits<double>::epsilon()) * scale +
        4.0L * static_cast<long double>(std::numeric_limits<double>::denorm_min());
    return std::isfinite(actual) && std::abs(static_cast<long double>(actual) - expected) <= tolerance;
}

void test_nominal_projection_and_identities() {
    const fusion_birth_floor_options_v1 o = ordinary_options();
    const double born_number[kSpecies] = {1.0, 2.5, 0.0, 3.0, 1000.0, 4.0};
    const double born_energy[kSpecies] = {0.25, 1.25, 0.0, 0.0, 375.0, 0.2};
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();

    expect_status(fusion_c_birth_floor_project(&o, born_number, born_energy, &out),
                  PB11_STATUS_OK, "nominal projection");

    long double returned_borrowed = 0.0L;
    for (std::size_t i = 0; i < kSpecies; ++i) {
        check(same_bits(out.born_number_m3[i], born_number[i]),
              "nominal species " + std::to_string(i) + " preserves born N exactly");
        check(same_bits(out.mapped_number_m3[i], born_number[i]),
              "nominal species " + std::to_string(i) + " preserves mapped N exactly");
        check(same_bits(out.born_energy_J_m3[i], born_energy[i]),
              "nominal species " + std::to_string(i) + " preserves born U exactly");

        const long double expected_mapped =
            static_cast<long double>(born_number[i]) * static_cast<long double>(o.first_center_J);
        check(close_enough(out.mapped_energy_J_m3[i], expected_mapped),
              "nominal species " + std::to_string(i) + " mapped energy");

        /* The identity is formed from the returned mapped double, as required by the API. */
        const long double expected_correction =
            static_cast<long double>(born_energy[i]) -
            static_cast<long double>(out.mapped_energy_J_m3[i]);
        check(close_enough(out.ion_energy_correction_J_m3[i], expected_correction),
              "nominal species " + std::to_string(i) + " signed correction U-Umap");
        const long double residual =
            static_cast<long double>(out.mapped_energy_J_m3[i]) +
            static_cast<long double>(out.ion_energy_correction_J_m3[i]) -
            static_cast<long double>(born_energy[i]);
        check(std::abs(residual) <=
                  8.0L * static_cast<long double>(std::numeric_limits<double>::epsilon()) *
                      (std::abs(static_cast<long double>(born_energy[i])) + 1.0L) +
                  4.0L * static_cast<long double>(std::numeric_limits<double>::denorm_min()),
              "nominal species " + std::to_string(i) + " energy identity residual");
        returned_borrowed -= static_cast<long double>(out.ion_energy_correction_J_m3[i]);
    }

    check(close_enough(out.remaining_ion_energy_J_m3,
                       static_cast<long double>(o.ion_energy_J_m3) - returned_borrowed),
          "nominal remaining reservoir follows returned correction");
}

void test_d102_ratio_fixture() {
    /* These are the D102 numerical values, in joules.  N and U remain per-trial amounts. */
    const fusion_birth_floor_options_v1 o = options(
        kD102FirstCenterKeV * kKiloelectronVoltJ,
        kD102IonTemperatureKeV * kKiloelectronVoltJ,
        1.0e-9,
        0.5,
        1.0);
    const double born_number[kSpecies] = {1.0e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    double born_energy[kSpecies] = {};
    born_energy[0] = 0.25 * born_number[0] * o.first_center_J;
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();

    expect_status(fusion_c_birth_floor_project(&o, born_number, born_energy, &out),
                  PB11_STATUS_OK, "D102 ratio fixture");
    check(o.first_center_J / o.ion_kT_J < o.max_center_over_ion_kT,
          "D102 center/kT ratio is below the configured gate");
    check(same_bits(out.born_number_m3[0], born_number[0]),
          "D102 amount preserves per-trial N exactly");
    check(close_enough(out.ion_energy_correction_J_m3[0],
                       static_cast<long double>(born_energy[0]) -
                           static_cast<long double>(out.mapped_energy_J_m3[0])),
          "D102 amount uses returned mapped energy in signed correction");
}

void test_aggregate_six_species_reservoir_gate() {
    const double born_number[kSpecies] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const double born_energy[kSpecies] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    const fusion_birth_floor_options_v1 accepted = options(1.0, 2.0, 10.0, 0.5, 0.4);
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&accepted, born_number, born_energy, &out),
                  PB11_STATUS_OK, "six-species aggregate reservoir under limit");
    check(same_bits(out.remaining_ion_energy_J_m3, 7.0),
          "six-species aggregate borrowed energy leaves 7 J/m3");

    const fusion_birth_floor_options_v1 rejected = options(1.0, 2.0, 10.0, 0.5, 0.25);
    out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&rejected, born_number, born_energy, &out),
                  PB11_STATUS_OUT_OF_RANGE, "six-species aggregate reservoir over limit");
    check(ledger_is_zero(out), "aggregate reservoir rejection is atomic");
}

void test_center_temperature_gate() {
    const double born_number[kSpecies] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double born_energy[kSpecies] = {};

    const fusion_birth_floor_options_v1 equality = options(1.0, 2.0, 10.0, 0.5, 1.0);
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&equality, born_number, born_energy, &out),
                  PB11_STATUS_OK, "center/kT equality boundary");

    const fusion_birth_floor_options_v1 rejected = options(1.0, 2.0, 10.0, 0.49, 1.0);
    out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&rejected, born_number, born_energy, &out),
                  PB11_STATUS_OUT_OF_RANGE, "center/kT gate rejection");
    check(ledger_is_zero(out), "center/kT rejection is atomic");
}

void test_empty_input() {
    const fusion_birth_floor_options_v1 o = options(1.0, 2.0, 42.0, 0.5, 0.0);
    const double born_number[kSpecies] = {};
    const double born_energy[kSpecies] = {};
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();

    expect_status(fusion_c_birth_floor_project(&o, born_number, born_energy, &out),
                  PB11_STATUS_OK, "empty input");
    for (std::size_t i = 0; i < kSpecies; ++i) {
        check(same_bits(out.born_number_m3[i], 0.0) &&
                  same_bits(out.mapped_number_m3[i], 0.0) &&
                  same_bits(out.born_energy_J_m3[i], 0.0) &&
                  same_bits(out.mapped_energy_J_m3[i], 0.0) &&
                  same_bits(out.ion_energy_correction_J_m3[i], 0.0) &&
                  same_bits(out.energy_residual_J_m3[i], 0.0),
              "empty input species ledger is zero");
    }
    check(same_bits(out.remaining_ion_energy_J_m3, o.ion_energy_J_m3),
          "empty input leaves ion reservoir unchanged");
}

void test_depleted_reservoir_boundaries() {
    const double born_number[kSpecies] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double born_energy[kSpecies] = {};
    /* This is the exact borrowed == ion-energy boundary. */
    const fusion_birth_floor_options_v1 exhausted = options(2.0, 2.0, 2.0, 1.0, 1.0);
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&exhausted, born_number, born_energy, &out),
                  PB11_STATUS_OUT_OF_RANGE, "exactly exhausted reservoir");
    check(ledger_is_zero(out), "exact reservoir exhaustion is atomic");

    const fusion_birth_floor_options_v1 positive = options(2.0, 2.0, 2.5, 1.0, 1.0);
    out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&positive, born_number, born_energy, &out),
                  PB11_STATUS_OK, "positive residual reservoir");
    check(same_bits(out.remaining_ion_energy_J_m3, 0.5),
          "positive reservoir boundary retains positive residual");
}

void expect_invalid_zero(const std::string& case_name,
                         fusion_birth_floor_options_v1* o,
                         const double* born_number,
                         const double* born_energy,
                         int expected_status = PB11_STATUS_INVALID_ARGUMENT) {
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(o, born_number, born_energy, &out),
                  expected_status, case_name);
    check(ledger_is_zero(out), case_name + " zeros every output field");
}

void test_invalid_nan_and_overflow_atomicity() {
    const double valid_number[kSpecies] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double valid_energy[kSpecies] = {};

    fusion_birth_floor_options_v1 o = ordinary_options();
    expect_invalid_zero("null options", nullptr, valid_number, valid_energy);
    expect_invalid_zero("null born numbers", &o, nullptr, valid_energy);
    expect_invalid_zero("null born energies", &o, valid_number, nullptr);

    o = ordinary_options();
    o.first_center_J = std::numeric_limits<double>::quiet_NaN();
    expect_invalid_zero("NaN option", &o, valid_number, valid_energy);

    o = ordinary_options();
    double nan_number[kSpecies] = {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0, 0.0, 0.0};
    expect_invalid_zero("NaN born number", &o, nan_number, valid_energy);

    o = ordinary_options();
    double nan_energy[kSpecies] = {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0, 0.0, 0.0};
    expect_invalid_zero("NaN born energy", &o, valid_number, nan_energy);

    o = ordinary_options();
    double negative_number[kSpecies] = {-1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    expect_invalid_zero("negative born number", &o, negative_number, valid_energy);

    o = ordinary_options();
    double negative_energy[kSpecies] = {-1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    expect_invalid_zero("negative born energy", &o, valid_number, negative_energy);

    o = ordinary_options();
    double infinite_number[kSpecies] = {std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0};
    expect_invalid_zero("infinite born number", &o, infinite_number, valid_energy);

    o = options(2.0, 2.0, std::numeric_limits<double>::max(), 1.0, 1.0);
    double huge_number[kSpecies] = {std::numeric_limits<double>::max(), 0.0, 0.0, 0.0, 0.0, 0.0};
    expect_invalid_zero("mapped energy overflow", &o, huge_number, valid_energy,
                        PB11_STATUS_NUMERICAL_FAILURE);

    o = ordinary_options();
    expect_status(fusion_c_birth_floor_project(&o, valid_number, valid_energy, nullptr),
                  PB11_STATUS_NULL_OUTPUT, "null output");
}

void test_subnormal_numbers_and_energy() {
    const double denorm = std::numeric_limits<double>::denorm_min();
    const fusion_birth_floor_options_v1 o = options(1.0, 1.0, 4.0, 1.0, 1.0);
    const double born_number[kSpecies] = {denorm, 1.0, 0.0, 0.0, 0.0, 0.0};
    const double born_energy[kSpecies] = {0.0, denorm, 0.0, 0.0, 0.0, 0.0};
    fusion_birth_floor_ledger_v1 out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&o, born_number, born_energy, &out),
                  PB11_STATUS_OK, "subnormal N/U");
    check(same_bits(out.born_number_m3[0], denorm) &&
              same_bits(out.mapped_number_m3[0], denorm),
          "subnormal N is copied exactly");
    check(same_bits(out.born_energy_J_m3[1], denorm),
          "subnormal U is copied exactly");
    check(close_enough(out.ion_energy_correction_J_m3[1],
                       static_cast<long double>(denorm) -
                           static_cast<long double>(out.mapped_energy_J_m3[1])),
          "subnormal U participates in signed correction");

    /* A positive mapped amount may underflow to returned double zero. */
    const fusion_birth_floor_options_v1 underflow =
        options(denorm, 1.0, 1.0, 1.0, 1.0);
    const double underflow_number[kSpecies] = {0.5, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double underflow_energy[kSpecies] = {};
    out = sentinel_ledger();
    expect_status(fusion_c_birth_floor_project(&underflow, underflow_number,
                                               underflow_energy, &out),
                  PB11_STATUS_OK, "positive mapped amount underflows to zero");
    check(same_bits(out.mapped_energy_J_m3[0], 0.0),
          "underflowed mapped energy is returned as zero");
    check(same_bits(out.ion_energy_correction_J_m3[0], 0.0),
          "underflow correction uses returned mapped zero");
}

}  // namespace

int main() {
    test_nominal_projection_and_identities();
    test_d102_ratio_fixture();
    test_aggregate_six_species_reservoir_gate();
    test_center_temperature_gate();
    test_empty_input();
    test_depleted_reservoir_boundaries();
    test_invalid_nan_and_overflow_atomicity();
    test_subnormal_numbers_and_energy();

    if (failures != 0) {
        std::cerr << "FAIL fusion_birth_floor: " << failures << " of " << checks
                  << " checks failed\n";
        return 1;
    }
    std::cout << "PASS fusion_birth_floor: " << checks << " checks\n";
    return 0;
}
