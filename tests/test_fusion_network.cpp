#include "fusion_network.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr std::array<double, FUSION_SPECIES_COUNT> kMassNumbers{{
    1.0, 2.0, 3.0, 3.0, 4.0, 11.0}};
constexpr std::array<double, FUSION_SPECIES_COUNT> kCharges{{
    1.0, 1.0, 1.0, 2.0, 2.0, 5.0}};

/* Independent expected stoichiometry, in p,D,T,He3,He4,B11 order. */
constexpr std::array<std::array<double, FUSION_SPECIES_COUNT>,
                     FUSION_CHANNEL_COUNT>
    kExpectedLoss{{
        {{1.0, 0.0, 0.0, 0.0, 0.0, 1.0}},
        {{0.0, 2.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.0, 2.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.0, 1.0, 1.0, 0.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0, 1.0, 0.0, 0.0}},
    }};
constexpr std::array<std::array<double, FUSION_SPECIES_COUNT>,
                     FUSION_CHANNEL_COUNT>
    kExpectedBirth{{
        {{0.0, 0.0, 0.0, 0.0, 3.0, 0.0}},
        {{1.0, 0.0, 1.0, 0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}},
        {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0}},
    }};
constexpr std::array<double, FUSION_CHANNEL_COUNT> kExpectedNeutrons{{
    0.0, 0.0, 1.0, 1.0, 0.0}};

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool all_zero(const fusion_particle_sources_v1& source) {
    for (int species = 0; species < FUSION_SPECIES_COUNT; ++species) {
        if (source.reactant_loss[species] != 0.0 ||
            source.product_birth[species] != 0.0 ||
            source.net_source[species] != 0.0) {
            return false;
        }
    }
    return source.neutron_birth == 0.0;
}

void check_closure(const std::array<double, FUSION_SPECIES_COUNT>& loss,
                   const std::array<double, FUSION_SPECIES_COUNT>& birth,
                   double neutron, const std::string& label) {
    double loss_nucleons = 0.0;
    double birth_nucleons = neutron;
    double loss_charge = 0.0;
    double birth_charge = 0.0;
    for (int species = 0; species < FUSION_SPECIES_COUNT; ++species) {
        loss_nucleons += kMassNumbers[species] * loss[species];
        birth_nucleons += kMassNumbers[species] * birth[species];
        loss_charge += kCharges[species] * loss[species];
        birth_charge += kCharges[species] * birth[species];
    }
    check(loss_nucleons == birth_nucleons,
          label + " conserves nucleon number including neutrons");
    check(loss_charge == birth_charge,
          label + " conserves electric charge");
}

void check_channel(int channel) {
    double rates[FUSION_CHANNEL_COUNT]{};
    rates[channel] = 1.0;
    fusion_particle_sources_v1 source{
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        -1.0};

    const int status = fusion_c_particle_sources(rates, &source);
    const std::string label = "unit channel " + std::to_string(channel);
    check(status == PB11_STATUS_OK, label + " returns OK");
    for (int species = 0; species < FUSION_SPECIES_COUNT; ++species) {
        check(source.reactant_loss[species] == kExpectedLoss[channel][species],
              label + " has the expected reactant loss");
        check(source.product_birth[species] == kExpectedBirth[channel][species],
              label + " has the expected product birth");
        check(source.net_source[species] ==
                  kExpectedBirth[channel][species] -
                      kExpectedLoss[channel][species],
              label + " separates birth and loss in net source");
    }
    check(source.neutron_birth == kExpectedNeutrons[channel],
          label + " has the expected neutron birth");
    check_closure(kExpectedLoss[channel], kExpectedBirth[channel],
                  kExpectedNeutrons[channel], label);
}

void check_failure(const double* rates, int expected_status,
                   const std::string& label) {
    fusion_particle_sources_v1 source{
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        -1.0};
    const int status = fusion_c_particle_sources(rates, &source);
    check(status == expected_status, label + " returns the expected status");
    check(all_zero(source), label + " clears every output field");
}

}  // namespace

int main() {
    for (int channel = 0; channel < FUSION_CHANNEL_COUNT; ++channel) {
        check_channel(channel);
    }

    /* Known integer mixture: DD values are already event rates, including 1/2. */
    const double mixed_rates[FUSION_CHANNEL_COUNT]{2.0, 7.0, 11.0, 13.0,
                                                    17.0};
    fusion_particle_sources_v1 mixed{};
    int status = fusion_c_particle_sources(mixed_rates, &mixed);
    check(status == PB11_STATUS_OK, "mixed integer channels return OK");
    const std::array<double, FUSION_SPECIES_COUNT> expected_loss{{
        2.0, 66.0, 13.0, 17.0, 0.0, 2.0}};
    const std::array<double, FUSION_SPECIES_COUNT> expected_birth{{
        24.0, 0.0, 7.0, 11.0, 36.0, 0.0}};
    for (int species = 0; species < FUSION_SPECIES_COUNT; ++species) {
        check(mixed.reactant_loss[species] == expected_loss[species],
              "mixed channels have known reactant losses");
        check(mixed.product_birth[species] == expected_birth[species],
              "mixed channels have known product births");
        check(mixed.net_source[species] ==
                  expected_birth[species] - expected_loss[species],
              "mixed channels keep birth and loss separate");
    }
    check(mixed.reactant_loss[FUSION_DEUTERON] == 66.0,
          "DD event rates are not halved a second time");
    check(mixed.neutron_birth == 24.0,
          "mixed channels count both neutron-producing branches");
    check_closure(expected_loss, expected_birth, mixed.neutron_birth,
                  "mixed channels");

    const double zero_rates[FUSION_CHANNEL_COUNT]{};
    fusion_particle_sources_v1 zero_source{};
    status = fusion_c_particle_sources(zero_rates, &zero_source);
    check(status == PB11_STATUS_OK && all_zero(zero_source),
          "zero event rates produce an all-zero source");

    double nan_rates[FUSION_CHANNEL_COUNT]{};
    nan_rates[FUSION_PB11_3ALPHA] =
        std::numeric_limits<double>::quiet_NaN();
    check_failure(nan_rates, PB11_STATUS_INVALID_ARGUMENT, "NaN event rate");

    double inf_rates[FUSION_CHANNEL_COUNT]{};
    inf_rates[FUSION_DD_TP] = std::numeric_limits<double>::infinity();
    check_failure(inf_rates, PB11_STATUS_INVALID_ARGUMENT,
                  "infinite event rate");

    double negative_rates[FUSION_CHANNEL_COUNT]{};
    negative_rates[FUSION_DT_ALPHAN] = -1.0;
    check_failure(negative_rates, PB11_STATUS_OUT_OF_RANGE,
                  "negative event rate");

    double overflow_rates[FUSION_CHANNEL_COUNT]{};
    overflow_rates[FUSION_PB11_3ALPHA] =
        std::numeric_limits<double>::max();
    check_failure(overflow_rates, PB11_STATUS_NUMERICAL_FAILURE,
                  "unrepresentable product birth");

    fusion_particle_sources_v1 null_rates_source{
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        -1.0};
    status = fusion_c_particle_sources(nullptr, &null_rates_source);
    check(status == PB11_STATUS_INVALID_ARGUMENT &&
              all_zero(null_rates_source),
          "null input rates are rejected and clear output");
    check(fusion_c_particle_sources(zero_rates, nullptr) ==
              PB11_STATUS_NULL_OUTPUT,
          "null output pointer is reported");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All fusion-network tests passed\n";
    return 0;
}
