#include "fusion_thermal_burn.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using Species = std::array<double, FUSION_SPECIES_COUNT>;
using Channels = std::array<double, FUSION_CHANNEL_COUNT>;

constexpr double kEnergyScale = 1.0e-16;

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
    }
    assert(condition);
    if (!condition) std::abort(); // checks remain active in Release builds
}

bool close_relative(double actual, double expected, double tolerance,
                    double absolute_floor = 0.0) {
    if (actual == expected) {
        return true;
    }
    if (!std::isfinite(actual) || !std::isfinite(expected)) {
        return false;
    }
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return std::abs(actual - expected) <=
           std::max(absolute_floor, tolerance * scale);
}

void clear_output(fusion_thermal_burn_v1& out, double value) {
    for (int c = 0; c < FUSION_CHANNEL_COUNT; ++c) {
        out.events_m3[c] = value;
    }
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        out.reactant_removed_m3[s] = value;
        out.reactant_removed_energy_J_m3[s] = value;
        out.fast_product_birth_m3[s] = value;
        out.number_residual_m3[s] = value;
        out.energy_residual_J_m3[s] = value;
    }
    out.neutron_birth_m3 = value;
}

bool output_is_zero(const fusion_thermal_burn_v1& out) {
    for (int c = 0; c < FUSION_CHANNEL_COUNT; ++c) {
        if (out.events_m3[c] != 0.0) {
            return false;
        }
    }
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        if (out.reactant_removed_m3[s] != 0.0 ||
            out.reactant_removed_energy_J_m3[s] != 0.0 ||
            out.fast_product_birth_m3[s] != 0.0 ||
            out.number_residual_m3[s] != 0.0 ||
            out.energy_residual_J_m3[s] != 0.0) {
            return false;
        }
    }
    return out.neutron_birth_m3 == 0.0;
}

void assert_finite_output(const fusion_thermal_burn_v1& out,
                          const std::string& label) {
    for (int c = 0; c < FUSION_CHANNEL_COUNT; ++c) {
        require(std::isfinite(out.events_m3[c]) && out.events_m3[c] >= 0.0,
                label + " event is finite and nonnegative");
    }
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        require(std::isfinite(out.reactant_removed_m3[s]) &&
                    out.reactant_removed_m3[s] >= 0.0,
                label + " removed number is finite and nonnegative");
        require(std::isfinite(out.reactant_removed_energy_J_m3[s]) &&
                    out.reactant_removed_energy_J_m3[s] >= 0.0,
                label + " removed energy is finite and nonnegative");
        require(std::isfinite(out.fast_product_birth_m3[s]) &&
                    out.fast_product_birth_m3[s] >= 0.0,
                label + " product birth is finite and nonnegative");
        require(std::isfinite(out.number_residual_m3[s]),
                label + " number residual is finite");
        require(std::isfinite(out.energy_residual_J_m3[s]),
                label + " energy residual is finite");
    }
    require(std::isfinite(out.neutron_birth_m3) && out.neutron_birth_m3 >= 0.0,
            label + " neutron birth is finite and nonnegative");
}

bool same_output(const fusion_thermal_burn_v1& lhs,
                 const fusion_thermal_burn_v1& rhs) {
    for (int c = 0; c < FUSION_CHANNEL_COUNT; ++c) {
        if (lhs.events_m3[c] != rhs.events_m3[c]) {
            return false;
        }
    }
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        if (lhs.reactant_removed_m3[s] != rhs.reactant_removed_m3[s] ||
            lhs.reactant_removed_energy_J_m3[s] !=
                rhs.reactant_removed_energy_J_m3[s] ||
            lhs.fast_product_birth_m3[s] != rhs.fast_product_birth_m3[s] ||
            lhs.number_residual_m3[s] != rhs.number_residual_m3[s] ||
            lhs.energy_residual_J_m3[s] != rhs.energy_residual_J_m3[s]) {
            return false;
        }
    }
    return lhs.neutron_birth_m3 == rhs.neutron_birth_m3;
}

int call_burn(double dt, const Species& old_number, const Species& old_energy,
              const Channels& reactivity, const Channels& mean_a,
              const Channels& mean_b, Species& trial_number,
              Species& trial_energy, fusion_thermal_burn_v1& out) {
    return fusion_c_thermal_burn_trial(
        dt, old_number.data(), old_energy.data(), reactivity.data(),
        mean_a.data(), mean_b.data(), trial_number.data(), trial_energy.data(),
        &out);
}

void assert_state_close(const Species& actual, const Species& expected,
                        double tolerance, const std::string& label) {
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        require(close_relative(actual[s], expected[s], tolerance),
                label + " species state");
    }
}

void assert_state_nonnegative(const Species& number, const Species& energy,
                             const std::string& label) {
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        require(std::isfinite(number[s]) && number[s] >= 0.0,
                label + " number is finite and nonnegative");
        require(std::isfinite(energy[s]) && energy[s] >= 0.0,
                label + " energy is finite and nonnegative");
    }
}

void assert_reactant_residuals(const Species& old_number,
                               const Species& old_energy,
                               const Species& trial_number,
                               const Species& trial_energy,
                               const fusion_thermal_burn_v1& out,
                               const std::string& label) {
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        const double number_scale =
            std::max(1.0, std::abs(old_number[s]) +
                               std::abs(trial_number[s]) +
                               std::abs(out.reactant_removed_m3[s]));
        const double energy_scale =
            std::max(1.0, std::abs(old_energy[s]) +
                               std::abs(trial_energy[s]) +
                               std::abs(out.reactant_removed_energy_J_m3[s]));
        require(std::abs(trial_number[s] - old_number[s] +
                         out.reactant_removed_m3[s]) <=
                    2.0e-11 * number_scale,
                label + " number balance closes");
        require(std::abs(trial_energy[s] - old_energy[s] +
                         out.reactant_removed_energy_J_m3[s]) <=
                    2.0e-11 * energy_scale,
                label + " energy balance closes");
    }
}

void assert_network_stoichiometry(double dt,
                                  const fusion_thermal_burn_v1& out,
                                  const std::string& label) {
    fusion_particle_sources_v1 sources{};
    Channels rates{};
    for(int c=0;c<FUSION_CHANNEL_COUNT;++c)rates[c]=out.events_m3[c]/dt;
    const int status = fusion_c_particle_sources(rates.data(), &sources);
    require(status == PB11_STATUS_OK, label + " network stoichiometry returns OK");
    for (int s = 0; s < FUSION_SPECIES_COUNT; ++s) {
        require(close_relative(out.reactant_removed_m3[s],
                               dt * sources.reactant_loss[s], 2.0e-12,
                               1.0e-30),
                label + " reactant stoichiometry");
        require(close_relative(out.fast_product_birth_m3[s],
                               dt * sources.product_birth[s], 2.0e-12,
                               1.0e-30),
                label + " fast-product stoichiometry");
    }
    require(close_relative(out.neutron_birth_m3, dt * sources.neutron_birth,
                           2.0e-12, 1.0e-30),
            label + " neutron stoichiometry");
}

void zero_inputs(Species& old_number, Species& old_energy,
                 Channels& reactivity, Channels& mean_a, Channels& mean_b) {
    old_number.fill(0.0);
    old_energy.fill(0.0);
    reactivity.fill(0.0);
    mean_a.fill(0.0);
    mean_b.fill(0.0);
}

void test_zero_reactivity_unchanged() {
    Species old_number{{1.0e20, 2.0e20, 3.0e20, 4.0e20, 5.0e20, 6.0e20}};
    Species old_energy{{1.0e5, 2.0e5, 3.0e5, 4.0e5, 5.0e5, 6.0e5}};
    const Species old_number_copy = old_number;
    const Species old_energy_copy = old_energy;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);

    const int status = call_burn(0.25, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_OK, "zero-K burn returns OK");
    assert_state_close(trial_number, old_number, 0.0,
                       "zero-K number state is unchanged");
    assert_state_close(trial_energy, old_energy, 0.0,
                       "zero-K energy state is unchanged");
    require(output_is_zero(out), "zero-K output ledger is zero");
    require(old_number == old_number_copy && old_energy == old_energy_copy,
            "zero-K old state remains unchanged");
}

void test_pure_pb11_quadratic() {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_PROTON] = 1.0e20;
    old_number[FUSION_BORON11] = 1.0e20;
    old_energy[FUSION_PROTON] = 7.0e5;
    old_energy[FUSION_BORON11] = 8.0e5;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_PB11_3ALPHA] = 1.0e-20;
    mean_a[FUSION_PB11_3ALPHA] = 0.7e-16;
    mean_b[FUSION_PB11_3ALPHA] = 2.6e-16;
    constexpr double dt = 0.3;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);

    const int status = call_burn(dt, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_OK, "pure pB11 burn returns OK");
    const long double stiffness =
        static_cast<long double>(dt) * reactivity[FUSION_PB11_3ALPHA];
    const long double n0 = old_number[FUSION_PROTON];
    const long double expected_new =
        2.0L * n0 /
        (1.0L + std::sqrt(1.0L + 4.0L * stiffness * n0));
    const long double expected_event =
        stiffness * expected_new * expected_new;
    require(close_relative(trial_number[FUSION_PROTON],
                           static_cast<double>(expected_new), 2.0e-12),
            "pure pB11 proton state matches quadratic root");
    require(close_relative(trial_number[FUSION_BORON11],
                           static_cast<double>(expected_new), 2.0e-12),
            "pure pB11 boron state matches quadratic root");
    require(close_relative(out.events_m3[FUSION_PB11_3ALPHA],
                           static_cast<double>(expected_event), 2.0e-12),
            "pure pB11 event matches quadratic root");
    require(close_relative(out.reactant_removed_m3[FUSION_PROTON],
                           static_cast<double>(expected_event), 2.0e-12),
            "pB11 consumes one proton per event");
    require(close_relative(out.reactant_removed_m3[FUSION_BORON11],
                           static_cast<double>(expected_event), 2.0e-12),
            "pB11 consumes one boron per event");
    require(close_relative(out.fast_product_birth_m3[FUSION_HELIUM4],
                           3.0 * static_cast<double>(expected_event),
                           2.0e-12),
            "pB11 births three fast alphas per event");
    require(out.neutron_birth_m3 == 0.0,
            "pB11 has no neutron birth in this channel");
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_PROTON],
                           static_cast<double>(expected_event) *
                               mean_a[FUSION_PB11_3ALPHA],
                           2.0e-12),
            "pB11 proton debit uses supplied conditional mean");
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_BORON11],
                           static_cast<double>(expected_event) *
                               mean_b[FUSION_PB11_3ALPHA],
                           2.0e-12),
            "pB11 boron debit uses supplied conditional mean");
    assert_state_nonnegative(trial_number, trial_energy, "pure pB11");
    assert_finite_output(out, "pure pB11");
    assert_network_stoichiometry(dt, out, "pure pB11");
    assert_reactant_residuals(old_number, old_energy, trial_number, trial_energy,
                              out, "pure pB11");
}

void test_pure_dd_two_branches() {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_DEUTERON] = 1.0e21;
    old_energy[FUSION_DEUTERON] = 1.0e7;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_DD_TP] = 1.0e-20;
    reactivity[FUSION_DD_HE3N] = 3.0e-20;
    mean_a[FUSION_DD_TP] = 0.8e-16;
    mean_b[FUSION_DD_TP] = 1.7e-16;
    mean_a[FUSION_DD_HE3N] = 2.1e-16;
    mean_b[FUSION_DD_HE3N] = 3.2e-16;
    constexpr double dt = 0.1;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);

    const int status = call_burn(dt, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_OK, "pure DD burn returns OK");
    const long double k_total = reactivity[FUSION_DD_TP] +
                                reactivity[FUSION_DD_HE3N];
    const long double a = static_cast<long double>(dt) * k_total;
    const long double d0 = old_number[FUSION_DEUTERON];
    const long double d_new =
        2.0L * d0 / (1.0L + std::sqrt(1.0L + 4.0L * a * d0));
    const long double e_tp =
        0.5L * dt * reactivity[FUSION_DD_TP] * d_new * d_new;
    const long double e_he3n =
        0.5L * dt * reactivity[FUSION_DD_HE3N] * d_new * d_new;
    require(close_relative(trial_number[FUSION_DEUTERON],
                           static_cast<double>(d_new), 2.0e-12),
            "pure DD deuteron state matches quadratic root");
    require(close_relative(out.events_m3[FUSION_DD_TP],
                           static_cast<double>(e_tp), 2.0e-12),
            "DD t+p event matches BE root");
    require(close_relative(out.events_m3[FUSION_DD_HE3N],
                           static_cast<double>(e_he3n), 2.0e-12),
            "DD He3+n event matches BE root");
    require(close_relative(out.events_m3[FUSION_DD_TP] /
                               out.events_m3[FUSION_DD_HE3N],
                           reactivity[FUSION_DD_TP] /
                               reactivity[FUSION_DD_HE3N],
                           2.0e-12),
            "DD branch event ratio follows frozen reactivities");
    const double total_dd = out.events_m3[FUSION_DD_TP] +
                            out.events_m3[FUSION_DD_HE3N];
    require(close_relative(out.reactant_removed_m3[FUSION_DEUTERON],
                           2.0 * total_dd, 2.0e-12),
            "DD consumes two deuterons per event");
    require(close_relative(out.fast_product_birth_m3[FUSION_PROTON],
                           out.events_m3[FUSION_DD_TP], 2.0e-12),
            "DD t+p branch births one proton");
    require(close_relative(out.fast_product_birth_m3[FUSION_TRITON],
                           out.events_m3[FUSION_DD_TP], 2.0e-12),
            "DD t+p branch births one triton");
    require(close_relative(out.fast_product_birth_m3[FUSION_HELIUM3],
                           out.events_m3[FUSION_DD_HE3N], 2.0e-12),
            "DD He3+n branch births one helium-3");
    require(close_relative(out.neutron_birth_m3,
                           out.events_m3[FUSION_DD_HE3N], 2.0e-12),
            "DD He3+n branch births one neutron");
    const double expected_deuteron_energy =
        out.events_m3[FUSION_DD_TP] *
            (mean_a[FUSION_DD_TP] + mean_b[FUSION_DD_TP]) +
        out.events_m3[FUSION_DD_HE3N] *
            (mean_a[FUSION_DD_HE3N] + mean_b[FUSION_DD_HE3N]);
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_DEUTERON],
                           expected_deuteron_energy, 2.0e-12),
            "both DD branch means debit the same deuteron pool");
    assert_state_nonnegative(trial_number, trial_energy, "pure DD");
    assert_finite_output(out, "pure DD");
    assert_network_stoichiometry(dt, out, "pure DD");
    assert_reactant_residuals(old_number, old_energy, trial_number, trial_energy,
                              out, "pure DD");
}

void test_dt_dhe3_shared_deuteron_competition() {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_DEUTERON] = 5.0e20;
    old_number[FUSION_TRITON] = 2.0e20;
    old_number[FUSION_HELIUM3] = 3.0e20;
    old_energy[FUSION_DEUTERON] = 8.0e6;
    old_energy[FUSION_TRITON] = 6.0e6;
    old_energy[FUSION_HELIUM3] = 7.0e6;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_DT_ALPHAN] = 2.0e-21;
    reactivity[FUSION_DHE3_ALPHAP] = 3.0e-21;
    mean_a[FUSION_DT_ALPHAN] = 0.9e-16;
    mean_b[FUSION_DT_ALPHAN] = 1.8e-16;
    mean_a[FUSION_DHE3_ALPHAP] = 2.2e-16;
    mean_b[FUSION_DHE3_ALPHAP] = 3.4e-16;
    constexpr double dt = 0.2;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);

    const int status = call_burn(dt, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_OK,
            "coupled DT/DHe3 burn returns OK");
    const double d_new = trial_number[FUSION_DEUTERON];
    const double t_new = trial_number[FUSION_TRITON];
    const double he3_new = trial_number[FUSION_HELIUM3];
    const double dt_event = out.events_m3[FUSION_DT_ALPHAN];
    const double dhe3_event = out.events_m3[FUSION_DHE3_ALPHAP];
    require(dt_event > 0.0 && dhe3_event > 0.0,
            "both coupled DT/DHe3 branches burn");
    require(close_relative(dt_event,
                           dt * reactivity[FUSION_DT_ALPHAN] * d_new * t_new,
                           2.0e-11),
            "DT event uses the common final deuteron and triton");
    require(close_relative(
                dhe3_event,
                dt * reactivity[FUSION_DHE3_ALPHAP] * d_new * he3_new,
                2.0e-11),
            "DHe3 event uses the common final deuteron and helium-3");
    require(close_relative(old_number[FUSION_DEUTERON] - d_new,
                           dt_event + dhe3_event, 2.0e-11),
            "DT and DHe3 compete for one shared deuteron pool");
    require(close_relative(old_number[FUSION_TRITON] - t_new, dt_event,
                           2.0e-11),
            "DT consumes one triton per event");
    require(close_relative(old_number[FUSION_HELIUM3] - he3_new, dhe3_event,
                           2.0e-11),
            "DHe3 consumes one helium-3 per event");
    require(close_relative(out.fast_product_birth_m3[FUSION_HELIUM4],
                           dt_event + dhe3_event, 2.0e-11),
            "both coupled branches birth one fast helium-4");
    require(close_relative(out.fast_product_birth_m3[FUSION_PROTON],
                           dhe3_event, 2.0e-11),
            "DHe3 branch births one fast proton");
    require(close_relative(out.neutron_birth_m3, dt_event, 2.0e-11),
            "DT branch births one neutron");
    const double expected_d_energy =
        dt_event * mean_a[FUSION_DT_ALPHAN] +
        dhe3_event * mean_a[FUSION_DHE3_ALPHAP];
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_DEUTERON],
                           expected_d_energy, 2.0e-11),
            "shared D energy debit sums both conditioned projectile means");
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_TRITON],
                           dt_event * mean_b[FUSION_DT_ALPHAN], 2.0e-11),
            "DT target energy debit uses supplied mean");
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_HELIUM3],
                           dhe3_event * mean_b[FUSION_DHE3_ALPHAP], 2.0e-11),
            "DHe3 target energy debit uses supplied mean");
    assert_state_nonnegative(trial_number, trial_energy, "coupled DT/DHe3");
    assert_finite_output(out, "coupled DT/DHe3");
    assert_network_stoichiometry(dt, out, "coupled DT/DHe3");
    assert_reactant_residuals(old_number, old_energy, trial_number, trial_energy,
                              out, "coupled DT/DHe3");
}

void test_small_dt_network_limit() {
    Species old_number{{1.4e20, 1.1e20, 0.9e20, 0.8e20, 0.0, 1.3e20}};
    Species old_energy{{3.0e6, 3.0e6, 3.0e6, 3.0e6, 0.0, 3.0e6}};
    Channels reactivity{{1.0e-20, 2.0e-20, 3.0e-20, 4.0e-20, 5.0e-20}};
    Channels mean_a{{0.4e-16, 0.7e-16, 1.0e-16, 1.3e-16, 1.6e-16}};
    Channels mean_b{{0.6e-16, 0.9e-16, 1.2e-16, 1.5e-16, 1.8e-16}};
    constexpr double dt = 1.0e-6;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);
    const int status = call_burn(dt, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_OK, "small-dt five-channel burn returns OK");

    Channels instantaneous{};
    instantaneous[FUSION_PB11_3ALPHA] =
        reactivity[FUSION_PB11_3ALPHA] * old_number[FUSION_PROTON] *
        old_number[FUSION_BORON11];
    instantaneous[FUSION_DD_TP] =
        0.5 * reactivity[FUSION_DD_TP] *
        old_number[FUSION_DEUTERON] * old_number[FUSION_DEUTERON];
    instantaneous[FUSION_DD_HE3N] =
        0.5 * reactivity[FUSION_DD_HE3N] *
        old_number[FUSION_DEUTERON] * old_number[FUSION_DEUTERON];
    instantaneous[FUSION_DT_ALPHAN] =
        reactivity[FUSION_DT_ALPHAN] * old_number[FUSION_DEUTERON] *
        old_number[FUSION_TRITON];
    instantaneous[FUSION_DHE3_ALPHAP] =
        reactivity[FUSION_DHE3_ALPHAP] * old_number[FUSION_DEUTERON] *
        old_number[FUSION_HELIUM3];
    fusion_particle_sources_v1 initial_sources{};
    require(fusion_c_particle_sources(instantaneous.data(), &initial_sources)==PB11_STATUS_OK,
            "instantaneous sources available");
    constexpr int reactant_a[5]={0,1,1,1,1},reactant_b[5]={5,1,1,2,3};
    double error=0;
    for (int c = 0; c < FUSION_CHANNEL_COUNT; ++c) {
        // Depletion only: implicit losses cannot exceed dt times initial
        // losses. Thus relative event-rate change is bounded by fa+fb.
        const int a=reactant_a[c],b=reactant_b[c];
        const double bound=dt*(initial_sources.reactant_loss[a]/old_number[a]
                            +initial_sources.reactant_loss[b]/old_number[b]);
        const double relative=std::abs(out.events_m3[c]/(dt*instantaneous[c])-1);
        require(relative<=bound+2e-15,"instantaneous-limit depletion bound");
        error=std::max(error,relative);
    }
    Species half_n,half_u;fusion_thermal_burn_v1 half{};
    require(call_burn(dt/2,old_number,old_energy,reactivity,mean_a,mean_b,
                     half_n,half_u,half)==PB11_STATUS_OK,"half-step limit call");
    double half_error=0;
    for(int c=0;c<FUSION_CHANNEL_COUNT;++c)
        half_error=std::max(half_error,std::abs(half.events_m3[c]/(dt/2*instantaneous[c])-1));
    require(error/half_error>1.99&&error/half_error<2.01,
            "instantaneous-rate error decreases linearly with dt");
    assert_network_stoichiometry(dt, out, "small-dt five-channel");
    assert_reactant_residuals(old_number, old_energy, trial_number, trial_energy,
                              out, "small-dt five-channel");
    assert_state_nonnegative(trial_number, trial_energy, "small-dt five-channel");
}

double run_pure_dd(double dt, int steps) {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_DEUTERON] = 1.0e20;
    old_energy[FUSION_DEUTERON] = 1.0e6;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_DD_TP] = 0.4e-20;
    reactivity[FUSION_DD_HE3N] = 0.6e-20;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    for (int step = 0; step < steps; ++step) {
        clear_output(out, -1.0);
        const int status = call_burn(dt, old_number, old_energy, reactivity,
                                     mean_a, mean_b, trial_number, trial_energy,
                                     out);
        require(status == PB11_STATUS_OK,
                "fixed-coefficient DD time step returns OK");
        assert_state_nonnegative(trial_number, trial_energy, "fixed DD step");
        assert_reactant_residuals(old_number, old_energy, trial_number,
                                  trial_energy, out, "fixed DD step");
        old_number = trial_number;
        old_energy = trial_energy;
    }
    return old_number[FUSION_DEUTERON];
}

void test_dd_time_refinement() {
    constexpr double d0 = 1.0e20;
    constexpr double k_total = 1.0e-20;
    constexpr double final_time = 1.0;
    const double exact = d0 / (1.0 + k_total * d0 * final_time);
    const double coarse = run_pure_dd(0.1, 10);
    const double fine = run_pure_dd(0.05, 20);
    const double coarse_error = std::abs(coarse - exact);
    const double fine_error = std::abs(fine - exact);
    require(fine_error < coarse_error,
            "DD implicit time refinement reduces the continuous-solution error");
    std::cout << "DD refinement errors " << coarse_error/d0 << " " << fine_error/d0
              << " ratio " << coarse_error/fine_error << "\n";
    require(coarse_error / fine_error > 1.7 && coarse_error / fine_error < 2.3,
            "DD fixed-coefficient convergence is first order");
}

void test_dynamic_range() {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_PROTON] = 1.0e20;
    old_number[FUSION_BORON11] = 1.0e20;
    old_energy[FUSION_PROTON] = 1.0e3;
    old_energy[FUSION_BORON11] = 1.0e3;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_PB11_3ALPHA] = 1.0e-40;
    mean_a[FUSION_PB11_3ALPHA] = 1.0e-16;
    mean_b[FUSION_PB11_3ALPHA] = 1.0e-16;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);
    int status = call_burn(1.0, old_number, old_energy, reactivity, mean_a,
                           mean_b, trial_number, trial_energy, out);
    require(status == PB11_STATUS_OK, "tiny representable burn returns OK");
    require(std::isfinite(out.events_m3[FUSION_PB11_3ALPHA]) &&
                out.events_m3[FUSION_PB11_3ALPHA] >= 0.0,
            "tiny burn event remains representable or rounds to zero");

    zero_inputs(old_number, old_energy, reactivity, mean_a, mean_b);
    old_number[FUSION_DEUTERON] = 1.0e20;
    old_energy[FUSION_DEUTERON] = 1.0e6;
    reactivity[FUSION_DD_TP] = 1.0e-18;
    reactivity[FUSION_DD_HE3N] = 1.0e-18;
    clear_output(out, -1.0);
    status = call_burn(1.0e4, old_number, old_energy, reactivity, mean_a,
                       mean_b, trial_number, trial_energy, out);
    require(status == PB11_STATUS_OK, "strong but representable burn returns OK");
    require(trial_number[FUSION_DEUTERON] > 0.0 &&
                trial_number[FUSION_DEUTERON] < old_number[FUSION_DEUTERON],
            "strong burn retains a positive representable survivor");
    assert_state_nonnegative(trial_number, trial_energy, "strong burn");
    assert_reactant_residuals(old_number, old_energy, trial_number, trial_energy,
                              out, "strong burn");
}

void test_energy_conditioned_means() {
    /* This is a focused check that no hidden 3*kT/2 substitution occurs. */
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_PROTON] = 1.0e20;
    old_number[FUSION_BORON11] = 1.0e20;
    old_energy[FUSION_PROTON] = 1.0e6;
    old_energy[FUSION_BORON11] = 1.0e6;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_PB11_3ALPHA] = 2.0e-20;
    mean_a[FUSION_PB11_3ALPHA] = 0.11e-16;
    mean_b[FUSION_PB11_3ALPHA] = 0.37e-16;
    constexpr double dt = 0.2;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);
    const int status = call_burn(dt, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_OK, "conditioned-mean burn returns OK");
    const double event = out.events_m3[FUSION_PB11_3ALPHA];
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_PROTON],
                           event * mean_a[FUSION_PB11_3ALPHA], 2.0e-12),
            "proton energy uses provided conditioned mean");
    require(close_relative(out.reactant_removed_energy_J_m3[FUSION_BORON11],
                           event * mean_b[FUSION_PB11_3ALPHA], 2.0e-12),
            "boron energy uses provided conditioned mean");
    require(std::abs(out.reactant_removed_energy_J_m3[FUSION_PROTON] -
                     event * 1.5 * kEnergyScale) >
                1.0e-4 * std::max(1.0, event * kEnergyScale),
            "conditioned mean is observably distinct from 1.5 kT");
    assert_reactant_residuals(old_number, old_energy, trial_number, trial_energy,
                              out, "conditioned means");
}

void test_energy_insufficient_rejects_whole_step() {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_PROTON] = 1.0e20;
    old_number[FUSION_BORON11] = 1.0e20;
    old_energy[FUSION_PROTON] = 1.0;
    old_energy[FUSION_BORON11] = 1.0;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_PB11_3ALPHA] = 1.0e-20;
    mean_a[FUSION_PB11_3ALPHA] = 1.0;
    mean_b[FUSION_PB11_3ALPHA] = 1.0;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);
    const int status = call_burn(1.0, old_number, old_energy, reactivity,
                                 mean_a, mean_b, trial_number, trial_energy,
                                 out);
    require(status == PB11_STATUS_NUMERICAL_FAILURE,
            "energy-underflow trial rejects as numerical failure");
    require(std::all_of(trial_number.begin(), trial_number.end(),
                        [](double value) { return value == 0.0; }) &&
                std::all_of(trial_energy.begin(), trial_energy.end(),
                            [](double value) { return value == 0.0; }) &&
                output_is_zero(out),
            "energy-underflow rejection clears the whole trial");
    require(old_number[FUSION_PROTON] == 1.0e20 &&
                old_number[FUSION_BORON11] == 1.0e20 &&
                old_energy[FUSION_PROTON] == 1.0 &&
                old_energy[FUSION_BORON11] == 1.0,
            "energy-underflow rejection leaves old state unchanged");
}

void test_empty_and_rollback() {
    Species old_number{};
    Species old_energy{};
    Channels reactivity{{1.0e-20, 2.0e-20, 3.0e-20, 4.0e-20, 5.0e-20}};
    Channels mean_a{{0.5e-16, 0.5e-16, 0.5e-16, 0.5e-16, 0.5e-16}};
    Channels mean_b{{0.6e-16, 0.6e-16, 0.6e-16, 0.6e-16, 0.6e-16}};
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    clear_output(out, -1.0);
    int status = call_burn(0.5, old_number, old_energy, reactivity, mean_a,
                           mean_b, trial_number, trial_energy, out);
    require(status == PB11_STATUS_OK, "empty fuel burn returns OK");
    require(std::all_of(trial_number.begin(), trial_number.end(),
                        [](double value) { return value == 0.0; }) &&
                std::all_of(trial_energy.begin(), trial_energy.end(),
                            [](double value) { return value == 0.0; }) &&
                output_is_zero(out),
            "empty fuel remains exactly empty");

    old_number = Species{{1.0e20, 2.0e20, 3.0e20, 4.0e20, 0.0, 5.0e20}};
    old_energy = Species{{5.0e5, 5.0e5, 5.0e5, 5.0e5, 0.0, 5.0e5}};
    const Species old_number_copy = old_number;
    const Species old_energy_copy = old_energy;
    reactivity = Channels{{1.0e-22, 2.0e-22, 3.0e-22, 4.0e-22, 5.0e-22}};
    mean_a = Channels{{0.3e-16, 0.4e-16, 0.5e-16, 0.6e-16, 0.7e-16}};
    mean_b = Channels{{0.8e-16, 0.9e-16, 1.0e-16, 1.1e-16, 1.2e-16}};
    Species first_trial_number;
    Species first_trial_energy;
    Species second_trial_number;
    Species second_trial_energy;
    fusion_thermal_burn_v1 first_out;
    fusion_thermal_burn_v1 second_out;
    clear_output(first_out, -1.0);
    clear_output(second_out, -1.0);
    const int first_status =
        call_burn(0.1, old_number, old_energy, reactivity, mean_a, mean_b,
                  first_trial_number, first_trial_energy, first_out);
    const int second_status =
        call_burn(0.1, old_number, old_energy, reactivity, mean_a, mean_b,
                  second_trial_number, second_trial_energy, second_out);
    require(first_status == PB11_STATUS_OK && second_status == PB11_STATUS_OK,
            "repeated burn trials return OK");
    require(first_trial_number == second_trial_number &&
                first_trial_energy == second_trial_energy &&
                same_output(first_out, second_out),
            "repeated burn trials are identical for rollback");
    require(old_number == old_number_copy && old_energy == old_energy_copy,
            "repeated rollback trials leave old state unchanged");
}

void assert_failed_and_cleared(int status, int expected_status,
                               const Species& trial_number,
                               const Species& trial_energy,
                               const fusion_thermal_burn_v1& out,
                               const std::string& label) {
    require(status == expected_status, label + " returns expected status");
    require(std::all_of(trial_number.begin(), trial_number.end(),
                        [](double value) { return value == 0.0; }) &&
                std::all_of(trial_energy.begin(), trial_energy.end(),
                            [](double value) { return value == 0.0; }) &&
                output_is_zero(out),
            label + " clears all nonnull outputs");
}

void test_invalid_nonfinite_negative_null_overflow() {
    Species old_number{};
    Species old_energy{};
    old_number[FUSION_PROTON] = 1.0e20;
    old_number[FUSION_BORON11] = 1.0e20;
    old_energy[FUSION_PROTON] = 1.0e6;
    old_energy[FUSION_BORON11] = 1.0e6;
    Channels reactivity{};
    Channels mean_a{};
    Channels mean_b{};
    reactivity[FUSION_PB11_3ALPHA] = 1.0e-20;
    mean_a[FUSION_PB11_3ALPHA] = 0.2e-16;
    mean_b[FUSION_PB11_3ALPHA] = 0.3e-16;
    Species trial_number;
    Species trial_energy;
    fusion_thermal_burn_v1 out;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    old_number[FUSION_PROTON] = nan;
    clear_output(out, -1.0);
    int status = call_burn(0.1, old_number, old_energy, reactivity, mean_a,
                           mean_b, trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_INVALID_ARGUMENT,
                               trial_number, trial_energy, out,
                               "nonfinite old number");
    old_number[FUSION_PROTON] = 1.0e20;

    old_energy[FUSION_PROTON] = -1.0;
    clear_output(out, -1.0);
    status = call_burn(0.1, old_number, old_energy, reactivity, mean_a,
                       mean_b, trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_INVALID_ARGUMENT,
                               trial_number, trial_energy, out,
                               "negative old energy");
    old_energy[FUSION_PROTON] = 1.0e6;

    reactivity[FUSION_PB11_3ALPHA] = -1.0;
    clear_output(out, -1.0);
    status = call_burn(0.1, old_number, old_energy, reactivity, mean_a, mean_b,
                       trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_INVALID_ARGUMENT,
                               trial_number, trial_energy, out,
                               "negative reactivity");
    reactivity[FUSION_PB11_3ALPHA] = 1.0e-20;

    mean_a[FUSION_PB11_3ALPHA] = nan;
    clear_output(out, -1.0);
    status = call_burn(0.1, old_number, old_energy, reactivity, mean_a, mean_b,
                       trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_INVALID_ARGUMENT,
                               trial_number, trial_energy, out,
                               "nonfinite conditional mean");
    mean_a[FUSION_PB11_3ALPHA] = 0.2e-16;

    reactivity[FUSION_PB11_3ALPHA] = 0.0;
    mean_a[FUSION_PB11_3ALPHA] = 1.0e-16;
    clear_output(out, -1.0);
    status = call_burn(0.1, old_number, old_energy, reactivity, mean_a, mean_b,
                       trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_INVALID_ARGUMENT,
                               trial_number, trial_energy, out,
                               "nonzero mean with zero reactivity");
    reactivity[FUSION_PB11_3ALPHA] = 1.0e-20;
    mean_a[FUSION_PB11_3ALPHA] = 0.2e-16;

    clear_output(out, -1.0);
    status = call_burn(nan, old_number, old_energy, reactivity, mean_a, mean_b,
                       trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_INVALID_ARGUMENT,
                               trial_number, trial_energy, out,
                               "nonfinite time step");

    clear_output(out, -1.0);
    status = fusion_c_thermal_burn_trial(
        0.1, old_number.data(), old_energy.data(), reactivity.data(),
        mean_a.data(), mean_b.data(), nullptr, trial_energy.data(), &out);
    require(status == PB11_STATUS_NULL_OUTPUT,
            "null number trial output returns NULL_OUTPUT");
    require(std::all_of(trial_energy.begin(), trial_energy.end(),
                        [](double value) { return value == 0.0; }) &&
                output_is_zero(out),
            "null number trial clears remaining outputs");

    Species null_trial_number;
    Species null_trial_energy;
    clear_output(out, -1.0);
    status = fusion_c_thermal_burn_trial(
        0.1, old_number.data(), old_energy.data(), reactivity.data(),
        mean_a.data(), mean_b.data(), null_trial_number.data(),
        null_trial_energy.data(), nullptr);
    require(status == PB11_STATUS_NULL_OUTPUT &&
                std::all_of(null_trial_number.begin(), null_trial_number.end(),
                            [](double value) { return value == 0.0; }) &&
                std::all_of(null_trial_energy.begin(), null_trial_energy.end(),
                            [](double value) { return value == 0.0; }),
            "null ledger returns NULL_OUTPUT and clears trial arrays");

    Species huge_number{};
    Species huge_energy{};
    huge_number[FUSION_PROTON] = DBL_MAX;
    huge_number[FUSION_BORON11] = DBL_MAX;
    huge_energy[FUSION_PROTON] = 1.0e6;
    huge_energy[FUSION_BORON11] = 1.0e6;
    Channels huge_k{};
    Channels huge_a{};
    Channels huge_b{};
    huge_k[FUSION_PB11_3ALPHA] = 1.0e-100;
    huge_a[FUSION_PB11_3ALPHA] = 1.0e-16;
    huge_b[FUSION_PB11_3ALPHA] = 1.0e-16;
    clear_output(out, -1.0);
    status = call_burn(1.0, huge_number, huge_energy, huge_k, huge_a, huge_b,
                       trial_number, trial_energy, out);
    assert_failed_and_cleared(status, PB11_STATUS_NUMERICAL_FAILURE,
                               trial_number, trial_energy, out,
                               "event overflow");
}

}  // namespace

int main() {
    test_zero_reactivity_unchanged();
    test_pure_pb11_quadratic();
    test_pure_dd_two_branches();
    test_dt_dhe3_shared_deuteron_competition();
    test_small_dt_network_limit();
    test_dd_time_refinement();
    test_dynamic_range();
    test_energy_conditioned_means();
    test_energy_insufficient_rejects_whole_step();
    test_empty_and_rollback();
    test_invalid_nonfinite_negative_null_overflow();
    std::cout << "All fusion thermal-burn tests passed\n";
    return 0;
}
