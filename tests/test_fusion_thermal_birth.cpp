#include "fusion_thermal_birth.h"
#include "fusion_nuclear_data.h"
#include "fusion_network.h"
#include "fusion_rate_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

constexpr double mev = 1.602176634e-13;
constexpr double kev = 1.602176634e-16;

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

bool close_scaled(double actual, double expected, double relative,
                  double absolute = 0.0) {
    if (!std::isfinite(actual) || !std::isfinite(expected))
        return false;
    return std::abs(actual - expected) <=
           absolute + relative * std::max(std::abs(actual), std::abs(expected));
}

std::vector<double> clipped_edges() {
    // The source energies extend above this grid for every two-body channel,
    // while pB also has a low-energy alpha tail.  The source portions outside
    // the center hull must stay in explicit below/above ledgers rather than
    // being clipped away.
    constexpr int cells = 8;
    std::vector<double> edges(cells + 1);
    for (int i = 0; i <= cells; ++i)
        edges[static_cast<std::size_t>(i)] = 2.0 * mev * i / cells;
    return edges;
}

fusion_thermal_birth_options_v1 options_for(double relative_max_MeV,
                                             int relative_order = 8,
                                             double cm_max_kT = 40.0) {
    fusion_thermal_birth_options_v1 options{};
    options.relative_max_J = relative_max_MeV * mev;
    options.cm_max_kT = cm_max_kT;
    options.ground_state_q_J = 91.84 * kev;
    // Use the exact lower-bound spelling accepted by the alpha-event API.
    options.cutoff_J = 0.001 * mev;
    options.l1_fraction = 0.76;
    options.relative_phase = 0.0;
    options.narrow_peak_fraction = 0.051;
    options.continuum_peak_scale = 1.0;
    options.continuation = FUSION_ENDPOINT_S;
    options.pb_low = FUSION_PB_LOW_TB;
    options.remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY;
    options.broad_mode = 13;
    options.fsci_policy = 0;
    options.relative_order = relative_order;
    options.cm_order = 8;
    // Keep the alpha event quadrature fixed.  Refinement below changes only
    // the reacting relative-energy quadrature.
    options.nq = 8;
    options.ncos = 8;
    return options;
}

bool finite_result(const fusion_thermal_birth_v1 &result) {
    const double scalar[] = {
        result.reactivity_m3_s,
        result.reference_reactivity_m3_s,
        result.product_energy_moment_J_m3_s,
        result.number_residual_m3_s,
        result.energy_residual_J_m3_s,
        result.relative_rate_discrepancy,
        result.relative_reactant_energy_discrepancy,
        result.cm_retained_probability,
        result.cm_tail_probability,
        result.cm_tail_energy_moment_J,
        result.max_cm_energy_shift_fraction,
        result.max_shell_remap_fraction,
    };
    for (double value : scalar)
        if (!std::isfinite(value))
            return false;
    for (double value : result.reactant_energy_moment_J_m3_s)
        if (!std::isfinite(value))
            return false;
    for (double value : result.reference_reactant_energy_moment_J_m3_s)
        if (!std::isfinite(value))
            return false;
    for (double value : result.below_number_m3_s)
        if (!std::isfinite(value))
            return false;
    for (double value : result.below_energy_J_m3_s)
        if (!std::isfinite(value))
            return false;
    for (double value : result.above_number_m3_s)
        if (!std::isfinite(value))
            return false;
    for (double value : result.above_energy_J_m3_s)
        if (!std::isfinite(value))
            return false;
    return true;
}

struct Evaluation {
    std::vector<double> edges;
    std::vector<double> birth;
    fusion_thermal_birth_v1 result{};
};

Evaluation evaluate(int channel, double temperature_keV, int relative_order = 8,
                    double cm_max_kT = 40.0) {
    Evaluation evaluation;
    evaluation.edges = clipped_edges();
    evaluation.birth.assign(7 * (evaluation.edges.size() - 1), 7.0);
    const double relative_max = channel == FUSION_PB11_3ALPHA ? 2.0 : 5.0;
    fusion_thermal_birth_options_v1 options =
        options_for(relative_max, relative_order, cm_max_kT);
    const int status = fusion_c_thermal_birth_grid(
        channel, temperature_keV * kev, &options,
        static_cast<int>(evaluation.edges.size() - 1), evaluation.edges.data(),
        evaluation.birth.data(), &evaluation.result);
    require(status == PB11_STATUS_OK, "thermal-birth evaluation status");
    return evaluation;
}

void check_conservation(const Evaluation &evaluation, int channel,
                        double temperature_keV) {
    const auto &edges = evaluation.edges;
    const auto &birth = evaluation.birth;
    const auto &result = evaluation.result;
    require(finite_result(result), "thermal-birth diagnostics are finite");
    require(result.reactivity_m3_s > 0.0 &&
                result.reference_reactivity_m3_s > 0.0,
            "thermal-birth rates are positive");
    require(result.cm_retained_probability > 0.0 &&
                result.cm_retained_probability <= 1.0 + 1e-12 &&
                result.cm_tail_probability > 0.0 &&
                result.cm_tail_probability < 1.0 &&
                close_scaled(result.cm_retained_probability +
                                 result.cm_tail_probability,
                             1.0, 3e-12, 2e-15),
            "finite CM cutoff probability closes");
    require(result.cm_tail_energy_moment_J > 0.0,
            "finite CM cutoff has a positive energy tail diagnostic");

    const double centers_count = static_cast<double>(edges.size() - 1);
    std::array<double, 7> number{};
    std::array<double, 7> energy{};
    for (std::size_t species = 0; species < 7; ++species) {
        for (std::size_t cell = 0; cell + 1 < edges.size(); ++cell) {
            const double value = birth[species * static_cast<std::size_t>(centers_count) + cell];
            require(std::isfinite(value) && value >= 0.0,
                    "birth coefficients are finite and nonnegative");
            number[species] += value;
            energy[species] +=
                0.5 * (edges[cell] + edges[cell + 1]) * value;
        }
        number[species] += result.below_number_m3_s[species] +
                           result.above_number_m3_s[species];
        energy[species] += result.below_energy_J_m3_s[species] +
                           result.above_energy_J_m3_s[species];
        require(std::isfinite(number[species]) && number[species] >= 0.0 &&
                    std::isfinite(energy[species]) && energy[species] >= 0.0,
                "mapped and spill ledgers are finite");
    }

    fusion_nuclear_channel_v1 reaction{};
    require(fusion_c_nuclear_channel(channel, &reaction) == PB11_STATUS_OK,
            "channel metadata status");
    std::array<int, 7> multiplicity{};
    for (int product = 0; product < reaction.product_count; ++product)
        ++multiplicity[static_cast<std::size_t>(reaction.product_ids[product])];
    double total_number = 0.0;
    for (int species = 0; species < 7; ++species) {
        total_number += number[static_cast<std::size_t>(species)];
        require(close_scaled(number[static_cast<std::size_t>(species)],
                             multiplicity[static_cast<std::size_t>(species)] *
                                 result.reactivity_m3_s,
                             3e-9, 1e-300),
                "product number closes including spills");
    }
    require(close_scaled(total_number,
                         reaction.product_count * result.reactivity_m3_s,
                         3e-9, 1e-300),
            "total product number closes including spills");

    double total_energy = 0.0;
    for (double value : energy)
        total_energy += value;
    require(close_scaled(total_energy, result.product_energy_moment_J_m3_s,
                         4e-9, 1e-300),
            "product energy closes including spills");
    require(close_scaled(result.number_residual_m3_s, 0.0, 0.0,
                         3e-9 * result.reactivity_m3_s),
            "reported product number residual is small");
    require(close_scaled(result.energy_residual_J_m3_s, 0.0, 0.0,
                         4e-9 * result.product_energy_moment_J_m3_s),
            "reported product energy residual is small");

    double spill_number = 0.0;
    for (int species = 0; species < 7; ++species)
        spill_number += result.below_number_m3_s[species] +
                        result.above_number_m3_s[species];
    require(spill_number > 0.0,
            "clipped grid exercises an explicit spill ledger");

    std::cout << "channel " << channel << " T/keV " << temperature_keV
              << " rate=" << result.reactivity_m3_s
              << " rate/reference="
              << result.reactivity_m3_s / result.reference_reactivity_m3_s
              << " debit discrepancy="
              << result.relative_reactant_energy_discrepancy << '\n';
}

void check_reference_parity() {
    // Equal-temperature source evaluations have a separately integrated
    // relative-rate/reference path.  A large but finite support leaves only
    // negligible physical tails here; compare both reactant energy debits to
    // those returned independent moments as well as comparing the rate.
    for (int channel = FUSION_DD_TP; channel < FUSION_CHANNEL_COUNT; ++channel) {
        for (double temperature : {0.1, 1.0, 20.0}) {
            const Evaluation evaluation = evaluate(channel, temperature, 8, 40.0);
            check_conservation(evaluation, channel, temperature);
            const auto &result = evaluation.result;
            require(close_scaled(result.reactivity_m3_s,
                                 result.reference_reactivity_m3_s,
                                 3e-7, 1e-300),
                    "non-pB rate agrees with returned independent reference");
            for (int reactant = 0; reactant < 2; ++reactant)
                require(close_scaled(
                            result.reactant_energy_moment_J_m3_s[reactant],
                            result.reference_reactant_energy_moment_J_m3_s[reactant],
                            3e-6, 1e-300),
                        "non-pB selected-reactant debit agrees with reference");
        }
    }
    for (double temperature : {1.0, 20.0}) {
        const Evaluation evaluation = evaluate(FUSION_PB11_3ALPHA, temperature,
                                                8, 40.0);
        check_conservation(evaluation, FUSION_PB11_3ALPHA, temperature);
        const auto &result = evaluation.result;
        require(close_scaled(result.reactivity_m3_s,
                             result.reference_reactivity_m3_s,
                             3e-7, 1e-300),
                "pB rate agrees with returned independent reference");
        for (int reactant = 0; reactant < 2; ++reactant)
            require(close_scaled(result.reactant_energy_moment_J_m3_s[reactant],
                                 result.reference_reactant_energy_moment_J_m3_s[reactant],
                                 3e-6, 1e-300),
                    "pB selected-reactant debit agrees with reference");
    }
}

void check_relative_refinement() {
    // Hold nq/ncos and CM order fixed.  This isolates relative-energy
    // quadrature convergence from the angular alpha-spectrum quadrature.
    for (int channel : {FUSION_DD_TP, FUSION_DD_HE3N,
                        FUSION_DT_ALPHAN, FUSION_DHE3_ALPHAP,
                        FUSION_PB11_3ALPHA}) {
        const Evaluation coarse = evaluate(channel, 20.0, 8, 40.0);
        const Evaluation refined = evaluate(channel, 20.0, 16, 40.0);
        require(close_scaled(coarse.result.reactivity_m3_s,
                             refined.result.reactivity_m3_s,
                             2e-6, 1e-300),
                "relative quadrature rate refinement converges");
        for (int reactant = 0; reactant < 2; ++reactant)
            require(close_scaled(
                        coarse.result.reactant_energy_moment_J_m3_s[reactant],
                        refined.result.reactant_energy_moment_J_m3_s[reactant],
                        2e-5, 1e-300),
                    "relative quadrature debit refinement converges");
        require(coarse.result.cm_retained_probability ==
                    refined.result.cm_retained_probability &&
                    coarse.result.cm_tail_probability ==
                    refined.result.cm_tail_probability,
                "refinement keeps CM cutoff diagnostics fixed");
    }
}

void check_finite_cutoff_diagnostic() {
    // The CM cutoff is deliberately made small for this diagnostic.  Its
    // omitted Maxwellian tail must be visible in the returned probability and
    // in the rate discrepancy; this is independent of relative quadrature.
    const Evaluation truncated = evaluate(FUSION_DT_ALPHAN, 20.0, 8, 8.0);
    const Evaluation nearly_full = evaluate(FUSION_DT_ALPHAN, 20.0, 8, 40.0);
    check_conservation(truncated, FUSION_DT_ALPHAN, 20.0);
    check_conservation(nearly_full, FUSION_DT_ALPHAN, 20.0);
    require(truncated.result.cm_tail_probability > 1e-4 &&
                truncated.result.cm_tail_probability < 0.01,
            "finite CM cutoff reports its non-negligible probability tail");
    require(truncated.result.relative_rate_discrepancy < -5e-4 &&
                truncated.result.reactivity_m3_s <
                    nearly_full.result.reactivity_m3_s * (1.0 - 5e-4),
            "finite CM cutoff lowers the rate and reports the discrepancy");
}

void check_bad_inputs_clear() {
    const std::vector<double> edges = clipped_edges();
    const int cells = static_cast<int>(edges.size() - 1);
    fusion_thermal_birth_options_v1 options = options_for(5.0, 8, 40.0);
    std::vector<double> birth(7 * static_cast<std::size_t>(cells), 7.0);
    fusion_thermal_birth_v1 out{};
    const auto poison = [&]() {
        for (double &value : birth)
            value = 7.0;
        out.reactivity_m3_s = 7.0;
        out.reference_reactivity_m3_s = 7.0;
        std::fill(std::begin(out.reactant_energy_moment_J_m3_s),
                  std::end(out.reactant_energy_moment_J_m3_s), 7.0);
        std::fill(std::begin(out.reference_reactant_energy_moment_J_m3_s),
                  std::end(out.reference_reactant_energy_moment_J_m3_s), 7.0);
        out.product_energy_moment_J_m3_s = 7.0;
        out.number_residual_m3_s = 7.0;
        out.energy_residual_J_m3_s = 7.0;
        out.relative_rate_discrepancy = 7.0;
        out.relative_reactant_energy_discrepancy = 7.0;
        out.cm_retained_probability = 7.0;
        out.cm_tail_probability = 7.0;
        out.cm_tail_energy_moment_J = 7.0;
        out.max_cm_energy_shift_fraction = 7.0;
        out.max_shell_remap_fraction = 7.0;
        std::fill(std::begin(out.below_number_m3_s),
                  std::end(out.below_number_m3_s), 7.0);
        std::fill(std::begin(out.below_energy_J_m3_s),
                  std::end(out.below_energy_J_m3_s), 7.0);
        std::fill(std::begin(out.above_number_m3_s),
                  std::end(out.above_number_m3_s), 7.0);
        std::fill(std::begin(out.above_energy_J_m3_s),
                  std::end(out.above_energy_J_m3_s), 7.0);
    };
    const auto cleared = [&]() {
        require(out.reactivity_m3_s == 0.0 &&
                    out.reference_reactivity_m3_s == 0.0 &&
                    std::all_of(std::begin(out.reactant_energy_moment_J_m3_s),
                                std::end(out.reactant_energy_moment_J_m3_s),
                                [](double value) { return value == 0.0; }) &&
                    std::all_of(std::begin(out.reference_reactant_energy_moment_J_m3_s),
                                std::end(out.reference_reactant_energy_moment_J_m3_s),
                                [](double value) { return value == 0.0; }) &&
                    out.product_energy_moment_J_m3_s == 0.0 &&
                    out.number_residual_m3_s == 0.0 &&
                    out.energy_residual_J_m3_s == 0.0 &&
                    out.relative_rate_discrepancy == 0.0 &&
                    out.relative_reactant_energy_discrepancy == 0.0 &&
                    out.cm_retained_probability == 0.0 &&
                    out.cm_tail_probability == 0.0,
                "bad input clears result diagnostics");
        require(out.cm_tail_energy_moment_J == 0.0 &&
                    out.max_cm_energy_shift_fraction == 0.0 &&
                    out.max_shell_remap_fraction == 0.0 &&
                    std::all_of(std::begin(out.below_number_m3_s),
                                std::end(out.below_number_m3_s),
                                [](double value) { return value == 0.0; }) &&
                    std::all_of(std::begin(out.below_energy_J_m3_s),
                                std::end(out.below_energy_J_m3_s),
                                [](double value) { return value == 0.0; }) &&
                    std::all_of(std::begin(out.above_number_m3_s),
                                std::end(out.above_number_m3_s),
                                [](double value) { return value == 0.0; }) &&
                    std::all_of(std::begin(out.above_energy_J_m3_s),
                                std::end(out.above_energy_J_m3_s),
                                [](double value) { return value == 0.0; }),
                "bad input clears spill diagnostics");
        for (double value : birth)
            require(value == 0.0, "bad input clears birth coefficients");
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();

    poison();
    require(fusion_c_thermal_birth_grid(99, kev, &options, cells, edges.data(),
                                        birth.data(), &out) ==
                PB11_STATUS_OUT_OF_RANGE,
            "bad channel is rejected");
    cleared();

    poison();
    require(fusion_c_thermal_birth_grid(FUSION_DT_ALPHAN, nan, &options, cells,
                                        edges.data(), birth.data(), &out) ==
                PB11_STATUS_INVALID_ARGUMENT,
            "NaN temperature is rejected");
    cleared();

    poison();
    require(fusion_c_thermal_birth_grid(FUSION_DT_ALPHAN, kev, &options, cells,
                                        nullptr, birth.data(), &out) ==
                PB11_STATUS_INVALID_ARGUMENT,
            "null energy grid is rejected");
    cleared();

    poison();
    std::vector<double> bad_edges = edges;
    bad_edges[3] = bad_edges[2];
    require(fusion_c_thermal_birth_grid(FUSION_DT_ALPHAN, kev, &options, cells,
                                        bad_edges.data(), birth.data(), &out) ==
                PB11_STATUS_INVALID_ARGUMENT,
            "non-increasing energy grid is rejected");
    cleared();

    poison();
    options.nq = 3;
    require(fusion_c_thermal_birth_grid(FUSION_DT_ALPHAN, kev, &options, cells,
                                        edges.data(), birth.data(), &out) ==
                PB11_STATUS_INVALID_ARGUMENT,
            "invalid angular quadrature order is rejected");
    cleared();
}

} // namespace

int main() {
    try {
        check_reference_parity();
        check_finite_cutoff_diagnostic();
        check_relative_refinement();
        check_bad_inputs_clear();
        std::cout << "PASS: thermal-birth five-channel conservation, finite-cutoff diagnostics, "
                     "reference parity, relative refinement and input clearing\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
