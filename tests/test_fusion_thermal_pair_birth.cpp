#include <fusion_thermal_birth.h>
#include <fusion_network.h>
#include <fusion_nuclear_data.h>
#include <fusion_rate_model.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kev = 1.602176634e-16;
constexpr double mev = 1.602176634e-13;
constexpr int cells = 16;
constexpr double rate_tolerance = 2.0e-4;
constexpr double debit_tolerance = 5.0e-4;
constexpr double number_tolerance = 3.0e-10;
constexpr double energy_tolerance = 3.0e-10;
constexpr double spectrum_tolerance = 2.0e-6;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

bool close_scaled(double actual, double expected, double relative,
                  double absolute = 0.0) {
    if (actual == expected)
        return true;
    if (!std::isfinite(actual) || !std::isfinite(expected))
        return false;
    return std::abs(actual - expected) <=
           absolute + relative * std::max(std::abs(actual), std::abs(expected));
}

std::vector<double> source_edges() {
    std::vector<double> edges(cells + 1);
    for (int i = 0; i <= cells; ++i)
        edges[static_cast<std::size_t>(i)] =
            2.0 * mev * static_cast<double>(i) / cells;
    return edges;
}

fusion_thermal_birth_options_v1 options_for(int channel) {
    fusion_thermal_birth_options_v1 options{};
    options.relative_max_J =
        (channel == FUSION_PB11_3ALPHA ? 2.0 : 5.0) * mev;
    options.cm_max_kT = 40.0;
    options.ground_state_q_J = 91.84 * kev;
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
    options.relative_order = 4;
    options.cm_order = 4;
    options.nq = 4;
    options.ncos = 4;
    return options;
}

struct Evaluation {
    std::vector<double> edges;
    std::vector<double> birth;
    fusion_thermal_birth_v1 result{};
};

Evaluation evaluate_pair(int channel, double ta, double tb,
                         int correlation_order = 4) {
    Evaluation evaluation;
    evaluation.edges = source_edges();
    evaluation.birth.assign(7 * cells, 0.0);
    fusion_thermal_birth_options_v1 options = options_for(channel);
    const int status = fusion_c_thermal_pair_birth_grid(
        channel, ta, tb, correlation_order, &options, cells,
        evaluation.edges.data(), evaluation.birth.data(), &evaluation.result);
    require(status == PB11_STATUS_OK,
            "pair birth call failed for channel " + std::to_string(channel) +
                " (status " + std::to_string(status) + ")");
    return evaluation;
}

Evaluation evaluate_legacy(int channel, double temperature) {
    Evaluation evaluation;
    evaluation.edges = source_edges();
    evaluation.birth.assign(7 * cells, 0.0);
    fusion_thermal_birth_options_v1 options = options_for(channel);
    const int status = fusion_c_thermal_birth_grid(
        channel, temperature, &options, cells, evaluation.edges.data(),
        evaluation.birth.data(), &evaluation.result);
    require(status == PB11_STATUS_OK,
            "legacy birth call failed for channel " +
                std::to_string(channel) + " (status " +
                std::to_string(status) + ")");
    return evaluation;
}

void require_finite_result(const fusion_thermal_birth_v1 &result,
                           const std::string &label) {
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
        require(std::isfinite(value), label + " has non-finite diagnostics");
    for (double value : result.reactant_energy_moment_J_m3_s)
        require(std::isfinite(value), label + " has non-finite reactant debit");
    for (double value : result.reference_reactant_energy_moment_J_m3_s)
        require(std::isfinite(value), label +
                                             " has non-finite reference debit");
    for (double value : result.below_number_m3_s)
        require(std::isfinite(value), label + " has non-finite lower spill");
    for (double value : result.below_energy_J_m3_s)
        require(std::isfinite(value), label + " has non-finite lower energy");
    for (double value : result.above_number_m3_s)
        require(std::isfinite(value), label + " has non-finite upper spill");
    for (double value : result.above_energy_J_m3_s)
        require(std::isfinite(value), label + " has non-finite upper energy");
}

void check_equal_temperature_parity() {
    // The pair API's equal-T branch is required to preserve the old five-
    // channel evaluator, including the source coefficient bytes.
    for (int channel = 0; channel < FUSION_CHANNEL_COUNT; ++channel) {
        const double temperature = 3.0 * kev;
        const Evaluation legacy = evaluate_legacy(channel, temperature);
        const Evaluation pair = evaluate_pair(channel, temperature, temperature);
        require_finite_result(legacy.result, "legacy equal-T result");
        require_finite_result(pair.result, "pair equal-T result");
        require(legacy.birth.size() == pair.birth.size(),
                "equal-T parity array size");
        require(std::memcmp(legacy.birth.data(), pair.birth.data(),
                            legacy.birth.size() * sizeof(double)) == 0,
                "equal-T pair source differs byte-for-byte from legacy");
        require(std::memcmp(&legacy.result, &pair.result,
                            sizeof(fusion_thermal_birth_v1)) == 0,
                "equal-T pair diagnostics differ byte-for-byte from legacy");
        for (std::size_t i = 0; i < legacy.birth.size(); ++i)
            require(close_scaled(legacy.birth[i], pair.birth[i], 0.0),
                    "equal-T pair source numeric parity");
        require(pair.result.reactivity_m3_s > 0.0,
                "equal-T pair has positive rate");
    }
}

void check_number_and_energy_closure(const Evaluation &evaluation, int channel,
                                     const std::string &label) {
    const auto &result = evaluation.result;
    require_finite_result(result, label);
    require(result.reactivity_m3_s > 0.0 &&
                result.reference_reactivity_m3_s > 0.0,
            label + " has positive rate and reference");

    fusion_nuclear_channel_v1 reaction{};
    require(fusion_c_nuclear_channel(channel, &reaction) == PB11_STATUS_OK,
            label + " channel metadata");
    std::array<int, 7> multiplicity{};
    for (int product = 0; product < reaction.product_count; ++product)
        ++multiplicity[static_cast<std::size_t>(reaction.product_ids[product])];

    std::array<double, 7> number{};
    std::array<double, 7> energy{};
    for (int species = 0; species < 7; ++species) {
        for (int cell = 0; cell < cells; ++cell) {
            const double value = evaluation.birth[
                static_cast<std::size_t>(species * cells + cell)];
            require(std::isfinite(value) && value >= 0.0,
                    label + " has invalid source coefficient");
            number[static_cast<std::size_t>(species)] += value;
            energy[static_cast<std::size_t>(species)] +=
                0.5 * (evaluation.edges[static_cast<std::size_t>(cell)] +
                       evaluation.edges[static_cast<std::size_t>(cell + 1)]) *
                value;
        }
        number[static_cast<std::size_t>(species)] +=
            result.below_number_m3_s[species] + result.above_number_m3_s[species];
        energy[static_cast<std::size_t>(species)] +=
            result.below_energy_J_m3_s[species] +
            result.above_energy_J_m3_s[species];
        require(close_scaled(
                    number[static_cast<std::size_t>(species)],
                    multiplicity[static_cast<std::size_t>(species)] *
                        result.reactivity_m3_s,
                    number_tolerance, 1.0e-300),
                label + " product number does not close with spills");
    }

    double total_energy = 0.0;
    double total_spill_number = 0.0;
    for (int species = 0; species < 7; ++species) {
        total_energy += energy[static_cast<std::size_t>(species)];
        total_spill_number += result.below_number_m3_s[species] +
                              result.above_number_m3_s[species];
    }
    require(total_spill_number > 0.0,
            label + " does not exercise the explicit spill ledger");

    const double expected_energy =
        reaction.q_J * result.reactivity_m3_s +
        result.reactant_energy_moment_J_m3_s[0] +
        result.reactant_energy_moment_J_m3_s[1];
    require(close_scaled(result.product_energy_moment_J_m3_s, expected_energy,
                         2.0e-9, 1.0e-300),
            label + " Q plus debit does not close in product moment");
    require(close_scaled(total_energy, expected_energy, energy_tolerance,
                         1.0e-300),
            label + " Q plus debit does not close with mapped spills");
    require(std::abs(result.number_residual_m3_s) <=
                2.0e-9 * result.reactivity_m3_s,
            label + " number residual is too large");
    require(std::abs(result.energy_residual_J_m3_s) <=
                2.0e-9 * expected_energy,
            label + " energy residual is too large");
}

void check_unequal_reference(int channel, double ta, double tb,
                             const std::string &label) {
    const Evaluation evaluation = evaluate_pair(channel, ta, tb);
    check_number_and_energy_closure(evaluation, channel, label);
    const fusion_thermal_birth_options_v1 options = options_for(channel);

    fusion_nuclear_channel_v1 reaction{};
    require(fusion_c_nuclear_channel(channel, &reaction) == PB11_STATUS_OK,
            label + " channel metadata for reference");
    fusion_nuclear_mass_v1 mass_a{}, mass_b{};
    require(fusion_c_nuclear_mass(reaction.reactant_ids[0], &mass_a) ==
                PB11_STATUS_OK,
            label + " reactant-a mass");
    require(fusion_c_nuclear_mass(reaction.reactant_ids[1], &mass_b) ==
                PB11_STATUS_OK,
            label + " reactant-b mass");

    fusion_rate_model_v1 reference{};
    require(fusion_c_thermal_pair_maxwellian_model(
                channel, options.continuation, options.pb_low, mass_a.mass_kg,
                mass_b.mass_kg, ta, tb, &reference) == PB11_STATUS_OK,
            label + " unequal-temperature reference model");

    const auto &result = evaluation.result;
    require(close_scaled(result.reference_reactivity_m3_s,
                         reference.total.resolved_reactivity_m3_s, 0.0),
            label + " returned reference rate differs from rate model");
    for (int reactant = 0; reactant < 2; ++reactant) {
        require(close_scaled(
                    result.reference_reactant_energy_moment_J_m3_s[reactant],
                    reactant == 0
                        ? reference.total.projectile_energy_reactivity_J_m3_s
                        : reference.total.target_energy_reactivity_J_m3_s,
                    0.0),
                label + " returned reference debit differs from rate model");
        require(close_scaled(result.reactant_energy_moment_J_m3_s[reactant],
                             result.reference_reactant_energy_moment_J_m3_s[
                                 reactant],
                             debit_tolerance, 1.0e-300),
                label + " reacting debit differs from unequal-T reference");
    }
    require(close_scaled(result.reactivity_m3_s,
                         reference.total.resolved_reactivity_m3_s,
                         rate_tolerance, 1.0e-300),
            label + " rate differs from unequal-T reference");

    const double expected_rate_discrepancy =
        result.reactivity_m3_s / result.reference_reactivity_m3_s - 1.0;
    require(close_scaled(result.relative_rate_discrepancy,
                         expected_rate_discrepancy, 1.0e-12, 1.0e-14),
            label + " reports the wrong relative rate discrepancy");
    double expected_debit_discrepancy = 0.0;
    for (int reactant = 0; reactant < 2; ++reactant)
        expected_debit_discrepancy = std::max(
            expected_debit_discrepancy,
            std::abs(result.reactant_energy_moment_J_m3_s[reactant] /
                         result.reference_reactant_energy_moment_J_m3_s[
                             reactant] -
                     1.0));
    require(close_scaled(result.relative_reactant_energy_discrepancy,
                         expected_debit_discrepancy, 1.0e-12, 1.0e-14),
            label + " reports the wrong debit discrepancy");
    require(std::abs(result.relative_rate_discrepancy) <= rate_tolerance &&
                result.relative_reactant_energy_discrepancy <= debit_tolerance,
            label + " quadrature discrepancy exceeds bounded tolerance");
}

void check_dd_swap() {
    const int channel = FUSION_DD_TP;
    const Evaluation forward = evaluate_pair(channel, 5.0 * kev, 12.0 * kev);
    const Evaluation swapped = evaluate_pair(channel, 12.0 * kev, 5.0 * kev);
    require(close_scaled(forward.result.reactivity_m3_s,
                         swapped.result.reactivity_m3_s, 2.0e-8, 1.0e-300),
            "DD temperature swap changes the rate");
    require(close_scaled(forward.result.reference_reactivity_m3_s,
                         swapped.result.reference_reactivity_m3_s, 2.0e-8,
                         1.0e-300),
            "DD temperature swap changes the reference rate");
    require(close_scaled(forward.result.reactant_energy_moment_J_m3_s[0],
                         swapped.result.reactant_energy_moment_J_m3_s[1],
                         2.0e-7, 1.0e-300),
            "DD temperature swap does not exchange debit zero with one");
    require(close_scaled(forward.result.reactant_energy_moment_J_m3_s[1],
                         swapped.result.reactant_energy_moment_J_m3_s[0],
                         2.0e-7, 1.0e-300),
            "DD temperature swap does not exchange debit one with zero");
    require(forward.birth.size() == swapped.birth.size(),
            "DD temperature swap source size");
    for (std::size_t i = 0; i < forward.birth.size(); ++i)
        require(close_scaled(forward.birth[i], swapped.birth[i],
                             spectrum_tolerance, 1.0e-300),
                "DD temperature swap changes the product spectrum");
    for (int species = 0; species < 7; ++species) {
        require(close_scaled(forward.result.below_number_m3_s[species],
                             swapped.result.below_number_m3_s[species],
                             spectrum_tolerance, 1.0e-300),
                "DD temperature swap changes lower spill spectrum");
        require(close_scaled(forward.result.above_number_m3_s[species],
                             swapped.result.above_number_m3_s[species],
                             spectrum_tolerance, 1.0e-300),
                "DD temperature swap changes upper spill spectrum");
    }
}

void check_near_equal_temperature_continuity() {
    // The unequal branch should approach the legacy/equal-T spectrum smoothly
    // as its two pool temperatures merge.  Use a relative 1e-6 perturbation
    // and compare a normalized source L1 difference.
    const double temperature = 20.0 * kev;
    const Evaluation equal =
        evaluate_pair(FUSION_DT_ALPHAN, temperature, temperature);
    const Evaluation near_equal = evaluate_pair(
        FUSION_DT_ALPHAN, temperature, temperature * (1.0 + 1.0e-6));
    double difference = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < equal.birth.size(); ++i) {
        difference += std::abs(equal.birth[i] - near_equal.birth[i]);
        scale += std::abs(equal.birth[i]);
    }
    for (int i=0;i<7;++i) {
        difference += std::abs(equal.result.below_number_m3_s[i]-near_equal.result.below_number_m3_s[i]);
        difference += std::abs(equal.result.above_number_m3_s[i]-near_equal.result.above_number_m3_s[i]);
        scale += equal.result.below_number_m3_s[i]+equal.result.above_number_m3_s[i];
    }
    require(scale > 0.0, "near-equal DT continuity has no source");
    require(difference <= 1.0e-5 * scale,
            "near-equal DT temperatures cause a discontinuous spectrum");
}

void poison_result(fusion_thermal_birth_v1 &result) {
    result.reactivity_m3_s = 7.0;
    result.reference_reactivity_m3_s = 7.0;
    for (double &value : result.reactant_energy_moment_J_m3_s)
        value = 7.0;
    for (double &value : result.reference_reactant_energy_moment_J_m3_s)
        value = 7.0;
    result.product_energy_moment_J_m3_s = 7.0;
    result.number_residual_m3_s = 7.0;
    result.energy_residual_J_m3_s = 7.0;
    result.relative_rate_discrepancy = 7.0;
    result.relative_reactant_energy_discrepancy = 7.0;
    result.cm_retained_probability = 7.0;
    result.cm_tail_probability = 7.0;
    result.cm_tail_energy_moment_J = 7.0;
    result.max_cm_energy_shift_fraction = 7.0;
    result.max_shell_remap_fraction = 7.0;
    for (double &value : result.below_number_m3_s)
        value = 7.0;
    for (double &value : result.below_energy_J_m3_s)
        value = 7.0;
    for (double &value : result.above_number_m3_s)
        value = 7.0;
    for (double &value : result.above_energy_J_m3_s)
        value = 7.0;
}

bool result_is_zero(const fusion_thermal_birth_v1 &result) {
    if (result.reactivity_m3_s != 0.0 ||
        result.reference_reactivity_m3_s != 0.0 ||
        result.product_energy_moment_J_m3_s != 0.0 ||
        result.number_residual_m3_s != 0.0 ||
        result.energy_residual_J_m3_s != 0.0 ||
        result.relative_rate_discrepancy != 0.0 ||
        result.relative_reactant_energy_discrepancy != 0.0 ||
        result.cm_retained_probability != 0.0 ||
        result.cm_tail_probability != 0.0 ||
        result.cm_tail_energy_moment_J != 0.0 ||
        result.max_cm_energy_shift_fraction != 0.0 ||
        result.max_shell_remap_fraction != 0.0)
        return false;
    for (double value : result.reactant_energy_moment_J_m3_s)
        if (value != 0.0) return false;
    for (double value : result.reference_reactant_energy_moment_J_m3_s)
        if (value != 0.0) return false;
    for (double value : result.below_number_m3_s)
        if (value != 0.0) return false;
    for (double value : result.below_energy_J_m3_s)
        if (value != 0.0) return false;
    for (double value : result.above_number_m3_s)
        if (value != 0.0) return false;
    for (double value : result.above_energy_J_m3_s)
        if (value != 0.0) return false;
    return true;
}

void require_cleared(int status, int expected_status,
                     const std::vector<double> &birth,
                     const fusion_thermal_birth_v1 &result,
                     const std::string &label) {
    require(status == expected_status,
            label + " returned status " + std::to_string(status));
    require(result_is_zero(result), label + " did not clear diagnostics");
    for (double value : birth)
        require(value == 0.0, label + " did not clear source coefficients");
}

void check_bad_inputs_clear() {
    const std::vector<double> edges = source_edges();
    fusion_thermal_birth_options_v1 options = options_for(FUSION_DT_ALPHAN);
    const double valid_temperature = 5.0 * kev;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    const auto check_temperature = [&](double ta, double tb, int expected,
                                       const std::string &label) {
        std::vector<double> birth(7 * cells, 7.0);
        fusion_thermal_birth_v1 result{};
        poison_result(result);
        const int status = fusion_c_thermal_pair_birth_grid(
            FUSION_DT_ALPHAN, ta, tb, 4, &options, cells, edges.data(),
            birth.data(), &result);
        require_cleared(status, expected, birth, result, label);
    };

    check_temperature(-valid_temperature, valid_temperature,
                      PB11_STATUS_OUT_OF_RANGE, "negative Ta");
    check_temperature(valid_temperature, -valid_temperature,
                      PB11_STATUS_OUT_OF_RANGE, "negative Tb");
    check_temperature(0.0, valid_temperature, PB11_STATUS_OUT_OF_RANGE,
                      "zero Ta");
    check_temperature(valid_temperature, 0.0, PB11_STATUS_OUT_OF_RANGE,
                      "zero Tb");
    check_temperature(nan, valid_temperature, PB11_STATUS_INVALID_ARGUMENT,
                      "NaN Ta");
    check_temperature(valid_temperature, nan, PB11_STATUS_INVALID_ARGUMENT,
                      "NaN Tb");

    for (int order : {3, 33}) {
        std::vector<double> birth(7 * cells, 7.0);
        fusion_thermal_birth_v1 result{};
        poison_result(result);
        const int status = fusion_c_thermal_pair_birth_grid(
            FUSION_DT_ALPHAN, valid_temperature, valid_temperature, order,
            &options, cells, edges.data(), birth.data(), &result);
        require_cleared(status, PB11_STATUS_INVALID_ARGUMENT, birth, result,
                        "invalid polar order " + std::to_string(order));
    }

    {
        fusion_thermal_birth_v1 result{};
        poison_result(result);
        const int status = fusion_c_thermal_pair_birth_grid(
            FUSION_DT_ALPHAN, valid_temperature, valid_temperature, 4,
            &options, cells, edges.data(), nullptr, &result);
        require(status == PB11_STATUS_NULL_OUTPUT,
                "null birth pointer status");
        require(result_is_zero(result), "null birth pointer clears result");
    }
    {
        std::vector<double> birth(7 * cells, 7.0);
        const int status = fusion_c_thermal_pair_birth_grid(
            FUSION_DT_ALPHAN, valid_temperature, valid_temperature, 4,
            &options, cells, edges.data(), birth.data(), nullptr);
        require(status == PB11_STATUS_NULL_OUTPUT,
                "null result pointer status");
        for (double value : birth)
            require(value == 0.0, "null result pointer clears birth");
    }
    {
        std::vector<double> birth(7 * cells, 7.0);
        fusion_thermal_birth_v1 result{};
        poison_result(result);
        const int status = fusion_c_thermal_pair_birth_grid(
            FUSION_DT_ALPHAN, valid_temperature, valid_temperature, 4,
            nullptr, cells, edges.data(), birth.data(), &result);
        require_cleared(status, PB11_STATUS_INVALID_ARGUMENT, birth, result,
                        "null options pointer");
    }
    {
        std::vector<double> birth(7 * cells, 7.0);
        fusion_thermal_birth_v1 result{};
        poison_result(result);
        const int status = fusion_c_thermal_pair_birth_grid(
            FUSION_DT_ALPHAN, valid_temperature, valid_temperature, 4,
            &options, cells, nullptr, birth.data(), &result);
        require_cleared(status, PB11_STATUS_INVALID_ARGUMENT, birth, result,
                        "null edges pointer");
    }
}

} // namespace

int main() {
    try {
        check_equal_temperature_parity();
        check_unequal_reference(FUSION_DT_ALPHAN, 5.0 * kev, 12.0 * kev,
                                "unequal DT");
        check_unequal_reference(FUSION_DD_TP, 5.0 * kev, 12.0 * kev,
                                "unequal DD-Tp");
        check_unequal_reference(FUSION_DD_HE3N, 5.0 * kev, 12.0 * kev,
                                "unequal DD-He3n");
        check_unequal_reference(FUSION_DHE3_ALPHAP, 5.0 * kev, 12.0 * kev,
                                "unequal DHe3");
        check_dd_swap();
        check_near_equal_temperature_continuity();
        check_bad_inputs_clear();
        std::cout << "PASS: pair thermal birth equal-T parity, unequal-T reference "
                     "debits, conservation, DD swap, and clearing\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
