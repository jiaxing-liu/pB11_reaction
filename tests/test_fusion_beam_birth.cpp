#include "fusion_beam_birth.h"
#include "fusion_network.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double MeV = 1.602176634e-13;
constexpr double keV = 1.602176634e-16;
constexpr int cells = 16;
constexpr double conservation_tolerance = 3.0e-10;
constexpr double reference_tolerance = 1.0e-4;

void require(bool condition, const std::string &message) {
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

fusion_beam_birth_options_v1 options_for(double relative_max_MeV) {
    fusion_beam_birth_options_v1 options{};
    options.relative_max_J = relative_max_MeV * MeV;
    options.angular_max_exponent = 40.0;
    options.ground_state_q_J = 91.84e-3 * MeV;
    // Keep the exact SI spelling used by the alpha-event lower-bound check.
    options.cutoff_J = 0.001 * MeV;
    options.l1_fraction = 0.76;
    options.relative_phase = 0.0;
    options.narrow_peak_fraction = 0.051;
    options.continuum_peak_scale = 1.0;
    options.continuation = FUSION_ENDPOINT_S;
    options.pb_low = FUSION_PB_LOW_TB;
    options.remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY;
    options.broad_mode = 13;
    options.fsci_policy = 0;
    options.relative_order = 16;
    options.angular_order = 8;
    options.nq = 4;
    options.ncos = 4;
    return options;
}

std::vector<double> source_edges() {
    // The four-MeV center grid retains the ordinary charged products while
    // forcing the DT neutron, and any high alpha tail, into explicit spills.
    std::vector<double> edges(cells + 1);
    for (int i = 0; i <= cells; ++i)
        edges[static_cast<std::size_t>(i)] = 4.0 * MeV * i / cells;
    return edges;
}

struct ChannelData {
    fusion_nuclear_channel_v1 channel{};
    std::array<fusion_nuclear_mass_v1, 2> reactant{};
};

ChannelData channel_data(int channel) {
    ChannelData data{};
    require(fusion_c_nuclear_channel(channel, &data.channel) == PB11_STATUS_OK,
            "channel metadata");
    for (int slot = 0; slot < 2; ++slot)
        require(fusion_c_nuclear_mass(data.channel.reactant_ids[slot],
                                      &data.reactant[slot]) == PB11_STATUS_OK,
                "reactant metadata");
    return data;
}

double exact_cold_rate(int channel, int projectile_slot,
                       double projectile_energy_J,
                       const fusion_beam_birth_options_v1 &options,
                       double *relative_energy_J = nullptr) {
    const ChannelData data = channel_data(channel);
    const double projectile_mass =
        data.reactant[static_cast<std::size_t>(projectile_slot)].mass_kg;
    const double target_mass =
        data.reactant[static_cast<std::size_t>(1 - projectile_slot)].mass_kg;
    const double reduced_mass = projectile_mass * target_mass /
                                (projectile_mass + target_mass);
    const double relative_energy = reduced_mass / projectile_mass *
                                   projectile_energy_J;
    if (relative_energy_J)
        *relative_energy_J = relative_energy;
    double sigma_m2 = 0.0;
    require(fusion_c_cross_section_model(
                channel, options.continuation, options.pb_low, relative_energy,
                &sigma_m2) == PB11_STATUS_OK,
            "cold cross section model");
    require(std::isfinite(sigma_m2) && sigma_m2 >= 0.0,
            "cold cross section is finite");
    const double speed = std::sqrt(2.0 * projectile_energy_J / projectile_mass);
    return sigma_m2 * speed;
}

void require_finite_nonnegative(const fusion_thermal_birth_v1 &spectrum) {
    const double scalar[] = {
        spectrum.reactivity_m3_s,
        spectrum.reference_reactivity_m3_s,
        spectrum.product_energy_moment_J_m3_s,
        spectrum.number_residual_m3_s,
        spectrum.energy_residual_J_m3_s,
        spectrum.relative_rate_discrepancy,
        spectrum.relative_reactant_energy_discrepancy,
        spectrum.cm_retained_probability,
        spectrum.cm_tail_probability,
        spectrum.cm_tail_energy_moment_J,
        spectrum.max_cm_energy_shift_fraction,
        spectrum.max_shell_remap_fraction,
    };
    for (double value : scalar)
        require(std::isfinite(value), "spectrum diagnostic is finite");
    for (double value : spectrum.reactant_energy_moment_J_m3_s)
        require(std::isfinite(value), "reactant moment is finite");
    for (double value : spectrum.reference_reactant_energy_moment_J_m3_s)
        require(std::isfinite(value), "reference reactant moment is finite");
    for (double value : spectrum.below_number_m3_s)
        require(std::isfinite(value) && value >= 0.0,
                "below number is finite and nonnegative");
    for (double value : spectrum.below_energy_J_m3_s)
        require(std::isfinite(value) && value >= 0.0,
                "below energy is finite and nonnegative");
    for (double value : spectrum.above_number_m3_s)
        require(std::isfinite(value) && value >= 0.0,
                "above number is finite and nonnegative");
    for (double value : spectrum.above_energy_J_m3_s)
        require(std::isfinite(value) && value >= 0.0,
                "above energy is finite and nonnegative");
}

struct Evaluation {
    int channel = 0;
    int projectile_slot = 0;
    double projectile_energy_J = 0.0;
    double target_kT_J = 0.0;
    std::vector<double> edges;
    std::vector<double> birth;
    fusion_beam_birth_v1 result{};
};

Evaluation evaluate(int channel, int projectile_slot, double projectile_energy_J,
                    double target_kT_J,
                    const fusion_beam_birth_options_v1 &options) {
    Evaluation evaluation;
    evaluation.channel = channel;
    evaluation.projectile_slot = projectile_slot;
    evaluation.projectile_energy_J = projectile_energy_J;
    evaluation.target_kT_J = target_kT_J;
    evaluation.edges = source_edges();
    evaluation.birth.assign(7 * static_cast<std::size_t>(cells), 7.0);
    const int status = fusion_c_beam_birth_grid(
        channel, projectile_slot, projectile_energy_J, target_kT_J, &options,
        cells, evaluation.edges.data(), evaluation.birth.data(),
        &evaluation.result);
    require(status == PB11_STATUS_OK, "beam birth evaluation status");
    require_finite_nonnegative(evaluation.result.spectrum);
    require(std::isfinite(evaluation.result.relative_retained_probability) &&
                std::isfinite(evaluation.result.retained_pair_probability) &&
                std::isfinite(evaluation.result.angular_omitted_pair_probability),
            "beam birth probabilities are finite");
    return evaluation;
}

void check_number_and_energy_closure(const Evaluation &evaluation) {
    const fusion_thermal_birth_v1 &spectrum = evaluation.result.spectrum;
    const fusion_nuclear_channel_v1 &reaction =
        channel_data(evaluation.channel).channel;
    std::array<double, 7> number{};
    std::array<double, 7> energy{};
    for (int species = 0; species < 7; ++species) {
        for (int cell = 0; cell < cells; ++cell) {
            const std::size_t index = static_cast<std::size_t>(species) * cells +
                                      static_cast<std::size_t>(cell);
            const double value = evaluation.birth[index];
            require(std::isfinite(value) && value >= 0.0,
                    "birth coefficient is finite and nonnegative");
            number[static_cast<std::size_t>(species)] += value;
            energy[static_cast<std::size_t>(species)] +=
                0.5 * (evaluation.edges[static_cast<std::size_t>(cell)] +
                       evaluation.edges[static_cast<std::size_t>(cell + 1)]) *
                value;
        }
        number[static_cast<std::size_t>(species)] +=
            spectrum.below_number_m3_s[species] +
            spectrum.above_number_m3_s[species];
        energy[static_cast<std::size_t>(species)] +=
            spectrum.below_energy_J_m3_s[species] +
            spectrum.above_energy_J_m3_s[species];
    }

    std::array<int, 7> multiplicity{};
    for (int product = 0; product < reaction.product_count; ++product) {
        require(reaction.product_ids[product] >= 0 &&
                    reaction.product_ids[product] < 7,
                "product ID fits explicit seven-slot output");
        ++multiplicity[static_cast<std::size_t>(reaction.product_ids[product])];
    }

    double total_number = 0.0;
    double total_energy = 0.0;
    double total_spill = 0.0;
    for (int species = 0; species < 7; ++species) {
        const double expected_number =
            multiplicity[static_cast<std::size_t>(species)] *
            spectrum.reactivity_m3_s;
        require(close_scaled(number[static_cast<std::size_t>(species)],
                             expected_number, conservation_tolerance, 1.0e-300),
                "product number closes including spills");
        total_number += number[static_cast<std::size_t>(species)];
        total_energy += energy[static_cast<std::size_t>(species)];
        total_spill += spectrum.below_number_m3_s[species] +
                       spectrum.above_number_m3_s[species];
    }
    require(close_scaled(total_number,
                         reaction.product_count * spectrum.reactivity_m3_s,
                         conservation_tolerance, 1.0e-300),
            "total product number closes");
    require(close_scaled(total_energy,
                         spectrum.product_energy_moment_J_m3_s,
                         conservation_tolerance, 1.0e-300),
            "mapped product energy closes");

    const double expected_energy =
        reaction.q_J * spectrum.reactivity_m3_s +
        spectrum.reactant_energy_moment_J_m3_s[0] +
        spectrum.reactant_energy_moment_J_m3_s[1];
    require(close_scaled(total_energy, expected_energy,
                         conservation_tolerance, 1.0e-300),
            "Q plus debited reactant energy closes");
    require(close_scaled(spectrum.number_residual_m3_s, 0.0, 0.0,
                         conservation_tolerance *
                             std::max(spectrum.reactivity_m3_s, 1.0e-300)),
            "reported product number residual closes");
    require(close_scaled(spectrum.energy_residual_J_m3_s, 0.0, 0.0,
                         conservation_tolerance *
                             std::max(expected_energy, 1.0e-300)),
            "reported energy residual closes");
    require(total_spill > 0.0, "explicit source spill is populated");

    if (evaluation.channel == FUSION_DT_ALPHAN) {
        require(multiplicity[FUSION_MASS_NEUTRON] == 1,
                "DT has one explicit neutron product");
        require(number[FUSION_MASS_NEUTRON] > 0.0,
                "explicit neutron number is retained");
        require(spectrum.above_number_m3_s[FUSION_MASS_NEUTRON] > 0.0,
                "explicit neutron is present in the spill ledger");
    }
}

void check_cold_rate_and_canonical_debits(int channel, int projectile_slot,
                                          double projectile_energy_J,
                                          const fusion_beam_birth_options_v1 &options) {
    double relative_energy_J = 0.0;
    const double exact = exact_cold_rate(channel, projectile_slot,
                                         projectile_energy_J, options,
                                         &relative_energy_J);
    require(exact > 0.0, "cold exact sigma*v is positive");
    const Evaluation evaluation =
        evaluate(channel, projectile_slot, projectile_energy_J, 0.0, options);
    const fusion_thermal_birth_v1 &spectrum = evaluation.result.spectrum;
    require(spectrum.reactivity_m3_s > 0.0,
            "cold beam birth rate is positive");
    require(close_scaled(spectrum.reactivity_m3_s, exact, 3.0e-12,
                         1.0e-300),
            "cold beam birth rate equals exact sigma*v");
    require(close_scaled(spectrum.reference_reactivity_m3_s, exact,
                         3.0e-12, 1.0e-300),
            "cold reference rate equals exact sigma*v");

    const ChannelData data = channel_data(channel);
    // Both debit arrays are indexed by canonical channel reactant slot.  The
    // projectile_slot argument selects which canonical reactant carries the
    // beam energy; it does not reorder either returned array.
    const double expected_debit[2] = {
        projectile_slot == 0 ? projectile_energy_J * exact : 0.0,
        projectile_slot == 1 ? projectile_energy_J * exact : 0.0,
    };
    (void)data;
    for (int slot = 0; slot < 2; ++slot) {
        require(close_scaled(spectrum.reactant_energy_moment_J_m3_s[slot],
                             expected_debit[slot], 3.0e-12,
                             1.0e-300),
                "cold debits use canonical projectile slot order");
        require(close_scaled(
                    spectrum.reference_reactant_energy_moment_J_m3_s[slot],
                    expected_debit[slot], 3.0e-12, 1.0e-300),
                "cold reference debits use canonical slot order");
    }
    require(relative_energy_J <= options.relative_max_J,
            "cold event is inside the finite relative cutoff");
    require(close_scaled(evaluation.result.relative_retained_probability, 1.0,
                         3.0e-12, 1.0e-300),
            "cold in-cutoff event has full relative retention");
    check_number_and_energy_closure(evaluation);
}

void check_warm_dt_against_model() {
    const int channel = FUSION_DT_ALPHAN;
    const int projectile_slot = 0;
    const double projectile_energy = 1.0 * MeV;
    const double target_kT = 20.0 * keV;
    const fusion_beam_birth_options_v1 options = options_for(5.0);
    const Evaluation evaluation = evaluate(channel, projectile_slot,
                                           projectile_energy, target_kT,
                                           options);
    const ChannelData data = channel_data(channel);
    fusion_rate_model_v1 model{};
    require(fusion_c_beam_maxwellian_model(
                channel, options.continuation, options.pb_low,
                data.reactant[0].mass_kg, data.reactant[1].mass_kg,
                projectile_energy, target_kT, &model) == PB11_STATUS_OK,
            "warm DT beam model status");
    const fusion_thermal_birth_v1 &spectrum = evaluation.result.spectrum;
    require(model.total.resolved_reactivity_m3_s > 0.0,
            "warm DT model rate is positive");
    require(close_scaled(spectrum.reactivity_m3_s,
                         model.total.resolved_reactivity_m3_s,
                         reference_tolerance, 1.0e-300),
            "warm DT source rate matches beam model");
    for (int slot = 0; slot < 2; ++slot) {
        const double expected = slot == 0
                                    ? model.total.projectile_energy_reactivity_J_m3_s
                                    : model.total.target_energy_reactivity_J_m3_s;
        require(close_scaled(spectrum.reactant_energy_moment_J_m3_s[slot],
                             expected, reference_tolerance, 1.0e-300),
                "warm DT debits match beam model in canonical order");
    }
    require(spectrum.reactant_energy_moment_J_m3_s[1] > 0.0,
            "warm target debit is positive");
    const double naive_target = 1.5 * target_kT * spectrum.reactivity_m3_s;
    require(!close_scaled(spectrum.reactant_energy_moment_J_m3_s[1],
                          naive_target, reference_tolerance, 1.0e-300),
            "warm target debit is conditional rather than 3T/2");
    check_number_and_energy_closure(evaluation);
}

void check_piecewise_fit_boundaries() {
    auto options = options_for(2.5);
    options.angular_order = 16;
    // Captured BALDUR800-cell zone19 input. Rate accuracy is independent of
    // the outgoing mapper grid, whose complete number/energy budget is checked.
    auto actual = evaluate(FUSION_DT_ALPHAN, 1, 2.1399859254691007e-13,
                           1.4401072590598843e-17, options);
    require(std::abs(actual.result.spectrum.relative_rate_discrepancy) < 1e-8,
            "actual zone19 source resolves DT fit boundary");
    require(actual.result.spectrum.relative_reactant_energy_discrepancy < 1e-8,
            "actual zone19 debit resolves DT fit boundary");
    check_number_and_energy_closure(actual);
    for (int channel : {FUSION_DT_ALPHAN, FUSION_DHE3_ALPHAP}) {
        auto data = channel_data(channel);
        double boundary = (channel == FUSION_DT_ALPHAN ? 530. : 900.) * keV;
        double energy = boundary * (data.reactant[0].mass_kg +
            data.reactant[1].mass_kg) / data.reactant[1].mass_kg;
        for (double factor : {.995, 1., 1.005}) {
            auto value = evaluate(channel, 0, energy * factor, .1 * keV, options);
            require(std::abs(value.result.spectrum.relative_rate_discrepancy) < 1e-8,
                    "piecewise source rate matches independent reference");
            require(value.result.spectrum.relative_reactant_energy_discrepancy < 1e-8,
                    "piecewise source debit matches independent reference");
            check_number_and_energy_closure(value);
        }
    }
}

void check_zero_energy_cold_source() {
    const fusion_beam_birth_options_v1 options = options_for(5.0);
    const Evaluation evaluation =
        evaluate(FUSION_DT_ALPHAN, 0, 0.0, 0.0, options);
    const fusion_thermal_birth_v1 &spectrum = evaluation.result.spectrum;
    require(spectrum.reactivity_m3_s == 0.0 &&
                spectrum.reference_reactivity_m3_s == 0.0,
            "zero-energy cold projectile has zero source rate");
    require(spectrum.reactant_energy_moment_J_m3_s[0] == 0.0 &&
                spectrum.reactant_energy_moment_J_m3_s[1] == 0.0 &&
                spectrum.product_energy_moment_J_m3_s == 0.0,
            "zero-energy cold projectile has zero energy source");
    for (double value : evaluation.birth)
        require(value == 0.0, "zero-energy cold projectile has zero birth grid");
    for (int species = 0; species < 7; ++species)
        require(evaluation.result.spectrum.below_number_m3_s[species] == 0.0 &&
                    evaluation.result.spectrum.above_number_m3_s[species] == 0.0,
                "zero-energy cold projectile has zero spills");
}

void check_cold_cutoff_exclusion() {
    const fusion_beam_birth_options_v1 options = options_for(5.0);
    // For a D projectile, E_rel=(mu/m_D)E_beam, so 9.5 MeV is outside
    // the five-MeV finite source cutoff while the declared model remains
    // defined and positive.  The DT parent remains in the event domain.
    const Evaluation evaluation = evaluate(FUSION_DT_ALPHAN, 0, 9.5 * MeV,
                                           0.0, options);
    const fusion_thermal_birth_v1 &spectrum = evaluation.result.spectrum;
    require(spectrum.reactivity_m3_s == 0.0,
            "cold event outside relative cutoff is excluded");
    require(spectrum.reference_reactivity_m3_s > 0.0,
            "cold event outside cutoff retains positive full reference");
    require(evaluation.result.relative_retained_probability == 0.0,
            "cold event outside cutoff has zero retained probability");
    require(evaluation.result.retained_pair_probability == 0.0,
            "cold event outside cutoff has zero retained pair probability");
    require(evaluation.result.angular_omitted_pair_probability == 0.0,
            "cold excluded event has no angular probability contribution");
    for (double value : evaluation.birth)
        require(value == 0.0, "excluded cold event has zero birth grid");
}

void poison(fusion_beam_birth_v1 &out, std::vector<double> &birth) {
    std::fill(birth.begin(), birth.end(), 7.0);
    fusion_thermal_birth_v1 &s = out.spectrum;
    s.reactivity_m3_s = 7.0;
    s.reference_reactivity_m3_s = 7.0;
    std::fill(std::begin(s.reactant_energy_moment_J_m3_s),
              std::end(s.reactant_energy_moment_J_m3_s), 7.0);
    std::fill(std::begin(s.reference_reactant_energy_moment_J_m3_s),
              std::end(s.reference_reactant_energy_moment_J_m3_s), 7.0);
    s.product_energy_moment_J_m3_s = 7.0;
    s.number_residual_m3_s = 7.0;
    s.energy_residual_J_m3_s = 7.0;
    s.relative_rate_discrepancy = 7.0;
    s.relative_reactant_energy_discrepancy = 7.0;
    s.cm_retained_probability = 7.0;
    s.cm_tail_probability = 7.0;
    s.cm_tail_energy_moment_J = 7.0;
    s.max_cm_energy_shift_fraction = 7.0;
    s.max_shell_remap_fraction = 7.0;
    std::fill(std::begin(s.below_number_m3_s), std::end(s.below_number_m3_s),
              7.0);
    std::fill(std::begin(s.below_energy_J_m3_s), std::end(s.below_energy_J_m3_s),
              7.0);
    std::fill(std::begin(s.above_number_m3_s), std::end(s.above_number_m3_s),
              7.0);
    std::fill(std::begin(s.above_energy_J_m3_s), std::end(s.above_energy_J_m3_s),
              7.0);
    out.relative_retained_probability = 7.0;
    out.retained_pair_probability = 7.0;
    out.angular_omitted_pair_probability = 7.0;
}

void require_cleared(const fusion_beam_birth_v1 &out,
                     const std::vector<double> &birth) {
    const fusion_thermal_birth_v1 &s = out.spectrum;
    const bool scalars_clear =
        s.reactivity_m3_s == 0.0 && s.reference_reactivity_m3_s == 0.0 &&
                s.product_energy_moment_J_m3_s == 0.0 &&
                s.number_residual_m3_s == 0.0 &&
                s.energy_residual_J_m3_s == 0.0 &&
                s.relative_rate_discrepancy == 0.0 &&
                s.relative_reactant_energy_discrepancy == 0.0 &&
                s.cm_retained_probability == 0.0 &&
                s.cm_tail_probability == 0.0 && s.cm_tail_energy_moment_J == 0.0 &&
                s.max_cm_energy_shift_fraction == 0.0 &&
                s.max_shell_remap_fraction == 0.0;
    require(scalars_clear, "scalar failure outputs clear");
    for (double value : s.reactant_energy_moment_J_m3_s)
        require(value == 0.0, "reactant failure outputs clear");
    for (double value : s.reference_reactant_energy_moment_J_m3_s)
        require(value == 0.0, "reference reactant failure outputs clear");
    for (double value : s.below_number_m3_s)
        require(value == 0.0, "below number failure outputs clear");
    for (double value : s.below_energy_J_m3_s)
        require(value == 0.0, "below energy failure outputs clear");
    for (double value : s.above_number_m3_s)
        require(value == 0.0, "above number failure outputs clear");
    for (double value : s.above_energy_J_m3_s)
        require(value == 0.0, "above energy failure outputs clear");
    require(out.relative_retained_probability == 0.0 &&
                out.retained_pair_probability == 0.0 &&
                out.angular_omitted_pair_probability == 0.0,
            "beam probability failure outputs clear");
    for (double value : birth)
        require(value == 0.0, "birth failure outputs clear");
}

template <typename Call>
void expect_failure_and_clear(Call call, const char *label) {
    std::vector<double> birth = source_edges();
    birth.assign(7 * static_cast<std::size_t>(cells), 0.0);
    fusion_beam_birth_v1 out{};
    poison(out, birth);
    const int status = call(birth, out);
    require(status != PB11_STATUS_OK, std::string(label) + " is rejected");
    require_cleared(out, birth);
}

void check_invalid_inputs_clear_all_outputs() {
    const std::vector<double> edges = source_edges();
    const fusion_beam_birth_options_v1 good = options_for(5.0);
    const double projectile_energy = 1.0 * MeV;
    const double target_kT = 20.0 * keV;

    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, projectile_energy, target_kT, nullptr,
                cells, edges.data(), birth.data(), &out);
        },
        "null options");

    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, projectile_energy, target_kT, &good,
                cells, nullptr, birth.data(), &out);
        },
        "null edges");

    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, std::numeric_limits<double>::quiet_NaN(),
                target_kT, &good, cells, edges.data(), birth.data(), &out);
        },
        "NaN projectile energy");

    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, projectile_energy,
                std::numeric_limits<double>::quiet_NaN(), &good, cells,
                edges.data(), birth.data(), &out);
        },
        "NaN target temperature");

    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, -projectile_energy, target_kT, &good,
                cells, edges.data(), birth.data(), &out);
        },
        "negative projectile energy");

    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, projectile_energy, -target_kT, &good,
                cells, edges.data(), birth.data(), &out);
        },
        "negative target temperature");

    for (const char *label : {"relative order", "angular order", "nq",
                              "ncos"}) {
        fusion_beam_birth_options_v1 bad = good;
        if (std::string(label) == "relative order")
            bad.relative_order = 3;
        else if (std::string(label) == "angular order")
            bad.angular_order = 3;
        else if (std::string(label) == "nq")
            bad.nq = 3;
        else
            bad.ncos = 3;
        expect_failure_and_clear(
            [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
                return fusion_c_beam_birth_grid(
                    FUSION_DT_ALPHAN, 0, projectile_energy, target_kT, &bad,
                    cells, edges.data(), birth.data(), &out);
            },
            label);
    }

    fusion_beam_birth_options_v1 nan_options = good;
    nan_options.angular_max_exponent = std::numeric_limits<double>::quiet_NaN();
    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, projectile_energy, target_kT, &nan_options,
                cells, edges.data(), birth.data(), &out);
        },
        "NaN angular cutoff");

    fusion_beam_birth_options_v1 bad_angular = good;
    bad_angular.angular_max_exponent = 7.0;
    expect_failure_and_clear(
        [&](std::vector<double> &birth, fusion_beam_birth_v1 &out) {
            return fusion_c_beam_birth_grid(
                FUSION_DT_ALPHAN, 0, projectile_energy, target_kT,
                &bad_angular, cells, edges.data(), birth.data(), &out);
        },
        "angular cutoff below range");

    {
        std::vector<double> birth(7 * static_cast<std::size_t>(cells), 7.0);
        const int status = fusion_c_beam_birth_grid(
            FUSION_DT_ALPHAN, 0, projectile_energy, target_kT, &good, cells,
            edges.data(), birth.data(), nullptr);
        require(status != PB11_STATUS_OK, "null result is rejected");
        for (double value : birth)
            require(value == 0.0, "null result clears non-null birth output");
    }

    {
        std::vector<double> birth(7 * static_cast<std::size_t>(cells), 7.0);
        fusion_beam_birth_v1 out{};
        poison(out, birth);
        const int status = fusion_c_beam_birth_grid(
            FUSION_DT_ALPHAN, 0, projectile_energy, target_kT, &good, cells,
            edges.data(), nullptr, &out);
        require(status != PB11_STATUS_OK, "null birth is rejected");
        // There is no writable birth array, but every non-null output must
        // still be cleared on this failure path.
        std::vector<double> no_birth;
        require_cleared(out, no_birth);
    }
}

} // namespace

int main() {
    try {
        const fusion_beam_birth_options_v1 dt_options = options_for(5.0);
        const fusion_beam_birth_options_v1 pb_options = options_for(2.5);
        check_cold_rate_and_canonical_debits(FUSION_DT_ALPHAN, 0, 1.0 * MeV,
                                             dt_options);
        check_cold_rate_and_canonical_debits(FUSION_DT_ALPHAN, 1, 1.0 * MeV,
                                             dt_options);
        check_cold_rate_and_canonical_debits(FUSION_PB11_3ALPHA, 0, 1.0 * MeV,
                                             pb_options);
        check_cold_rate_and_canonical_debits(FUSION_PB11_3ALPHA, 1, 5.0 * MeV,
                                             pb_options);
        check_warm_dt_against_model();
        check_piecewise_fit_boundaries();
        check_zero_energy_cold_source();
        check_cold_cutoff_exclusion();
        check_invalid_inputs_clear_all_outputs();
        std::cout << "PASS: bounded beam birth cold/warm rates, canonical debits, "
                     "product and Q-energy closure, cutoff, zero-source, and clearing\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
