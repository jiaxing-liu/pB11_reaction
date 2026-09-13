#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
constexpr double kev = 1.602176634e-16;
constexpr double model_tolerance = 2e-8;

void require(bool condition, const char *label) {
    if (!condition) throw std::runtime_error(label);
}

bool near(double a, double b, double tolerance = model_tolerance) {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    // Deliberately use no absolute floor: a nonzero value must agree
    // relatively even when the reference is a very small tail contribution.
    return std::abs(a - b) <= tolerance * std::max(std::abs(a), std::abs(b));
}

void require_near(double actual, double expected, const char *label,
                  double tolerance = model_tolerance) {
    if (!near(actual, expected, tolerance)) {
        std::cerr << "FAIL: " << label << " actual=" << actual
                  << " expected=" << expected << '\n';
        throw std::runtime_error(label);
    }
}

bool zero_window(const fusion_beam_window_v1 &x) {
    return x.resolved_reactivity_m3_s == 0 &&
           x.projectile_energy_reactivity_J_m3_s == 0 &&
           x.target_energy_reactivity_J_m3_s == 0 &&
           x.relative_energy_reactivity_J_m3_s == 0 &&
           x.cm_energy_reactivity_J_m3_s == 0 &&
           x.resolved_pair_probability == 0 &&
           x.unresolved_pair_probability == 0 &&
           x.unresolved_relative_speed_m_s == 0 &&
           x.quadrature_error_m3_s == 0 &&
           x.energy_identity_error_J_m3_s == 0 &&
           x.domain_incomplete == 0;
}

bool zero_model(const fusion_rate_model_v1 &x) {
    return zero_window(x.total) && zero_window(x.fit) &&
           zero_window(x.below) && zero_window(x.above);
}

void check_window(const fusion_beam_window_v1 &x, const char *label) {
    const double nonnegative[] = {
        x.resolved_reactivity_m3_s,
        x.projectile_energy_reactivity_J_m3_s,
        x.target_energy_reactivity_J_m3_s,
        x.relative_energy_reactivity_J_m3_s,
        x.cm_energy_reactivity_J_m3_s,
        x.resolved_pair_probability,
        x.unresolved_pair_probability,
        x.unresolved_relative_speed_m_s,
        x.quadrature_error_m3_s,
    };
    for (double value : nonnegative)
        require(std::isfinite(value) && value >= 0, label);
    require(std::isfinite(x.energy_identity_error_J_m3_s), label);
    require(x.domain_incomplete == 0 || x.domain_incomplete == 1, label);
}

void check_energy_identity(const fusion_beam_window_v1 &x, const char *label) {
    check_window(x, label);
    require_near(x.projectile_energy_reactivity_J_m3_s +
                     x.target_energy_reactivity_J_m3_s,
                 x.relative_energy_reactivity_J_m3_s +
                     x.cm_energy_reactivity_J_m3_s,
                 label);
}

void check_model(const fusion_rate_model_v1 &x) {
    check_energy_identity(x.total, "model total energy identity");
    check_energy_identity(x.fit, "model fit energy identity");
    check_energy_identity(x.below, "model below energy identity");
    check_energy_identity(x.above, "model above energy identity");
    require(x.total.domain_incomplete == 0, "model total defines all energies");
    require(x.total.resolved_pair_probability >= 0 &&
                x.total.resolved_pair_probability <= 1 + model_tolerance,
            "model total probability");
}

void check_segment_sum(const fusion_rate_model_v1 &x) {
    require_near(x.total.resolved_reactivity_m3_s,
                 x.fit.resolved_reactivity_m3_s +
                     x.below.resolved_reactivity_m3_s +
                     x.above.resolved_reactivity_m3_s,
                 "model total rate is segment sum");
    require_near(x.total.projectile_energy_reactivity_J_m3_s,
                 x.fit.projectile_energy_reactivity_J_m3_s +
                     x.below.projectile_energy_reactivity_J_m3_s +
                     x.above.projectile_energy_reactivity_J_m3_s,
                 "model projectile moment is segment sum");
    require_near(x.total.target_energy_reactivity_J_m3_s,
                 x.fit.target_energy_reactivity_J_m3_s +
                     x.below.target_energy_reactivity_J_m3_s +
                     x.above.target_energy_reactivity_J_m3_s,
                 "model target moment is segment sum");
    require_near(x.total.relative_energy_reactivity_J_m3_s,
                 x.fit.relative_energy_reactivity_J_m3_s +
                     x.below.relative_energy_reactivity_J_m3_s +
                     x.above.relative_energy_reactivity_J_m3_s,
                 "model relative moment is segment sum");
    require_near(x.total.cm_energy_reactivity_J_m3_s,
                 x.fit.cm_energy_reactivity_J_m3_s +
                     x.below.cm_energy_reactivity_J_m3_s +
                     x.above.cm_energy_reactivity_J_m3_s,
                 "model CM moment is segment sum");
    require_near(x.total.resolved_pair_probability,
                 x.fit.resolved_pair_probability +
                     x.below.resolved_pair_probability +
                     x.above.resolved_pair_probability,
                 "model resolved probability is segment sum");
    require_near(x.total.quadrature_error_m3_s,
                 x.fit.quadrature_error_m3_s +
                     x.below.quadrature_error_m3_s +
                     x.above.quadrature_error_m3_s,
                 "model quadrature error is segment sum");
}

std::array<double, FUSION_SPECIES_COUNT> canonical_masses() {
    std::array<double, FUSION_SPECIES_COUNT> masses{};
    for (int id = 0; id < FUSION_SPECIES_COUNT; ++id) {
        fusion_nuclear_mass_v1 record{};
        require(fusion_c_nuclear_mass(id, &record) == PB11_STATUS_OK,
                "canonical nuclear mass status");
        require(std::isfinite(record.mass_kg) && record.mass_kg > 0,
                "canonical nuclear mass is positive");
        require(std::isfinite(record.rest_energy_J) && record.rest_energy_J > 0,
                "canonical rest energy is positive");
        masses[id] = record.mass_kg;
    }
    return masses;
}

struct ReactantMasses {
    double a;
    double b;
};

ReactantMasses channel_masses(
    const std::array<double, FUSION_SPECIES_COUNT> &masses, int channel) {
    fusion_nuclear_channel_v1 record{};
    require(fusion_c_nuclear_channel(channel, &record) == PB11_STATUS_OK,
            "canonical channel record status");
    require(record.reactant_ids[0] >= 0 &&
                record.reactant_ids[0] < FUSION_SPECIES_COUNT &&
                record.reactant_ids[1] >= 0 &&
                record.reactant_ids[1] < FUSION_SPECIES_COUNT,
            "canonical channel reactant IDs");
    return {masses[record.reactant_ids[0]], masses[record.reactant_ids[1]]};
}

struct ThermalReference {
    int channel;
    double temperature_keV;
    double k_fit;
    double k_low;
    double k_high;
    double k_high_flat;
    double mean_relative_fit_keV;
    double mean_relative_low_keV;
    double mean_relative_high_keV;
};

// Values transcribed from the versioned
// docs/validation/nuclear-windows/nuclear-window-study-refined.csv.  Keep
// this table in the test so a runtime checkout of the validation artifact is
// not required.
constexpr ThermalReference thermal_references[] = {
    {0, .03, 1.963384623079805e-90, 0, 0, 0,
     1.7444063507149148, 0, 0},
    {0, 3, 3.7163515171599253e-33, 0, 0, 0,
     39.678501672330242, 0, 0},
    {0, 500, 5.6042653004648081e-22, 0, 9.9478443011483924e-31,
     1.0078723976068963e-30, 927.42203783063167, 0, 10277.477368273931},
    {1, .03, 6.6814862890433519e-46, 8.8261366296951495e-47, 0, 0,
     .65328877054852263, .45742305421038232, 0},
    {1, 3, 1.6019568268440368e-26, 4.9597506859418059e-43, 0, 0,
     15.690832890919889, .47997223901065095, 0},
    {1, 500, 8.705215809072417e-23, 2.7024733422681864e-46,
     3.8929005969748553e-26, 4.1988639793522838e-26,
     1076.0520034330116, .48009013398147132, 5508.7596485479517},
    {2, .03, 6.4660003949774518e-46, 8.5384125857996043e-47, 0, 0,
     .65331202952023359, .45742305421038226, 0},
    {2, 3, 1.6009601112408587e-26, 4.798067314842419e-43, 0, 0,
     15.757566236579143, .47997223901065117, 0},
    {2, 500, 1.0373050381587637e-22, 2.6143751639620376e-46,
     4.962825904807019e-26, 5.3598008304076533e-26,
     1051.1155052219879, .48009013398147132, 5408.9928012985165},
    {3, .03, 3.1692853641367208e-45, 2.0118784332273387e-46, 0, 0,
     .68192283053891756, .4630110121432916, 0},
    {3, 3, 1.8742675712238179e-24, 1.2511251194131282e-42, 0, 0,
     17.189379321489419, .4815374705271307, 0},
    {3, 500, 3.694833424110943e-22, 6.8206589582466453e-46,
     5.93060420253042e-26, 6.4092802291801366e-26,
     450.2462958916891, .48163855557542423, 5210.4092091217781},
    {4, .03, 8.4305640877198062e-62, 3.830757867320745e-78, 0, 0,
     1.0458526366945746, .29466501844825049, 0},
    {4, 3, 2.6383113724428645e-28, 6.4875834224708537e-77, 0, 0,
     24.75693394956501, .2954290005825736, 0},
    {4, 500, 2.5202450868684062e-22, 3.3252279990153315e-80,
     4.0390796253262967e-26, 4.2593193278156832e-26,
     612.18238679945102, .2954356110988498, 5320.5833390204507},
};

void thermal_endpoint_reference(
    const std::array<double, FUSION_SPECIES_COUNT> &masses) {
    for (const ThermalReference &reference : thermal_references) {
        const ReactantMasses reactants = channel_masses(masses, reference.channel);
        fusion_rate_model_v1 result{};
        require(fusion_c_thermal_pair_maxwellian_model(
                    reference.channel, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                    reactants.a, reactants.b,
                    reference.temperature_keV * kev,
                    reference.temperature_keV * kev, &result) == PB11_STATUS_OK,
                "endpoint-S thermal model status");
        check_model(result);
        check_segment_sum(result);

        require_near(result.fit.resolved_reactivity_m3_s, reference.k_fit,
                     "endpoint-S fit K reference");
        require_near(result.below.resolved_reactivity_m3_s, reference.k_low,
                     "endpoint-S low K reference");
        require_near(result.above.resolved_reactivity_m3_s, reference.k_high,
                     "endpoint-S high K reference");
        require_near(result.total.resolved_reactivity_m3_s,
                     reference.k_fit + reference.k_low + reference.k_high,
                     "endpoint-S total K reference");

        require_near(result.fit.relative_energy_reactivity_J_m3_s,
                     reference.k_fit * reference.mean_relative_fit_keV * kev,
                     "endpoint-S fit relative moment reference");
        require_near(result.below.relative_energy_reactivity_J_m3_s,
                     reference.k_low * reference.mean_relative_low_keV * kev,
                     "endpoint-S low relative moment reference");
        require_near(result.above.relative_energy_reactivity_J_m3_s,
                     reference.k_high * reference.mean_relative_high_keV * kev,
                     "endpoint-S high relative moment reference");
        require_near(result.total.relative_energy_reactivity_J_m3_s,
                     reference.k_fit * reference.mean_relative_fit_keV * kev +
                         reference.k_low * reference.mean_relative_low_keV * kev +
                         reference.k_high * reference.mean_relative_high_keV * kev,
                     "endpoint-S total relative moment reference");

    }
}

void thermal_high_flat_reference(
    const std::array<double, FUSION_SPECIES_COUNT> &masses) {
    for (const ThermalReference &reference : thermal_references) {
        const ReactantMasses reactants = channel_masses(masses, reference.channel);
        fusion_rate_model_v1 result{};
        require(fusion_c_thermal_pair_maxwellian_model(
                    reference.channel, FUSION_HIGH_FLAT, FUSION_PB_LOW_TB,
                    reactants.a, reactants.b,
                    reference.temperature_keV * kev,
                    reference.temperature_keV * kev, &result) == PB11_STATUS_OK,
                "high-flat thermal model status");
        check_model(result);
        check_segment_sum(result);
        require_near(result.above.resolved_reactivity_m3_s,
                     reference.k_high_flat, "high-flat K reference");
        require_near(result.total.resolved_reactivity_m3_s,
                     reference.k_fit + reference.k_low +
                         reference.k_high_flat,
                     "high-flat total K reference");
    }
}

struct PbLowReference {
    double temperature_keV;
    double k_tb_low400;
    double k_ns_low400;
    double k_c0_plus12_low400;
    double k_c0_minus12_low400;
};

// These are the pB low400 columns transcribed from the same refined CSV.
constexpr PbLowReference pb_low_references[] = {
    {3, 3.7163515171599253e-33, 3.6953466296433803e-33,
     3.9290786625287101e-33, 3.5036243717911282e-33},
    {100, 2.7788496162201461e-23, 2.6824619530693336e-23,
     2.8836143417892723e-23, 2.6740848906510159e-23},
};

fusion_rate_model_v1 pB_low_result(
    const std::array<double, FUSION_SPECIES_COUNT> &masses, int low,
    double temperature_keV) {
    fusion_rate_model_v1 result{};
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_PB11_3ALPHA, FUSION_ENDPOINT_S, low,
                masses[FUSION_PROTON], masses[FUSION_BORON11],
                temperature_keV * kev, temperature_keV * kev, &result) ==
                PB11_STATUS_OK,
            "pB low-piece thermal model status");
    check_model(result);
    check_segment_sum(result);
    return result;
}

void pB_low_piece_references(
    const std::array<double, FUSION_SPECIES_COUNT> &masses) {
    for (const PbLowReference &reference : pb_low_references) {
        const fusion_rate_model_v1 tb =
            pB_low_result(masses, FUSION_PB_LOW_TB, reference.temperature_keV);
        const fusion_rate_model_v1 ns =
            pB_low_result(masses, FUSION_PB_LOW_NS, reference.temperature_keV);
        const fusion_rate_model_v1 plus = pB_low_result(
            masses, FUSION_PB_LOW_C0_PLUS12, reference.temperature_keV);
        const fusion_rate_model_v1 minus = pB_low_result(
            masses, FUSION_PB_LOW_C0_MINUS12, reference.temperature_keV);

        // The API exposes the complete fit segment, while the source CSV
        // isolates the only changed piece as 0..400 keV.  Thus each complete
        // fit-segment difference must equal the corresponding low400 delta.
        require_near(ns.fit.resolved_reactivity_m3_s -
                         tb.fit.resolved_reactivity_m3_s,
                     reference.k_ns_low400 - reference.k_tb_low400,
                     "pB NS low400 rate difference");
        require_near(plus.fit.resolved_reactivity_m3_s -
                         tb.fit.resolved_reactivity_m3_s,
                     reference.k_c0_plus12_low400 - reference.k_tb_low400,
                     "pB C0 plus low400 rate difference");
        require_near(minus.fit.resolved_reactivity_m3_s -
                         tb.fit.resolved_reactivity_m3_s,
                     reference.k_c0_minus12_low400 - reference.k_tb_low400,
                     "pB C0 minus low400 rate difference");
    }
}

void cold_beam_boundaries() {
    constexpr double ma = 1e-27;
    constexpr double mb = 1e-27;
    for (int channel = 0; channel < FUSION_CHANNEL_COUNT; ++channel) {
        double minimum = 0, maximum = 0;
        require(fusion_c_cross_section_domain(channel, &minimum, &maximum) ==
                    PB11_STATUS_OK,
                "cold boundary domain status");
        const double relative_energies[] = {
            .5 * minimum, minimum, maximum, 2 * maximum};
        for (double relative_energy : relative_energies) {
            fusion_rate_model_v1 result{};
            const double projectile_energy = 2 * relative_energy;
            require(fusion_c_beam_maxwellian_model(
                        channel, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB, ma, mb,
                        projectile_energy, 0, &result) == PB11_STATUS_OK,
                    "cold boundary beam model status");
            check_model(result);
            check_segment_sum(result);

            const int selected = relative_energy == 0 ? 0 :
                relative_energy < minimum ? 1 :
                relative_energy > maximum ? 2 : 0;
            const fusion_beam_window_v1 *segments[] = {
                &result.fit, &result.below, &result.above};
            for (int segment = 0; segment < 3; ++segment) {
                require(segments[segment]->resolved_pair_probability ==
                            (segment == selected ? 1. : 0.),
                        "cold boundary has one resolved segment");
                require(segments[segment]->target_energy_reactivity_J_m3_s ==
                            0,
                        "cold beam target moment is zero");
            }
            require(result.total.target_energy_reactivity_J_m3_s == 0,
                    "cold beam total target moment is zero");

            double sigma = 0;
            require(fusion_c_cross_section_model(
                        channel, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                        relative_energy, &sigma) == PB11_STATUS_OK,
                    "cold boundary model cross section status");
            const double speed = std::sqrt(2 * projectile_energy / ma);
            require_near(result.total.resolved_reactivity_m3_s, sigma * speed,
                         "cold boundary sigma-v");
        }
    }
}

void unequal_temperature_swap(
    const std::array<double, FUSION_SPECIES_COUNT> &masses) {
    const double ma = masses[FUSION_DEUTERON];
    const double mb = masses[FUSION_TRITON];
    const double ta = 5 * kev;
    const double tb = 12 * kev;
    fusion_rate_model_v1 original{};
    fusion_rate_model_v1 swapped{};
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                ma, mb, ta, tb, &original) == PB11_STATUS_OK,
            "unequal-temperature model status");
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mb, ma, tb, ta, &swapped) == PB11_STATUS_OK,
            "swapped unequal-temperature model status");
    check_model(original);
    check_model(swapped);
    check_segment_sum(original);
    check_segment_sum(swapped);
    require_near(swapped.total.resolved_reactivity_m3_s,
                 original.total.resolved_reactivity_m3_s,
                 "temperature/mass swap preserves total K");
    require_near(swapped.total.projectile_energy_reactivity_J_m3_s,
                 original.total.target_energy_reactivity_J_m3_s,
                 "temperature/mass swap exchanges projectile moment");
    require_near(swapped.total.target_energy_reactivity_J_m3_s,
                 original.total.projectile_energy_reactivity_J_m3_s,
                 "temperature/mass swap exchanges target moment");
}

void fast_beam_check(
    const std::array<double, FUSION_SPECIES_COUNT> &masses) {
    fusion_rate_model_v1 result{};
    require(fusion_c_beam_maxwellian_model(
                FUSION_PB11_3ALPHA, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                masses[FUSION_PROTON], masses[FUSION_BORON11],
                100 * kev, 3 * kev, &result) == PB11_STATUS_OK,
            "fast beam model status");
    check_model(result);
    check_segment_sum(result);
    require(result.total.resolved_reactivity_m3_s >= 0 &&
                std::isfinite(result.total.resolved_reactivity_m3_s),
            "fast beam K is nonnegative");
    require(result.total.resolved_reactivity_m3_s > 0,
            "fast beam has a resolved reaction rate");
}

void fill_model(fusion_rate_model_v1 &x, double value) {
    x.total.resolved_reactivity_m3_s = value;
    x.fit.resolved_reactivity_m3_s = value;
    x.below.resolved_reactivity_m3_s = value;
    x.above.resolved_reactivity_m3_s = value;
    x.total.domain_incomplete = 1;
    x.fit.domain_incomplete = 1;
    x.below.domain_incomplete = 1;
    x.above.domain_incomplete = 1;
}

void invalid_inputs() {
    constexpr double mass = 1e-27;
    constexpr double energy = 100 * kev;
    constexpr double temperature = 3 * kev;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double sigma = -1;

    require(fusion_c_cross_section_model(
                FUSION_DT_ALPHAN, 99, FUSION_PB_LOW_TB, energy, &sigma) ==
                PB11_STATUS_INVALID_ARGUMENT && sigma == 0,
            "cross-section invalid continuation clears output");
    sigma = -1;
    require(fusion_c_cross_section_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, 99, energy, &sigma) ==
                PB11_STATUS_INVALID_ARGUMENT && sigma == 0,
            "cross-section invalid low policy clears output");
    sigma = -1;
    require(fusion_c_cross_section_model(
                -1, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB, energy, &sigma) ==
                PB11_STATUS_INVALID_ARGUMENT && sigma == 0,
            "cross-section invalid channel clears output");
    sigma = -1;
    require(fusion_c_cross_section_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                -energy, &sigma) == PB11_STATUS_OUT_OF_RANGE && sigma == 0,
            "cross-section negative energy clears output");
    sigma = -1;
    require(fusion_c_cross_section_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                nan, &sigma) == PB11_STATUS_INVALID_ARGUMENT && sigma == 0,
            "cross-section NaN energy clears output");
    require(fusion_c_cross_section_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                energy, nullptr) == PB11_STATUS_NULL_OUTPUT,
            "cross-section null output");

    fusion_rate_model_v1 result{};
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, 99, FUSION_PB_LOW_TB,
                mass, 3 * mass, energy, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "beam invalid continuation clears output");
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, 99,
                mass, 3 * mass, energy, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "beam invalid low policy clears output");
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                99, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, energy, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "beam invalid channel clears output");
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, -energy, temperature, &result) ==
                PB11_STATUS_OUT_OF_RANGE && zero_model(result),
            "beam negative projectile energy clears output");
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, energy, nan, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "beam NaN target temperature clears output");
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                -mass, 3 * mass, energy, temperature, &result) ==
                PB11_STATUS_OUT_OF_RANGE && zero_model(result),
            "beam negative mass clears output");
    fill_model(result, -7);
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                nan, 3 * mass, energy, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "beam NaN mass clears output");
    require(fusion_c_beam_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, energy, temperature, nullptr) ==
                PB11_STATUS_NULL_OUTPUT,
            "beam null output");

    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, 99, FUSION_PB_LOW_TB,
                mass, 3 * mass, temperature, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "thermal invalid continuation clears output");
    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, 99,
                mass, 3 * mass, temperature, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "thermal invalid low policy clears output");
    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                99, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, temperature, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "thermal invalid channel clears output");
    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, -temperature, temperature, &result) ==
                PB11_STATUS_OUT_OF_RANGE && zero_model(result),
            "thermal negative temperature clears output");
    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, temperature, nan, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "thermal NaN temperature clears output");
    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                -mass, 3 * mass, temperature, temperature, &result) ==
                PB11_STATUS_OUT_OF_RANGE && zero_model(result),
            "thermal negative mass clears output");
    fill_model(result, -7);
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                nan, 3 * mass, temperature, temperature, &result) ==
                PB11_STATUS_INVALID_ARGUMENT && zero_model(result),
            "thermal NaN mass clears output");
    require(fusion_c_thermal_pair_maxwellian_model(
                FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB,
                mass, 3 * mass, temperature, temperature, nullptr) ==
                PB11_STATUS_NULL_OUTPUT,
            "thermal null output");
}
}

int main() {
    try {
        const auto masses = canonical_masses();
        thermal_endpoint_reference(masses);
        thermal_high_flat_reference(masses);
        pB_low_piece_references(masses);
        cold_beam_boundaries();
        unequal_temperature_swap(masses);
        fast_beam_check(masses);
        invalid_inputs();
        std::cout << "All fusion rate-model tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
