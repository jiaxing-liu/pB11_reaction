#include "fusion_birth_table.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

/* Keep the canonical decimal spellings used by the direct source tests.  A
 * one-ulp lower product of 1e3*eV can fall just outside a model endpoint. */
constexpr double kKeVJ = 1.602176634e-16;
constexpr double kMeVJ = 1.602176634e-13;
constexpr int kBirthSpecies = 7;
constexpr int kCells = 128;

constexpr std::array<int, 3> kRepresentativeChannels{{
    FUSION_DT_ALPHAN, FUSION_DD_TP, FUSION_DHE3_ALPHAP}};

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void expect_ok(int status, const std::string& message) {
    check(status == PB11_STATUS_OK,
          message + " (status=" + std::to_string(status) + ")");
}

void expect_rejected(int status, const std::string& message) {
    check(status != PB11_STATUS_OK,
          message + " unexpectedly returned OK");
}

fusion_thermal_birth_options_v1 source_options() {
    fusion_thermal_birth_options_v1 options{};
    options.relative_max_J = 5.0 * kMeVJ;
    options.cm_max_kT = 40.0;
    options.ground_state_q_J = 91.84 * kKeVJ;
    options.cutoff_J = 0.001 * kMeVJ;
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
    options.cm_order = 12;
    options.nq = 4;
    options.ncos = 4;
    return options;
}

fusion_birth_table_control_v1 table_control() {
    fusion_birth_table_control_v1 control{};
    control.max_rate_error = 1.0e-2;
    control.max_debit_error = 1.0e-2;
    control.max_number_L1 = 1.0e-2;
    control.max_energy_L1 = 1.0e-2;
    control.max_direct_rate_discrepancy = 1.0e-5;
    control.max_direct_debit_discrepancy = 1.0e-5;
    control.max_knots = 32;
    control.max_evaluations = 256;
    control.max_depth = 10;
    return control;
}

std::vector<double> wide_energy_edges() {
    std::vector<double> edges(static_cast<std::size_t>(kCells) + 1);
    for (int i = 0; i <= kCells; ++i) {
        edges[static_cast<std::size_t>(i)] =
            25.0 * kMeVJ * static_cast<double>(i) /
            static_cast<double>(kCells);
    }
    return edges;
}

struct TableOwner {
    fusion_birth_table_v1* table = nullptr;

    ~TableOwner() { fusion_c_birth_table_destroy(table); }

    TableOwner() = default;
    TableOwner(const TableOwner&) = delete;
    TableOwner& operator=(const TableOwner&) = delete;
};

struct DirectEvaluation {
    std::vector<double> birth;
    fusion_thermal_birth_v1 result{};
};

struct TableEvaluation {
    std::vector<double> birth;
    fusion_birth_coefficients_v1 result{};
};

DirectEvaluation evaluate_direct(int channel, double kT_J,
                                 const fusion_thermal_birth_options_v1& source,
                                 const std::vector<double>& edges) {
    DirectEvaluation evaluation;
    evaluation.birth.assign(kBirthSpecies * static_cast<std::size_t>(kCells),
                            3.0);
    expect_ok(fusion_c_thermal_birth_grid(
                  channel, kT_J, &source, kCells, edges.data(),
                  evaluation.birth.data(), &evaluation.result),
              "evaluate direct thermal-birth source");
    return evaluation;
}

TableEvaluation evaluate_table(const fusion_birth_table_v1* table,
                                double kT_J, int cells,
                                const std::vector<double>& edges) {
    (void)edges;
    TableEvaluation evaluation;
    evaluation.birth.assign(
        kBirthSpecies * static_cast<std::size_t>(cells), 3.0);
    expect_ok(fusion_c_birth_table_evaluate(
                  table, kT_J, cells, evaluation.birth.data(),
                  &evaluation.result),
              "evaluate birth table");
    return evaluation;
}

void create_table(TableOwner& owner, int channel, double lower_kT_J,
                  double upper_kT_J,
                  const fusion_thermal_birth_options_v1& source,
                  const fusion_birth_table_control_v1& control,
                  const std::vector<double>& edges) {
    expect_ok(fusion_c_birth_table_create(
                  channel, lower_kT_J, upper_kT_J, &source, &control, kCells,
                  edges.data(), &owner.table),
              "construct representative birth table");
    check(owner.table != nullptr, "birth table constructor returns a table");
}

bool close_scaled(double actual, double expected, double relative,
                  double absolute = 0.0) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    return std::abs(actual - expected) <=
           absolute + relative * std::max(std::abs(actual), std::abs(expected));
}

void check_source_options(const fusion_thermal_birth_options_v1& actual,
                          const fusion_thermal_birth_options_v1& expected) {
    check(actual.relative_max_J == expected.relative_max_J,
          "table info copies relative_max_J");
    check(actual.cm_max_kT == expected.cm_max_kT,
          "table info copies cm_max_kT");
    check(actual.ground_state_q_J == expected.ground_state_q_J,
          "table info copies ground_state_q_J");
    check(actual.cutoff_J == expected.cutoff_J,
          "table info copies cutoff_J");
    check(actual.l1_fraction == expected.l1_fraction,
          "table info copies l1_fraction");
    check(actual.relative_phase == expected.relative_phase,
          "table info copies relative_phase");
    check(actual.narrow_peak_fraction == expected.narrow_peak_fraction,
          "table info copies narrow_peak_fraction");
    check(actual.continuum_peak_scale == expected.continuum_peak_scale,
          "table info copies continuum_peak_scale");
    check(actual.continuation == expected.continuation,
          "table info copies continuation");
    check(actual.pb_low == expected.pb_low, "table info copies pb_low");
    check(actual.remainder_policy == expected.remainder_policy,
          "table info copies remainder policy");
    check(actual.broad_mode == expected.broad_mode,
          "table info copies broad mode");
    check(actual.fsci_policy == expected.fsci_policy,
          "table info copies fsci policy");
    check(actual.relative_order == expected.relative_order,
          "table info copies relative order");
    check(actual.cm_order == expected.cm_order,
          "table info copies CM order");
    check(actual.nq == expected.nq, "table info copies nq");
    check(actual.ncos == expected.ncos, "table info copies ncos");
}

void check_table_control(const fusion_birth_table_control_v1& actual,
                         const fusion_birth_table_control_v1& expected) {
    check(actual.max_rate_error == expected.max_rate_error,
          "table info copies max_rate_error");
    check(actual.max_debit_error == expected.max_debit_error,
          "table info copies max_debit_error");
    check(actual.max_number_L1 == expected.max_number_L1,
          "table info copies max_number_L1");
    check(actual.max_energy_L1 == expected.max_energy_L1,
          "table info copies max_energy_L1");
    check(actual.max_direct_rate_discrepancy ==
              expected.max_direct_rate_discrepancy,
          "table info copies max_direct_rate_discrepancy");
    check(actual.max_direct_debit_discrepancy ==
              expected.max_direct_debit_discrepancy,
          "table info copies max_direct_debit_discrepancy");
    check(actual.max_knots == expected.max_knots,
          "table info copies max_knots");
    check(actual.max_evaluations == expected.max_evaluations,
          "table info copies max_evaluations");
    check(actual.max_depth == expected.max_depth,
          "table info copies max_depth");
}

void check_info(const fusion_birth_table_v1* table, int channel,
                double lower_kT_J, double upper_kT_J,
                const fusion_thermal_birth_options_v1& source,
                const fusion_birth_table_control_v1& control) {
    fusion_birth_table_info_v1 info{};
    expect_ok(fusion_c_birth_table_info(table, &info),
              "query birth table info");
    check(info.lower_kT_J == lower_kT_J && info.upper_kT_J == upper_kT_J,
          "table info preserves temperature bounds");
    check(info.channel == channel && info.cells == kCells,
          "table info preserves channel and grid dimensions");
    check(info.knots >= 2 && info.knots <= control.max_knots,
          "table uses bounded adaptive knots");
    check(info.direct_evaluations >= 5 &&
              info.direct_evaluations <= control.max_evaluations,
          "table uses bounded direct evaluations");
    check(std::isfinite(info.max_validated_rate_error) &&
              std::isfinite(info.max_validated_debit_error) &&
              std::isfinite(info.max_validated_number_L1) &&
              std::isfinite(info.max_validated_energy_L1),
          "table validation metrics are finite");
    check(info.max_validated_rate_error >= 0.0 &&
              info.max_validated_rate_error <= control.max_rate_error + 1e-10,
          "table rate gate is met");
    check(info.max_validated_debit_error >= 0.0 &&
              info.max_validated_debit_error <= control.max_debit_error + 1e-10,
          "table debit gate is met");
    check(info.max_validated_number_L1 >= 0.0 &&
              info.max_validated_number_L1 <= control.max_number_L1 + 1e-10,
          "table number L1 gate is met");
    check(info.max_validated_energy_L1 >= 0.0 &&
              info.max_validated_energy_L1 <= control.max_energy_L1 + 1e-10,
          "table energy L1 gate is met");
    check(std::isfinite(info.max_sampled_direct_rate_discrepancy) &&
              std::isfinite(info.max_sampled_direct_debit_discrepancy) &&
              info.max_sampled_direct_rate_discrepancy >= 0.0 &&
              info.max_sampled_direct_rate_discrepancy <=
                  control.max_direct_rate_discrepancy + 1e-10 &&
              info.max_sampled_direct_debit_discrepancy >= 0.0 &&
              info.max_sampled_direct_debit_discrepancy <=
                  control.max_direct_debit_discrepancy + 1e-10,
          "sampled direct source gates are recorded and met");
    check_source_options(info.source, source);
    check_table_control(info.control, control);
}

void check_exact_endpoint(const TableEvaluation& table,
                          const DirectEvaluation& direct) {
    check(table.result.reactivity_m3_s == direct.result.reactivity_m3_s,
          "table endpoint rate is exactly direct");
    for (int i = 0; i < 2; ++i) {
        check(table.result.reactant_energy_moment_J_m3_s[i] ==
                  direct.result.reactant_energy_moment_J_m3_s[i],
              "table endpoint reactant debit is exactly direct");
    }
    for (int id = 0; id < kBirthSpecies; ++id) {
        check(table.result.below_number_m3_s[id] ==
                  direct.result.below_number_m3_s[id],
              "table endpoint below number is exactly direct");
        check(table.result.below_energy_J_m3_s[id] ==
                  direct.result.below_energy_J_m3_s[id],
              "table endpoint below energy is exactly direct");
        check(table.result.above_number_m3_s[id] ==
                  direct.result.above_number_m3_s[id],
              "table endpoint above number is exactly direct");
        check(table.result.above_energy_J_m3_s[id] ==
                  direct.result.above_energy_J_m3_s[id],
              "table endpoint above energy is exactly direct");
    }
    check(table.birth == direct.birth,
          "table endpoint product grid is exactly direct");
}

struct SourceView {
    double reactivity = 0.0;
    std::array<double, 2> reactant_energy{};
    std::array<double, kBirthSpecies> below_number{};
    std::array<double, kBirthSpecies> below_energy{};
    std::array<double, kBirthSpecies> above_number{};
    std::array<double, kBirthSpecies> above_energy{};
};

SourceView view(const fusion_birth_coefficients_v1& result) {
    SourceView output;
    output.reactivity = result.reactivity_m3_s;
    for (int i = 0; i < 2; ++i)
        output.reactant_energy[static_cast<std::size_t>(i)] =
            result.reactant_energy_moment_J_m3_s[i];
    for (int id = 0; id < kBirthSpecies; ++id) {
        output.below_number[static_cast<std::size_t>(id)] =
            result.below_number_m3_s[id];
        output.below_energy[static_cast<std::size_t>(id)] =
            result.below_energy_J_m3_s[id];
        output.above_number[static_cast<std::size_t>(id)] =
            result.above_number_m3_s[id];
        output.above_energy[static_cast<std::size_t>(id)] =
            result.above_energy_J_m3_s[id];
    }
    return output;
}

SourceView view(const fusion_thermal_birth_v1& result) {
    SourceView output;
    output.reactivity = result.reactivity_m3_s;
    for (int i = 0; i < 2; ++i)
        output.reactant_energy[static_cast<std::size_t>(i)] =
            result.reactant_energy_moment_J_m3_s[i];
    for (int id = 0; id < kBirthSpecies; ++id) {
        output.below_number[static_cast<std::size_t>(id)] =
            result.below_number_m3_s[id];
        output.below_energy[static_cast<std::size_t>(id)] =
            result.below_energy_J_m3_s[id];
        output.above_number[static_cast<std::size_t>(id)] =
            result.above_number_m3_s[id];
        output.above_energy[static_cast<std::size_t>(id)] =
            result.above_energy_J_m3_s[id];
    }
    return output;
}

std::array<double, kBirthSpecies> mapped_number(
    const std::vector<double>& birth, const SourceView& source) {
    std::array<double, kBirthSpecies> number{};
    for (int id = 0; id < kBirthSpecies; ++id) {
        for (int cell = 0; cell < kCells; ++cell) {
            number[static_cast<std::size_t>(id)] +=
                birth[static_cast<std::size_t>(id) * kCells + cell];
        }
        number[static_cast<std::size_t>(id)] +=
            source.below_number[static_cast<std::size_t>(id)] +
            source.above_number[static_cast<std::size_t>(id)];
    }
    return number;
}

std::array<double, kBirthSpecies> mapped_energy(
    const std::vector<double>& birth, const SourceView& source,
    const std::vector<double>& edges) {
    std::array<double, kBirthSpecies> energy{};
    for (int id = 0; id < kBirthSpecies; ++id) {
        for (int cell = 0; cell < kCells; ++cell) {
            const double center =
                0.5 * (edges[static_cast<std::size_t>(cell)] +
                       edges[static_cast<std::size_t>(cell + 1)]);
            energy[static_cast<std::size_t>(id)] +=
                center * birth[static_cast<std::size_t>(id) * kCells + cell];
        }
        energy[static_cast<std::size_t>(id)] +=
            source.below_energy[static_cast<std::size_t>(id)] +
            source.above_energy[static_cast<std::size_t>(id)];
    }
    return energy;
}

void check_nonnegative_source(const std::vector<double>& birth,
                              const SourceView& source) {
    for (double value : birth)
        check(std::isfinite(value) && value >= 0.0,
              "birth table grid coefficients are finite and nonnegative");
    for (int id = 0; id < kBirthSpecies; ++id) {
        check(std::isfinite(source.below_number[static_cast<std::size_t>(id)]) &&
                  source.below_number[static_cast<std::size_t>(id)] >= 0.0 &&
                  std::isfinite(source.below_energy[static_cast<std::size_t>(id)]) &&
                  source.below_energy[static_cast<std::size_t>(id)] >= 0.0 &&
                  std::isfinite(source.above_number[static_cast<std::size_t>(id)]) &&
                  source.above_number[static_cast<std::size_t>(id)] >= 0.0 &&
                  std::isfinite(source.above_energy[static_cast<std::size_t>(id)]) &&
                  source.above_energy[static_cast<std::size_t>(id)] >= 0.0,
              "birth table spill coefficients are finite and nonnegative");
    }
}

void check_product_closure(int channel, const std::vector<double>& birth,
                           const SourceView& source,
                           const std::vector<double>& edges) {
    check(std::isfinite(source.reactivity) && source.reactivity > 0.0,
          "representative source reactivity is positive");
    check_nonnegative_source(birth, source);
    fusion_nuclear_channel_v1 reaction{};
    expect_ok(fusion_c_nuclear_channel(channel, &reaction),
              "query representative channel metadata");

    std::array<int, kBirthSpecies> multiplicity{};
    for (int i = 0; i < reaction.product_count; ++i) {
        const int product = reaction.product_ids[i];
        check(product >= 0 && product < kBirthSpecies,
              "channel product ID fits birth output");
        ++multiplicity[static_cast<std::size_t>(product)];
    }
    const auto number = mapped_number(birth, source);
    const auto energy = mapped_energy(birth, source, edges);
    for (int id = 0; id < kBirthSpecies; ++id) {
        const double expected =
            static_cast<double>(multiplicity[static_cast<std::size_t>(id)]) *
            source.reactivity;
        check(close_scaled(number[static_cast<std::size_t>(id)], expected,
                           3.0e-8, 1.0e-300),
              "independent product stoichiometry closes");
    }

    double total_energy = 0.0;
    for (double value : energy) total_energy += value;
    const double expected_energy =
        reaction.q_J * source.reactivity + source.reactant_energy[0] +
        source.reactant_energy[1];
    check(close_scaled(total_energy, expected_energy, 3.0e-8, 1.0e-300),
          "independent Q plus reactant kinetic energy closes products");
}

std::pair<double, double> max_species_L1(const TableEvaluation& table,
                                         const DirectEvaluation& direct,
                                         const std::vector<double>& edges) {
    const SourceView table_source = view(table.result);
    const SourceView direct_source = view(direct.result);
    double maximum_number_L1 = 0.0;
    double maximum_energy_L1 = 0.0;
    for (int id = 0; id < kBirthSpecies; ++id) {
        double number_reference = 0.0;
        double energy_reference = 0.0;
        double number_error = 0.0;
        double energy_error = 0.0;
        for (int cell = 0; cell < kCells; ++cell) {
            const std::size_t index = static_cast<std::size_t>(id) * kCells +
                                      static_cast<std::size_t>(cell);
            const double center =
                0.5 * (edges[static_cast<std::size_t>(cell)] +
                       edges[static_cast<std::size_t>(cell + 1)]);
            number_reference += std::abs(direct.birth[index]);
            number_error += std::abs(table.birth[index] - direct.birth[index]);
            energy_reference += center * std::abs(direct.birth[index]);
            energy_error +=
                center * std::abs(table.birth[index] - direct.birth[index]);
        }
        number_reference += direct_source.below_number[static_cast<std::size_t>(id)] +
                            direct_source.above_number[static_cast<std::size_t>(id)];
        number_error += std::abs(table_source.below_number[static_cast<std::size_t>(id)] -
                                 direct_source.below_number[static_cast<std::size_t>(id)]);
        number_error += std::abs(table_source.above_number[static_cast<std::size_t>(id)] -
                                 direct_source.above_number[static_cast<std::size_t>(id)]);
        energy_reference += direct_source.below_energy[static_cast<std::size_t>(id)] +
                            direct_source.above_energy[static_cast<std::size_t>(id)];
        energy_error += std::abs(table_source.below_energy[static_cast<std::size_t>(id)] -
                                 direct_source.below_energy[static_cast<std::size_t>(id)]);
        energy_error += std::abs(table_source.above_energy[static_cast<std::size_t>(id)] -
                                 direct_source.above_energy[static_cast<std::size_t>(id)]);
        if (number_reference == 0.0) {
            check(number_error == 0.0,
                  "zero-reference product number has exact table agreement");
        } else {
            maximum_number_L1 =
                std::max(maximum_number_L1, number_error / number_reference);
        }
        if (energy_reference == 0.0) {
            check(energy_error == 0.0,
                  "zero-reference product energy has exact table agreement");
        } else {
            maximum_energy_L1 =
                std::max(maximum_energy_L1, energy_error / energy_reference);
        }
    }
    return {maximum_number_L1, maximum_energy_L1};
}

void test_representative_tables() {
    const double lower_kT_J = 20.0 * kKeVJ;
    const double upper_kT_J = 20.5 * kKeVJ;
    const double sample_kT_J =
        std::exp(std::log(lower_kT_J) + 0.37 *
                                      (std::log(upper_kT_J) -
                                       std::log(lower_kT_J)));
    const auto source = source_options();
    const auto control = table_control();
    const auto edges = wide_energy_edges();
    std::array<TableOwner, kRepresentativeChannels.size()> tables;

    for (std::size_t i = 0; i < kRepresentativeChannels.size(); ++i) {
        const int channel = kRepresentativeChannels[i];
        create_table(tables[i], channel, lower_kT_J, upper_kT_J, source,
                     control, edges);
        check_info(tables[i].table, channel, lower_kT_J, upper_kT_J, source,
                   control);

        for (double temperature : {lower_kT_J, upper_kT_J}) {
            const DirectEvaluation direct =
                evaluate_direct(channel, temperature, source, edges);
            check(std::abs(direct.result.relative_rate_discrepancy) <=
                      control.max_direct_rate_discrepancy + 1.0e-10,
                  "direct source rate gate is met");
            check(direct.result.relative_reactant_energy_discrepancy >= 0.0 &&
                      direct.result.relative_reactant_energy_discrepancy <=
                          control.max_direct_debit_discrepancy + 1.0e-10,
                  "direct source debit gate is met");
            const TableEvaluation table =
                evaluate_table(tables[i].table, temperature, kCells, edges);
            check_exact_endpoint(table, direct);
            check_product_closure(channel, table.birth, view(table.result),
                                  edges);
            check_product_closure(channel, direct.birth, view(direct.result),
                                  edges);
        }

        const DirectEvaluation direct =
            evaluate_direct(channel, sample_kT_J, source, edges);
        const TableEvaluation table =
            evaluate_table(tables[i].table, sample_kT_J, kCells, edges);
        check(std::abs(table.result.reactivity_m3_s -
                       direct.result.reactivity_m3_s) /
                  std::max(std::abs(direct.result.reactivity_m3_s), 1.0e-300) <=
                  control.max_rate_error + 1.0e-10,
              "non-knot table rate follows direct source gate");
        for (int reactant = 0; reactant < 2; ++reactant) {
            const double table_debit =
                table.result.reactant_energy_moment_J_m3_s[reactant];
            const double direct_debit =
                direct.result.reactant_energy_moment_J_m3_s[reactant];
            check(std::abs(table_debit - direct_debit) /
                          std::max(std::abs(direct_debit), 1.0e-300) <=
                      control.max_debit_error + 1.0e-10,
                  "non-knot table debit follows direct source gate");
        }
        const auto l1 = max_species_L1(table, direct, edges);
        check(l1.first <= control.max_number_L1 + 1.0e-10,
              "non-knot table maximum per-species number L1 follows gate");
        check(l1.second <= control.max_energy_L1 + 1.0e-10,
              "non-knot table maximum per-species energy L1 follows gate");
        check_product_closure(channel, table.birth, view(table.result), edges);
        check_product_closure(channel, direct.birth, view(direct.result), edges);
    }
}

void poison_coefficients(fusion_birth_coefficients_v1& result) {
    result.reactivity_m3_s = 7.0;
    for (double& value : result.reactant_energy_moment_J_m3_s) value = 7.0;
    for (double& value : result.below_number_m3_s) value = 7.0;
    for (double& value : result.below_energy_J_m3_s) value = 7.0;
    for (double& value : result.above_number_m3_s) value = 7.0;
    for (double& value : result.above_energy_J_m3_s) value = 7.0;
}

void check_zero_coefficients(const fusion_birth_coefficients_v1& result) {
    fusion_birth_coefficients_v1 zero{};
    check(result.reactivity_m3_s == zero.reactivity_m3_s,
          "failed table evaluation clears rate");
    for (int i = 0; i < 2; ++i)
        check(result.reactant_energy_moment_J_m3_s[i] == 0.0,
              "failed table evaluation clears reactant moments");
    for (int id = 0; id < kBirthSpecies; ++id) {
        check(result.below_number_m3_s[id] == 0.0 &&
                  result.below_energy_J_m3_s[id] == 0.0 &&
                  result.above_number_m3_s[id] == 0.0 &&
                  result.above_energy_J_m3_s[id] == 0.0,
              "failed table evaluation clears spill coefficients");
    }
}

void check_zero_birth(const std::vector<double>& birth) {
    for (double value : birth)
        check(value == 0.0, "failed table evaluation clears birth grid");
}

void test_bad_evaluation_inputs(const fusion_birth_table_v1* table,
                                const std::vector<double>& edges) {
    for (double temperature : {19.99 * kKeVJ, 20.51 * kKeVJ}) {
        std::vector<double> birth(kBirthSpecies * static_cast<std::size_t>(kCells),
                                  9.0);
        fusion_birth_coefficients_v1 result{};
        poison_coefficients(result);
        expect_rejected(fusion_c_birth_table_evaluate(
                            table, temperature, kCells, birth.data(), &result),
                        "out-of-range table temperature is rejected");
        check_zero_birth(birth);
        check_zero_coefficients(result);
    }

    const int wrong_cells = kCells - 1;
    std::vector<double> wrong_birth(
        kBirthSpecies * static_cast<std::size_t>(wrong_cells), 9.0);
    fusion_birth_coefficients_v1 wrong_result{};
    poison_coefficients(wrong_result);
    expect_rejected(fusion_c_birth_table_evaluate(
                        table, 20.25 * kKeVJ, wrong_cells, wrong_birth.data(),
                        &wrong_result),
                    "mismatched table grid size is rejected");
    check_zero_birth(wrong_birth);
    check_zero_coefficients(wrong_result);

    std::vector<double> null_table_birth(
        kBirthSpecies * static_cast<std::size_t>(kCells), 9.0);
    fusion_birth_coefficients_v1 null_table_result{};
    poison_coefficients(null_table_result);
    expect_rejected(fusion_c_birth_table_evaluate(
                        nullptr, 20.25 * kKeVJ, kCells,
                        null_table_birth.data(), &null_table_result),
                    "null table is rejected");
    check_zero_birth(null_table_birth);
    check_zero_coefficients(null_table_result);

    fusion_birth_table_info_v1 info{};
    info.channel = 99;
    info.cells = 99;
    expect_rejected(fusion_c_birth_table_info(nullptr, &info),
                    "null table info query is rejected");
    check(info.channel == 0 && info.cells == 0,
          "failed table info query clears output");

    (void)edges;
}

void test_bad_construction_inputs(const std::vector<double>& edges) {
    const auto source = source_options();
    const auto control = table_control();
    const double lower_kT_J = 20.0 * kKeVJ;
    const double upper_kT_J = 20.5 * kKeVJ;
    fusion_birth_table_v1* table = nullptr;

    auto too_few_knots = control;
    too_few_knots.max_knots = 2;
    too_few_knots.max_rate_error = 1.0e-12;
    too_few_knots.max_debit_error = 1.0e-12;
    too_few_knots.max_number_L1 = 1.0e-12;
    too_few_knots.max_energy_L1 = 1.0e-12;
    expect_rejected(fusion_c_birth_table_create(
                        FUSION_DT_ALPHAN, 19.5 * kKeVJ, 23.0 * kKeVJ, &source,
                        &too_few_knots, kCells, edges.data(), &table),
                    "too-small knot budget is rejected");
    check(table == nullptr, "failed small-knot construction returns no table");

    auto too_few_evaluations = control;
    too_few_evaluations.max_evaluations = 4;
    expect_rejected(fusion_c_birth_table_create(
                        FUSION_DT_ALPHAN, lower_kT_J, upper_kT_J, &source,
                        &too_few_evaluations, kCells, edges.data(), &table),
                    "too-small evaluation budget is rejected");
    check(table == nullptr,
          "failed small-evaluation construction returns no table");

    auto strict_zero = control;
    strict_zero.max_rate_error = 0.0;
    strict_zero.max_debit_error = 0.0;
    strict_zero.max_number_L1 = 0.0;
    strict_zero.max_energy_L1 = 0.0;
    strict_zero.max_direct_rate_discrepancy = 0.0;
    strict_zero.max_direct_debit_discrepancy = 0.0;
    strict_zero.max_depth = 0;
    expect_rejected(fusion_c_birth_table_create(
                        FUSION_DT_ALPHAN, lower_kT_J, upper_kT_J, &source,
                        &strict_zero, kCells, edges.data(), &table),
                    "strict zero gates with no refinement are rejected");
    check(table == nullptr,
          "strict zero-gate construction does not leak a partial table");

    auto bad_source = source;
    bad_source.relative_order = 3;
    expect_rejected(fusion_c_birth_table_create(
                        FUSION_DT_ALPHAN, lower_kT_J, upper_kT_J, &bad_source,
                        &control, kCells, edges.data(), &table),
                    "invalid source options are rejected");
    check(table == nullptr, "bad source construction returns no table");

    std::vector<double> bad_edges = edges;
    bad_edges[10] = bad_edges[9];
    expect_rejected(fusion_c_birth_table_create(
                        FUSION_DT_ALPHAN, lower_kT_J, upper_kT_J, &source,
                        &control, kCells, bad_edges.data(), &table),
                    "non-increasing construction grid is rejected");
    check(table == nullptr, "bad grid construction returns no table");
}

void test_independent_tables(const std::vector<double>& edges) {
    const auto source = source_options();
    const auto control = table_control();
    const double lower_kT_J = 20.0 * kKeVJ;
    const double upper_kT_J = 20.5 * kKeVJ;
    const double sample_kT_J =
        std::exp(std::log(lower_kT_J) + 0.37 *
                                      (std::log(upper_kT_J) -
                                       std::log(lower_kT_J)));
    std::array<TableOwner, kRepresentativeChannels.size()> tables;
    for (std::size_t i = 0; i < kRepresentativeChannels.size(); ++i) {
        create_table(tables[i], kRepresentativeChannels[i], lower_kT_J,
                     upper_kT_J, source, control, edges);
    }

    const TableEvaluation before =
        evaluate_table(tables[0].table, sample_kT_J, kCells, edges);
    const TableEvaluation dd =
        evaluate_table(tables[1].table, sample_kT_J, kCells, edges);
    const TableEvaluation dhe3 =
        evaluate_table(tables[2].table, sample_kT_J, kCells, edges);
    check(dd.result.reactivity_m3_s > 0.0 &&
              dhe3.result.reactivity_m3_s > 0.0,
          "independent DD and DHe3 tables remain usable");
    const TableEvaluation after =
        evaluate_table(tables[0].table, sample_kT_J, kCells, edges);
    check(after.birth == before.birth &&
              after.result.reactivity_m3_s == before.result.reactivity_m3_s &&
              after.result.reactant_energy_moment_J_m3_s[0] ==
                  before.result.reactant_energy_moment_J_m3_s[0] &&
              after.result.reactant_energy_moment_J_m3_s[1] ==
                  before.result.reactant_energy_moment_J_m3_s[1],
          "independent table construction does not perturb DT evaluation");

    fusion_c_birth_table_destroy(tables[1].table);
    tables[1].table = nullptr;
    const TableEvaluation after_destroy =
        evaluate_table(tables[0].table, sample_kT_J, kCells, edges);
    check(after_destroy.birth == before.birth &&
              after_destroy.result.reactivity_m3_s ==
                  before.result.reactivity_m3_s,
          "destroying one table does not perturb another table");
}

}  // namespace

int main() {
    try {
        test_representative_tables();
        const auto edges = wide_energy_edges();

        /* Build one valid representative table for the failure-path checks. */
        const auto source = source_options();
        const auto control = table_control();
        TableOwner valid;
        create_table(valid, FUSION_DT_ALPHAN, 20.0 * kKeVJ, 20.5 * kKeVJ,
                     source, control, edges);
        test_bad_evaluation_inputs(valid.table, edges);
        test_bad_construction_inputs(edges);
        test_independent_tables(edges);
        std::cout << "All fusion birth-table tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
