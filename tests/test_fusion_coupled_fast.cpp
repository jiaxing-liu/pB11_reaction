#include "fusion_coupled_fast.h"
#include "fusion_coupled_sources.h"
#include "fusion_network.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"
#include "fusion_source_state.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kSpecies = FUSION_SPECIES_COUNT;
constexpr int kChannels = FUSION_CHANNEL_COUNT;
constexpr int kBaths = 7;
constexpr double kElectronVoltJ = 1.602176634e-19;
constexpr double kKeVJ = 1.0e3 * kElectronVoltJ;
constexpr double kMeVJ = 1.0e6 * kElectronVoltJ;

using Species = std::array<double, kSpecies>;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool close_scaled(double actual, double expected, double relative,
                  double absolute = 0.0) {
    if (actual == expected) return true;
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return std::abs(actual - expected) <=
           std::max(absolute, relative * scale);
}

struct Grid {
    std::vector<double> edges;

    int cells() const { return static_cast<int>(edges.size()) - 1; }
};

Grid make_grid() {
    // A small but useful grid: a real zero lower edge, an eV-scale first
    // center for the FP operator, and a linear tail through all charged DT
    // products.  No pB channel is enabled in these tests.
    constexpr int log_cells = 48;
    constexpr int linear_cells = 64;
    constexpr int cells = 1 + log_cells + linear_cells;
    const double first_center = 1.0 * kElectronVoltJ;
    const double first_edge = 2.0 * first_center;
    const double split = 0.5 * kMeVJ;
    const double maximum = 25.0 * kMeVJ;

    Grid grid;
    grid.edges.resize(static_cast<std::size_t>(cells) + 1);
    grid.edges[0] = 0.0;
    grid.edges[1] = first_edge;
    const double log_first = std::log(first_edge);
    const double log_split = std::log(split);
    for (int i = 1; i <= log_cells; ++i) {
        const double fraction = static_cast<double>(i) / log_cells;
        grid.edges[static_cast<std::size_t>(1 + i)] =
            std::exp(log_first + fraction * (log_split - log_first));
    }
    for (int i = 1; i <= linear_cells; ++i) {
        const double fraction = static_cast<double>(i) / linear_cells;
        grid.edges[static_cast<std::size_t>(1 + log_cells + i)] =
            split + fraction * (maximum - split);
    }
    require(grid.cells() == cells, "fast test grid cell count");
    require(grid.edges.front() == 0.0 &&
                grid.edges.back() == 25.0 * kMeVJ,
            "fast test grid bounds");
    for (std::size_t i = 1; i < grid.edges.size(); ++i)
        require(grid.edges[i] > grid.edges[i - 1],
                "fast test grid is strictly increasing");
    return grid;
}

Grid make_subnormal_spill_grid() {
    // Keep a real zero lower boundary and use the same 257-cell geometric
    // grid as the physical fast-D spill fixture.  Its upper edge is 25 MeV;
    // only the far-tail alpha coefficient is subnormal.
    constexpr int cells = 257;
    const double first_edge = 2.0e-6 * kElectronVoltJ;
    const double maximum = 25.0 * kMeVJ;
    const double ratio =
        std::pow(maximum / first_edge, 1.0 / static_cast<double>(cells - 1));

    Grid grid;
    grid.edges.resize(static_cast<std::size_t>(cells) + 1);
    grid.edges[0] = 0.0;
    grid.edges[1] = first_edge;
    for (int edge = 2; edge <= cells; ++edge)
        grid.edges[static_cast<std::size_t>(edge)] =
            first_edge * std::pow(ratio, static_cast<double>(edge - 1));
    grid.edges.back() = maximum;
    require(grid.cells() == cells, "subnormal spill grid cell count");
    require(grid.edges.front() == 0.0 && grid.edges.back() == maximum,
            "subnormal spill grid bounds");
    for (std::size_t edge = 1; edge < grid.edges.size(); ++edge)
        require(grid.edges[edge] > grid.edges[edge - 1],
                "subnormal spill grid is strictly increasing");
    return grid;
}

double cell_center(const Grid& grid, int cell) {
    return 0.5 * (grid.edges[static_cast<std::size_t>(cell)] +
                  grid.edges[static_cast<std::size_t>(cell + 1)]);
}

int nearest_cell(const Grid& grid, double energy) {
    int result = 0;
    double best = std::numeric_limits<double>::infinity();
    for (int cell = 0; cell < grid.cells(); ++cell) {
        const double distance = std::abs(cell_center(grid, cell) - energy);
        if (distance < best) {
            best = distance;
            result = cell;
        }
    }
    return result;
}

fusion_coupled_thermal_options_v1 make_options(
    const std::array<int, kChannels>& channels) {
    fusion_coupled_thermal_options_v1 options{};
    options.birth.relative_max_J = 5.0 * kMeVJ;
    options.birth.cm_max_kT = 40.0;
    options.birth.ground_state_q_J = 91.84 * kKeVJ;
    options.birth.cutoff_J = 0.002 * kMeVJ;
    options.birth.l1_fraction = 0.76;
    options.birth.relative_phase = 0.0;
    options.birth.narrow_peak_fraction = 0.051;
    options.birth.continuum_peak_scale = 1.0;
    options.birth.continuation = FUSION_ENDPOINT_S;
    options.birth.pb_low = FUSION_PB_LOW_TB;
    options.birth.remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY;
    options.birth.broad_mode = 13;
    options.birth.fsci_policy = 0;
    options.birth.relative_order = 8;
    options.birth.cm_order = 8;
    options.birth.nq = 8;
    options.birth.ncos = 8;
    options.max_source_rate_error = 2.0e-3;
    options.max_source_debit_error = 2.0e-3;
    options.handoff_max_L1 = 1.0e-3;
    options.handoff_max_mean_error = 1.0e-3;
    options.handoff_enabled = 0;
    for (int channel = 0; channel < kChannels; ++channel)
        options.channels[channel] = channels[static_cast<std::size_t>(channel)];
    return options;
}

fusion_fast_target_options_v1 make_fast_options(
    const std::array<int, kChannels>& channels) {
    fusion_fast_target_options_v1 options{};
    for (int channel = 0; channel < kChannels; ++channel)
        options.channels[channel] = channels[static_cast<std::size_t>(channel)];
    options.angular_order = 8;
    options.angular_max_exponent = 40.0;
    return options;
}

struct Inputs {
    Grid grid = make_grid();
    Species thermal_number{{0.0, 0.0, 5.0e19, 0.0, 0.0, 0.0}};
    double electron_energy_J_m3 = 0.0;
    double ion_energy_J_m3 = 0.0;
    double electron_density_m3 = 1.0e20;
    Species thermal_charge_squared{{1.0, 1.0, 1.0, 4.0, 4.0, 25.0}};
    int inert_count = 0;
    std::vector<double> coulomb_logs;
    std::vector<double> old_s;
    std::vector<double> old_t;
    std::vector<double> external_birth;
    std::vector<double> escape;
    fusion_coupled_thermal_options_v1 options{};
    fusion_fast_target_options_v1 fast_options{};

    explicit Inputs(bool thermal_dt = false) {
        const double ti = 10.0 * kKeVJ;
        electron_energy_J_m3 = 1.5 * electron_density_m3 * ti;
        double thermal_ions = thermal_number[FUSION_TRITON];
        if (thermal_dt) {
            thermal_number[FUSION_DEUTERON] = 5.0e19;
            thermal_ions += thermal_number[FUSION_DEUTERON];
        }
        ion_energy_J_m3 = 1.5 * thermal_ions * ti;

        const std::size_t fast_size =
            static_cast<std::size_t>(kSpecies) * grid.cells();
        coulomb_logs.assign(static_cast<std::size_t>(kSpecies) * kBaths,
                            15.0);
        old_s.assign(fast_size, 0.0);
        old_t.assign(fast_size, 0.0);
        external_birth.assign(fast_size, 0.0);
        escape.assign(fast_size, 0.0);
        options = make_options({{0, 0, 0, thermal_dt ? 1 : 0, 0}});
        fast_options = make_fast_options({{0, 0, 0, 0, 0}});

        // A representative fast deuteron in both components and one energy
        // cell.  Fast-target DT is then a distinct D+thermal-T pool reaction.
        const int cell = nearest_cell(grid, 100.0 * kKeVJ);
        const std::size_t index = static_cast<std::size_t>(FUSION_DEUTERON) *
                                      grid.cells() +
                                  static_cast<std::size_t>(cell);
        old_s[index] = 6.0e19;
        old_t[index] = 4.0e19;
    }
};

struct Trial {
    Species thermal_number{};
    std::vector<double> s;
    std::vector<double> t;
    fusion_coupled_thermal_v1 result{};
};

Trial make_trial(const Inputs& inputs) {
    Trial trial;
    const std::size_t fast_size =
        static_cast<std::size_t>(kSpecies) * inputs.grid.cells();
    trial.thermal_number.fill(-7.0);
    trial.s.assign(fast_size, -7.0);
    trial.t.assign(fast_size, -7.0);
    trial.result = {};
    return trial;
}

int call_fast(const Inputs& inputs, double dt_s, Trial& trial) {
    return fusion_c_coupled_fast_trial(
        dt_s, &inputs.options, &inputs.fast_options, inputs.grid.cells(),
        inputs.grid.edges.data(), inputs.thermal_number.data(),
        inputs.electron_energy_J_m3, inputs.ion_energy_J_m3,
        inputs.electron_density_m3, inputs.thermal_charge_squared.data(),
        inputs.inert_count, nullptr, inputs.coulomb_logs.data(),
        inputs.old_s.data(), inputs.old_t.data(), inputs.external_birth.data(),
        inputs.escape.data(), trial.thermal_number.data(), trial.s.data(),
        trial.t.data(), &trial.result);
}

int call_legacy(const Inputs& inputs, double dt_s, Trial& trial) {
    return fusion_c_coupled_thermal_trial(
        dt_s, &inputs.options, inputs.grid.cells(), inputs.grid.edges.data(),
        inputs.thermal_number.data(), inputs.electron_energy_J_m3,
        inputs.ion_energy_J_m3, inputs.electron_density_m3,
        inputs.thermal_charge_squared.data(), inputs.inert_count, nullptr,
        inputs.coulomb_logs.data(), inputs.old_s.data(), inputs.old_t.data(),
        inputs.external_birth.data(), inputs.escape.data(),
        trial.thermal_number.data(), trial.s.data(), trial.t.data(),
        &trial.result);
}

double fast_number(const Trial& trial, int cells, int species) {
    double number = 0.0;
    for (int cell = 0; cell < cells; ++cell) {
        const std::size_t index = static_cast<std::size_t>(species) * cells +
                                  static_cast<std::size_t>(cell);
        number += trial.s[index] + trial.t[index];
    }
    return number;
}

double fast_energy(const Grid& grid, const std::vector<double>& s,
                   const std::vector<double>& t, int species) {
    const int cells = grid.cells();
    double energy = 0.0;
    for (int cell = 0; cell < cells; ++cell) {
        const std::size_t index = static_cast<std::size_t>(species) * cells +
                                  static_cast<std::size_t>(cell);
        energy += cell_center(grid, cell) * (s[index] + t[index]);
    }
    return energy;
}

double input_number(const Inputs& inputs, int species) {
    double number = 0.0;
    const int cells = inputs.grid.cells();
    for (int cell = 0; cell < cells; ++cell) {
        const std::size_t index = static_cast<std::size_t>(species) * cells +
                                  static_cast<std::size_t>(cell);
        number += inputs.old_s[index] + inputs.old_t[index];
    }
    return number;
}

bool same_array(const double* left, const double* right, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i)
        if (left[i] != right[i]) return false;
    return true;
}

bool same_ledger(const fusion_source_ledger_v1& left,
                 const fusion_source_ledger_v1& right) {
    return same_array(left.events_m3, right.events_m3, kChannels) &&
           same_array(left.nuclear_born_number_m3,
                      right.nuclear_born_number_m3, kSpecies) &&
           same_array(left.nuclear_born_energy_J_m3,
                      right.nuclear_born_energy_J_m3, kSpecies) &&
           same_array(left.external_born_number_m3,
                      right.external_born_number_m3, kSpecies) &&
           same_array(left.external_born_energy_J_m3,
                      right.external_born_energy_J_m3, kSpecies) &&
           same_array(left.thermal_consumed_number_m3,
                      right.thermal_consumed_number_m3, kSpecies) &&
           same_array(left.thermal_consumed_energy_J_m3,
                      right.thermal_consumed_energy_J_m3, kSpecies) &&
           same_array(left.fast_consumed_number_m3,
                      right.fast_consumed_number_m3, kSpecies) &&
           same_array(left.fast_consumed_energy_J_m3,
                      right.fast_consumed_energy_J_m3, kSpecies) &&
           same_array(left.escaped_number_m3, right.escaped_number_m3,
                      kSpecies) &&
           same_array(left.escaped_energy_J_m3, right.escaped_energy_J_m3,
                      kSpecies) &&
           same_array(left.handed_off_number_m3,
                      right.handed_off_number_m3, kSpecies) &&
           same_array(left.handed_off_energy_J_m3,
                      right.handed_off_energy_J_m3, kSpecies) &&
           same_array(left.heat_to_bath_J_m3, right.heat_to_bath_J_m3,
                      kSpecies * kBaths) &&
           left.neutron_number_m3 == right.neutron_number_m3 &&
           left.neutron_energy_J_m3 == right.neutron_energy_J_m3;
}

bool same_result(const fusion_coupled_thermal_v1& left,
                 const fusion_coupled_thermal_v1& right) {
    return same_ledger(left.ledger, right.ledger) &&
           same_array(left.inert_ion_heat_J_m3, right.inert_ion_heat_J_m3,
                      kSpecies) &&
           left.electron_energy_J_m3 == right.electron_energy_J_m3 &&
           left.ion_energy_J_m3 == right.ion_energy_J_m3 &&
           same_array(left.particle_residual_m3, right.particle_residual_m3,
                      kSpecies) &&
           left.energy_residual_J_m3 == right.energy_residual_J_m3 &&
           left.max_source_rate_discrepancy ==
               right.max_source_rate_discrepancy &&
           left.max_source_debit_discrepancy ==
               right.max_source_debit_discrepancy &&
           same_array(left.handoff_L1, right.handoff_L1, kSpecies) &&
           same_array(left.handoff_mean_error, right.handoff_mean_error,
                      kSpecies) &&
           std::equal(std::begin(left.handoff_projected),
                      std::end(left.handoff_projected),
                      std::begin(right.handoff_projected));
}

bool same_diagnostics(const fusion_handoff_diagnostics_v1& left,
                      const fusion_handoff_diagnostics_v1& right) {
    return same_array(left.candidate_number_m3, right.candidate_number_m3,
                      kSpecies) &&
           same_array(left.candidate_energy_J_m3, right.candidate_energy_J_m3,
                      kSpecies) &&
           same_array(left.transferred_number_m3,
                      right.transferred_number_m3, kSpecies) &&
           same_array(left.transferred_energy_J_m3,
                      right.transferred_energy_J_m3, kSpecies) &&
           same_array(left.target_kT_J, right.target_kT_J, kSpecies) &&
           std::equal(std::begin(left.tested), std::end(left.tested),
                      std::begin(right.tested));
}

bool same_usage(const fusion_beam_table_usage_v1& left,
                const fusion_beam_table_usage_v1& right) {
    return left.direct_evaluations == right.direct_evaluations &&
           left.table_evaluations == right.table_evaluations &&
           left.max_validated_rate_error == right.max_validated_rate_error &&
           left.max_validated_debit_error == right.max_validated_debit_error &&
           left.max_validated_number_L1 == right.max_validated_number_L1 &&
           left.max_validated_energy_L1 == right.max_validated_energy_L1;
}

void require_finite_trial(const Trial& trial, const std::string& label) {
    for (double value : trial.thermal_number)
        require(std::isfinite(value) && value >= 0.0,
                label + " thermal output is finite/nonnegative");
    for (double value : trial.s)
        require(std::isfinite(value) && value >= 0.0,
                label + " S output is finite/nonnegative");
    for (double value : trial.t)
        require(std::isfinite(value) && value >= 0.0,
                label + " T output is finite/nonnegative");
    const fusion_source_ledger_v1& ledger = trial.result.ledger;
    for (double value : ledger.events_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " event ledger is finite/nonnegative");
    for (double value : ledger.nuclear_born_number_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " nuclear number ledger is finite/nonnegative");
    for (double value : ledger.nuclear_born_energy_J_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " nuclear energy ledger is finite/nonnegative");
    for (double value : ledger.thermal_consumed_number_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " thermal number ledger is finite/nonnegative");
    for (double value : ledger.thermal_consumed_energy_J_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " thermal energy ledger is finite/nonnegative");
    for (double value : ledger.fast_consumed_number_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " fast number ledger is finite/nonnegative");
    for (double value : ledger.fast_consumed_energy_J_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " fast energy ledger is finite/nonnegative");
    require(std::isfinite(ledger.neutron_number_m3) &&
                ledger.neutron_number_m3 >= 0.0,
            label + " neutron number is finite/nonnegative");
    require(std::isfinite(ledger.neutron_energy_J_m3) &&
                ledger.neutron_energy_J_m3 >= 0.0,
            label + " neutron energy is finite/nonnegative");
    require(std::isfinite(trial.result.electron_energy_J_m3) &&
                std::isfinite(trial.result.ion_energy_J_m3) &&
                std::isfinite(trial.result.energy_residual_J_m3),
            label + " thermal energy result is finite");
}

bool zero_result(const fusion_coupled_thermal_v1& result) {
    fusion_coupled_thermal_v1 zero{};
    return same_result(result, zero);
}

void require_cleared(const Trial& trial, const std::string& label) {
    for (double value : trial.thermal_number)
        require(value == 0.0, label + " clears thermal output");
    for (double value : trial.s)
        require(value == 0.0, label + " clears S output");
    for (double value : trial.t)
        require(value == 0.0, label + " clears T output");
    require(zero_result(trial.result), label + " clears result");
}

void test_fast_dt_accounting_and_shared_components() {
    Inputs inputs;
    inputs.fast_options = make_fast_options({{0, 0, 0, 1, 0}});
    const double dt = 1.0e-3;
    const Species old_thermal = inputs.thermal_number;
    const std::vector<double> old_s = inputs.old_s;
    const std::vector<double> old_t = inputs.old_t;
    const std::vector<double> old_external = inputs.external_birth;
    const std::vector<double> old_escape = inputs.escape;
    Trial trial = make_trial(inputs);
    require(call_fast(inputs, dt, trial) == PB11_STATUS_OK,
            "fast-only DT trial returns OK");
    require_finite_trial(trial, "fast-only DT");

    const fusion_source_ledger_v1& ledger = trial.result.ledger;
    const double events = ledger.events_m3[FUSION_DT_ALPHAN];
    require(events > 0.0, "fast-only DT has nonzero event amount");
    for (int channel = 0; channel < kChannels; ++channel)
        if (channel != FUSION_DT_ALPHAN)
            require(ledger.events_m3[channel] == 0.0,
                    "fast-only DT leaves disabled channels at zero");
    require(ledger.fast_consumed_number_m3[FUSION_DEUTERON] > 0.0 &&
                ledger.fast_consumed_energy_J_m3[FUSION_DEUTERON] > 0.0,
            "fast deuteron consumption is explicit");
    require(ledger.thermal_consumed_number_m3[FUSION_TRITON] > 0.0 &&
                ledger.thermal_consumed_energy_J_m3[FUSION_TRITON] > 0.0,
            "thermal triton consumption is explicit");
    require(ledger.nuclear_born_number_m3[FUSION_HELIUM4] > 0.0 &&
                ledger.nuclear_born_energy_J_m3[FUSION_HELIUM4] > 0.0 &&
                ledger.neutron_number_m3 > 0.0 &&
                ledger.neutron_energy_J_m3 > 0.0,
            "fast-only DT records alpha and neutron products");

    require(close_scaled(ledger.fast_consumed_number_m3[FUSION_DEUTERON],
                         events, 5.0e-10, 1.0e-20),
            "one fast deuteron is consumed per DT event");
    require(close_scaled(ledger.thermal_consumed_number_m3[FUSION_TRITON],
                         events, 5.0e-10, 1.0e-20),
            "one thermal triton is consumed per DT event");
    require(close_scaled(ledger.nuclear_born_number_m3[FUSION_HELIUM4],
                         events, 5.0e-10, 1.0e-20),
            "one alpha is born per DT event");
    require(close_scaled(ledger.neutron_number_m3, events, 5.0e-10,
                         1.0e-20),
            "one neutron is born per DT event");

    // The fast source is a nuclear birth and its external source ledger stays
    // independent.  In this fixture there is no external source at all.
    for (int species = 0; species < kSpecies; ++species)
        require(ledger.external_born_number_m3[species] == 0.0 &&
                    ledger.external_born_energy_J_m3[species] == 0.0,
                "fast nuclear birth is not aliased as external birth");

    fusion_nuclear_channel_v1 dt_channel{};
    require(fusion_c_nuclear_channel(FUSION_DT_ALPHAN, &dt_channel) ==
                PB11_STATUS_OK,
            "DT channel metadata is available");
    const double q_energy = events * dt_channel.q_J;
    const double product_energy =
        ledger.nuclear_born_energy_J_m3[FUSION_HELIUM4] +
        ledger.neutron_energy_J_m3;
    const double fuel_energy =
        ledger.fast_consumed_energy_J_m3[FUSION_DEUTERON] +
        ledger.thermal_consumed_energy_J_m3[FUSION_TRITON];
    require(close_scaled(product_energy, fuel_energy + q_energy, 5.0e-8,
                         1.0e-20),
            "fast DT product/fuel energy closes Q budget");

    const double old_fast_d = input_number(inputs, FUSION_DEUTERON);
    const double new_fast_d = fast_number(trial, inputs.grid.cells(),
                                           FUSION_DEUTERON);
    require(close_scaled(new_fast_d, old_fast_d -
                                      ledger.fast_consumed_number_m3[
                                          FUSION_DEUTERON],
                         5.0e-9, 1.0e-20),
            "fast deuteron number inventory closes");
    require(close_scaled(trial.thermal_number[FUSION_TRITON],
                         old_thermal[FUSION_TRITON] -
                             ledger.thermal_consumed_number_m3[
                                 FUSION_TRITON],
                         5.0e-9, 1.0e-20),
            "thermal triton number inventory closes");
    require(close_scaled(fast_number(trial, inputs.grid.cells(),
                                     FUSION_HELIUM4),
                         ledger.nuclear_born_number_m3[FUSION_HELIUM4],
                         5.0e-8, 1.0e-20),
            "alpha fast number inventory closes");

    double old_fast_energy = 0.0;
    double new_fast_energy = 0.0;
    for (int species = 0; species < kSpecies; ++species) {
        old_fast_energy +=
            fast_energy(inputs.grid, inputs.old_s, inputs.old_t, species);
        new_fast_energy +=
            fast_energy(inputs.grid, trial.s, trial.t, species);
    }
    const double total_energy_change =
        trial.result.electron_energy_J_m3 - inputs.electron_energy_J_m3 +
        trial.result.ion_energy_J_m3 - inputs.ion_energy_J_m3 +
        new_fast_energy - old_fast_energy + ledger.neutron_energy_J_m3;
    require(close_scaled(total_energy_change, q_energy, 5.0e-7, 1.0e-18),
            "fast DT total thermal/kinetic/neutron energy closes Q");
    require(std::abs(trial.result.energy_residual_J_m3) <=
                5.0e-8 * std::max(1.0, std::abs(q_energy)),
            "fast DT coupled energy residual is small");

    // The same hazard is applied to old S and T in a shared cell.  Both
    // components therefore retain positive deuterium after the fast burn;
    // this also guards the implementation against silently dropping T.
    const int d_cell = nearest_cell(inputs.grid, 100.0 * kKeVJ);
    const std::size_t d_index = static_cast<std::size_t>(FUSION_DEUTERON) *
                                    inputs.grid.cells() +
                                static_cast<std::size_t>(d_cell);
    require(trial.s[d_index] > 0.0 && trial.t[d_index] > 0.0,
            "shared fast D S/T components survive in their source cell");

    require(inputs.thermal_number == old_thermal && inputs.old_s == old_s &&
                inputs.old_t == old_t && inputs.external_birth == old_external &&
                inputs.escape == old_escape,
            "fast DT trial leaves all input populations unchanged");
}

void test_all_fast_disabled_exact_legacy_parity() {
    Inputs inputs(/*thermal_dt=*/true);
    const double dt = 1.0e-5;
    Trial legacy = make_trial(inputs);
    Trial additive = make_trial(inputs);
    require(call_legacy(inputs, dt, legacy) == PB11_STATUS_OK,
            "legacy DT fixture returns OK");
    require(call_fast(inputs, dt, additive) == PB11_STATUS_OK,
            "all-fast-disabled additive fixture returns OK");
    require(legacy.thermal_number == additive.thermal_number &&
                legacy.s == additive.s && legacy.t == additive.t &&
                same_result(legacy.result, additive.result),
            "all-fast-disabled API preserves exact legacy full output");
}

void test_invalid_fast_options_clear_outputs() {
    const std::array<int, kChannels> good_channels{{0, 0, 0, 1, 0}};
    for (int kind = 0; kind < 3; ++kind) {
        Inputs inputs;
        inputs.fast_options = make_fast_options(good_channels);
        if (kind == 0) {
            inputs.fast_options.channels[FUSION_DT_ALPHAN] = 2;
        } else if (kind == 1) {
            inputs.fast_options.angular_order = 3;
        } else {
            inputs.fast_options.angular_max_exponent =
                std::numeric_limits<double>::quiet_NaN();
        }
        Trial trial = make_trial(inputs);
        require(call_fast(inputs, 1.0e-4, trial) != PB11_STATUS_OK,
                "invalid fast options are rejected");
        require_cleared(trial, "invalid fast options");
    }
}

void test_deterministic_retry_and_input_immutability() {
    Inputs inputs;
    inputs.fast_options = make_fast_options({{0, 0, 0, 1, 0}});
    const Inputs before = inputs;
    Trial first = make_trial(inputs);
    Trial second = make_trial(inputs);
    require(call_fast(inputs, 1.0e-3, first) == PB11_STATUS_OK,
            "first fast retry trial returns OK");
    require(call_fast(inputs, 1.0e-3, second) == PB11_STATUS_OK,
            "second fast retry trial returns OK");
    require(first.thermal_number == second.thermal_number &&
                first.s == second.s && first.t == second.t &&
                same_result(first.result, second.result),
            "repeated fast trials are bitwise deterministic");
    require(inputs.thermal_number == before.thermal_number &&
                inputs.old_s == before.old_s && inputs.old_t == before.old_t &&
                inputs.external_birth == before.external_birth &&
                inputs.escape == before.escape &&
                inputs.electron_energy_J_m3 == before.electron_energy_J_m3 &&
                inputs.ion_energy_J_m3 == before.ion_energy_J_m3,
            "repeated fast trials do not mutate caller inputs");
}

void test_source_state_atomic_acceptance() {
    Inputs inputs;
    inputs.fast_options = make_fast_options({{0, 0, 0, 1, 0}});
    const double dt = 1.0e-3;
    Trial trial = make_trial(inputs);
    require(call_fast(inputs, dt, trial) == PB11_STATUS_OK,
            "source-state fast trial returns OK");

    fusion_source_state_v1* state = nullptr;
    const std::uint64_t tag = UINT64_C(0x4631544153543031);
    require(fusion_c_source_state_create(
                inputs.grid.cells(), inputs.grid.edges.data(),
                inputs.old_s.data(), inputs.old_t.data(), 2.0, tag, &state) ==
                PB11_STATUS_OK &&
                state != nullptr,
            "source state accepts fast fixture creation");

    std::vector<double> accepted_s(trial.s.size(), 0.0);
    std::vector<double> accepted_t(trial.t.size(), 0.0);
    fusion_source_ledger_v1 cumulative{};
    double accepted_time = 0.0;
    std::uint64_t epoch = 0;
    require(fusion_c_source_state_snapshot(
                state, accepted_s.data(), accepted_t.data(), &cumulative,
                &accepted_time, &epoch) == PB11_STATUS_OK,
            "source state snapshots initial fast state");
    const std::vector<double> initial_s = accepted_s;
    const std::vector<double> initial_t = accepted_t;
    const fusion_source_ledger_v1 initial_ledger = cumulative;
    const double initial_time = accepted_time;
    const std::uint64_t initial_epoch = epoch;

    // Stage is provisional: accepted arrays, cumulative ledger and epoch stay
    // unchanged while a caller decides whether to discard or publish it.
    std::uint64_t discard_ticket = 0;
    require(fusion_c_source_state_begin(state, dt, &discard_ticket) ==
                PB11_STATUS_OK,
            "source state begins provisional fast trial");
    require(fusion_c_source_state_stage(
                state, discard_ticket, trial.s.data(), trial.t.data(),
                &trial.result.ledger) == PB11_STATUS_OK,
            "source state stages provisional fast ledger");
    std::vector<double> pending_s(trial.s.size(), 0.0);
    std::vector<double> pending_t(trial.t.size(), 0.0);
    fusion_source_ledger_v1 pending_ledger{};
    double pending_time = 0.0;
    std::uint64_t pending_epoch = 0;
    require(fusion_c_source_state_snapshot(
                state, pending_s.data(), pending_t.data(), &pending_ledger,
                &pending_time, &pending_epoch) == PB11_STATUS_OK &&
                pending_s == initial_s && pending_t == initial_t &&
                same_ledger(pending_ledger, initial_ledger) &&
                pending_time == initial_time &&
                pending_epoch == initial_epoch,
            "staged fast trial leaves accepted source state unchanged");
    require(fusion_c_source_state_discard(state, discard_ticket) ==
                PB11_STATUS_OK,
            "source state discards provisional fast trial");
    require(fusion_c_source_state_snapshot(
                state, pending_s.data(), pending_t.data(), &pending_ledger,
                &pending_time, &pending_epoch) == PB11_STATUS_OK &&
                pending_s == initial_s && pending_t == initial_t &&
                same_ledger(pending_ledger, initial_ledger) &&
                pending_time == initial_time &&
                pending_epoch == initial_epoch,
            "discarded fast trial leaves accepted source state unchanged");

    std::uint64_t ticket = 0;
    require(fusion_c_source_state_begin(state, dt, &ticket) == PB11_STATUS_OK,
            "source state begins accepted fast trial");
    require(fusion_c_source_state_stage(
                state, ticket, trial.s.data(), trial.t.data(),
                &trial.result.ledger) == PB11_STATUS_OK,
            "source state stages complete fast ledger atomically");
    require(fusion_c_source_state_commit(state, ticket) == PB11_STATUS_OK,
            "source state commits complete fast ledger atomically");
    require(fusion_c_source_state_snapshot(
                state, accepted_s.data(), accepted_t.data(), &cumulative,
                &accepted_time, &epoch) == PB11_STATUS_OK,
            "source state snapshots accepted fast trial");
    require(accepted_s == trial.s && accepted_t == trial.t &&
                same_ledger(cumulative, trial.result.ledger),
            "source state accepts exact fast populations and ledger");
    require(accepted_time == 2.0 + dt && epoch == 1,
            "source state advances fast trial clock once");

    // A repeated commit is rejected and cannot advance the accepted state.
    const std::vector<double> committed_s = accepted_s;
    const std::vector<double> committed_t = accepted_t;
    const fusion_source_ledger_v1 committed_ledger = cumulative;
    const double committed_time = accepted_time;
    const std::uint64_t committed_epoch = epoch;
    require(fusion_c_source_state_commit(state, ticket) != PB11_STATUS_OK,
            "source state rejects repeated fast commit");
    require(fusion_c_source_state_snapshot(
                state, accepted_s.data(), accepted_t.data(), &cumulative,
                &accepted_time, &epoch) == PB11_STATUS_OK &&
                accepted_s == committed_s && accepted_t == committed_t &&
                same_ledger(cumulative, committed_ledger) &&
                accepted_time == committed_time &&
                epoch == committed_epoch,
            "repeated fast commit leaves accepted state unchanged");
    fusion_c_source_state_destroy(state);
}

void test_late_charged_spill_rejects_and_retains_inputs() {
    Inputs inputs;
    inputs.fast_options = make_fast_options({{0, 0, 0, 1, 0}});
    // Keep the lower edge at zero but clip the charged alpha source above
    // 2 MeV.  The DT alpha birth is around 3.5 MeV, so rejection occurs after
    // the fast target/birth path has found a charged spill.
    const int split_edge = 1 + 48;
    for (int edge = split_edge + 1; edge <= inputs.grid.cells(); ++edge) {
        const double fraction = static_cast<double>(edge - split_edge) /
                                (inputs.grid.cells() - split_edge);
        inputs.grid.edges[static_cast<std::size_t>(edge)] =
            0.5 * kMeVJ + fraction * (2.0 * kMeVJ - 0.5 * kMeVJ);
    }
    const Inputs before = inputs;
    Trial trial = make_trial(inputs);
    require(call_fast(inputs, 1.0e-3, trial) != PB11_STATUS_OK,
            "charged fast birth spill rejects complete trial");
    require_cleared(trial, "late charged-spill failure");
    require(inputs.thermal_number == before.thermal_number &&
                inputs.old_s == before.old_s && inputs.old_t == before.old_t &&
                inputs.external_birth == before.external_birth &&
                inputs.escape == before.escape,
            "late charged-spill failure retains caller inputs");
}

Inputs make_subnormal_spill_inputs(double old_fast_deuteron) {
    Inputs inputs;
    inputs.grid = make_subnormal_spill_grid();
    inputs.thermal_number.fill(0.0);
    inputs.thermal_number[FUSION_TRITON] = 1.0e19;
    inputs.electron_density_m3 = 1.0e19;
    const double kT = 10.0 * kKeVJ;
    inputs.electron_energy_J_m3 = 1.5 * inputs.electron_density_m3 * kT;
    inputs.ion_energy_J_m3 = 1.5 * inputs.thermal_number[FUSION_TRITON] * kT;
    inputs.thermal_charge_squared = {{1.0, 1.0, 1.0, 4.0, 4.0, 25.0}};

    const std::size_t fast_size =
        static_cast<std::size_t>(kSpecies) * inputs.grid.cells();
    inputs.coulomb_logs.assign(static_cast<std::size_t>(kSpecies) * kBaths,
                               15.0);
    inputs.old_s.assign(fast_size, 0.0);
    inputs.old_t.assign(fast_size, 0.0);
    inputs.external_birth.assign(fast_size, 0.0);
    inputs.escape.assign(fast_size, 0.0);
    inputs.options = make_options({{0, 0, 0, 0, 0}});
    inputs.fast_options = make_fast_options({{0, 0, 0, 1, 0}});
    inputs.fast_options.angular_order = 16;
    inputs.fast_options.angular_max_exponent = 40.0;

    constexpr double projectile_keV = 4538.1174509712828;
    const int d_cell = nearest_cell(inputs.grid, projectile_keV * kKeVJ);
    require(d_cell == 242, "subnormal spill fixture selects D cell 242");
    require(close_scaled(cell_center(inputs.grid, d_cell),
                         projectile_keV * kKeVJ, 1.0e-12),
            "subnormal spill fixture selects the diagnosed D center");
    const std::size_t d_index = static_cast<std::size_t>(FUSION_DEUTERON) *
                                    inputs.grid.cells() +
                                static_cast<std::size_t>(d_cell);
    inputs.old_s[d_index] = old_fast_deuteron;
    return inputs;
}

void test_subnormal_fast_dt_spill_bound() {
    constexpr double dt = 1.0e-4;
    constexpr double subnormal_old_s = 2.9601203997565862e-291;
    Inputs tiny = make_subnormal_spill_inputs(subnormal_old_s);
    Trial accepted = make_trial(tiny);
    require(call_fast(tiny, dt, accepted) == PB11_STATUS_OK,
            "double-zero fast DT spill is accepted");
    require_finite_trial(accepted, "double-zero fast DT spill");
    require(accepted.result.ledger.events_m3[FUSION_DT_ALPHAN] > 0.0,
            "double-zero fast DT spill still records an event");
    require(accepted.result.ledger.fast_consumed_number_m3[FUSION_DEUTERON] >
                0.0 &&
                accepted.result.ledger.thermal_consumed_number_m3[
                    FUSION_TRITON] > 0.0,
            "double-zero fast DT spill records both reactant debits");

    Inputs representable = make_subnormal_spill_inputs(1.0e15);
    Trial rejected = make_trial(representable);
    require(call_fast(representable, dt, rejected) == PB11_STATUS_OUT_OF_RANGE,
            "representable fast DT spill remains an out-of-range rejection");
    require_cleared(rejected, "representable fast DT spill");
}


void test_diagnosed_fast_interfaces() {
 Inputs inputs;inputs.fast_options=make_fast_options({{0,0,0,1,0}});inputs.fast_options.angular_order=16;
 const double dt_s=1e-3;
 Trial baseline=make_trial(inputs);
 require(call_fast(inputs,dt_s,baseline)==0,"undiagnosed fast reference succeeds");
 const fusion_birth_table_v1* tables[5]{};
 auto evaluate=[&](int mode,Trial&trial,fusion_handoff_diagnostics_v1*diagnostics)->int {
 if(mode==0){
    return fusion_c_coupled_fast_trial_diagnosed(
        dt_s, &inputs.options, &inputs.fast_options, inputs.grid.cells(),
        inputs.grid.edges.data(), inputs.thermal_number.data(),
        inputs.electron_energy_J_m3, inputs.ion_energy_J_m3,
        inputs.electron_density_m3, inputs.thermal_charge_squared.data(),
        inputs.inert_count, nullptr, inputs.coulomb_logs.data(),
        inputs.old_s.data(), inputs.old_t.data(), inputs.external_birth.data(),
        inputs.escape.data(), trial.thermal_number.data(), trial.s.data(),
        trial.t.data(), &trial.result,diagnostics);
 }
 if(mode==1){
    return fusion_c_coupled_fast_table_trial_diagnosed(
        dt_s, &inputs.options, &inputs.fast_options, tables, inputs.grid.cells(),
        inputs.grid.edges.data(), inputs.thermal_number.data(),
        inputs.electron_energy_J_m3, inputs.ion_energy_J_m3,
        inputs.electron_density_m3, inputs.thermal_charge_squared.data(),
        inputs.inert_count, nullptr, inputs.coulomb_logs.data(),
        inputs.old_s.data(), inputs.old_t.data(), inputs.external_birth.data(),
        inputs.escape.data(), trial.thermal_number.data(), trial.s.data(),
        trial.t.data(), &trial.result,diagnostics);
 }
 {
    return fusion_c_coupled_fast_table_trial_effective_charge_diagnosed(
        dt_s, &inputs.options, &inputs.fast_options, tables, inputs.grid.cells(),
        inputs.grid.edges.data(), inputs.thermal_number.data(),
        inputs.electron_energy_J_m3, inputs.ion_energy_J_m3,
        inputs.electron_density_m3, inputs.thermal_charge_squared.data(),
        inputs.inert_count, nullptr, inputs.coulomb_logs.data(),
        inputs.old_s.data(), inputs.old_t.data(), inputs.external_birth.data(),
        inputs.escape.data(), trial.thermal_number.data(), trial.s.data(),
        trial.t.data(), &trial.result,diagnostics);
 }
 };
 for(int mode=0;mode<3;++mode){
  Trial trial=make_trial(inputs);fusion_handoff_diagnostics_v1 d{};
  require(evaluate(mode,trial,&d)==0,"diagnosed fast interface succeeds");
  require(trial.s==baseline.s&&trial.t==baseline.t&&trial.thermal_number==baseline.thermal_number&&same_result(trial.result,baseline.result),"diagnosed fast output parity");
  for(int id=0;id<6;++id){long double sum=0;for(int j=0;j<inputs.grid.cells();++j)sum+=trial.t[id*inputs.grid.cells()+j];
   require(std::abs(sum-d.candidate_number_m3[id])<=1e-14L*std::max(sum,static_cast<long double>(d.candidate_number_m3[id])),"diagnostic candidate matches pre-handoff T");
   require(d.tested[id]==0&&d.target_kT_J[id]==0,"disabled handoff diagnostics remain explicit");
  }
  require(d.transferred_number_m3[FUSION_DEUTERON]>0,"S-to-T diagnostic retained");
  require(evaluate(mode,trial,nullptr)==PB11_STATUS_NULL_OUTPUT&&same_result(trial.result,fusion_coupled_thermal_v1{}),"null diagnostic rejects and clears");
  inputs.fast_options.angular_order=3;
  for(int id=0;id<6;++id){d.candidate_number_m3[id]=1;d.candidate_energy_J_m3[id]=1;d.transferred_number_m3[id]=1;d.transferred_energy_J_m3[id]=1;d.target_kT_J[id]=1;d.tested[id]=1;}
  require(evaluate(mode,trial,&d)!=0,"bad fast options reject diagnostics");
  for(int id=0;id<6;++id)require(d.candidate_number_m3[id]==0&&d.candidate_energy_J_m3[id]==0&&d.transferred_number_m3[id]==0&&d.transferred_energy_J_m3[id]==0&&d.target_kT_J[id]==0&&d.tested[id]==0,"all diagnostics clear on failure");
  inputs.fast_options.angular_order=16;
 }
}

fusion_beam_birth_options_v1 beam_options(const Inputs& in) {
    fusion_beam_birth_options_v1 b{};
    const auto& t=in.options.birth;
    b.relative_max_J=t.relative_max_J;b.angular_max_exponent=in.fast_options.angular_max_exponent;
    b.ground_state_q_J=t.ground_state_q_J;b.cutoff_J=t.cutoff_J;
    b.l1_fraction=t.l1_fraction;b.relative_phase=t.relative_phase;
    b.narrow_peak_fraction=t.narrow_peak_fraction;b.continuum_peak_scale=t.continuum_peak_scale;
    b.continuation=t.continuation;b.pb_low=t.pb_low;b.remainder_policy=t.remainder_policy;
    b.broad_mode=t.broad_mode;b.fsci_policy=t.fsci_policy;b.relative_order=t.relative_order;
    b.angular_order=in.fast_options.angular_order;b.nq=t.nq;b.ncos=t.ncos;
    return b;
}

void test_explicit_full_source_tables() {
    Inputs in;
    in.fast_options.channels[3]=1;
    const int n=in.grid.cells(),cell=nearest_cell(in.grid,100*kKeVJ);
    for(int j=cell-1;j<=cell+1;++j){in.old_s[n+j]=2e18;in.old_t[n+j]=1e18;}
    const Inputs original=in;
    auto options=beam_options(in);
    fusion_birth_table_control_v1 control{};
    control.max_rate_error=control.max_debit_error=1e-5;
    control.max_number_L1=control.max_energy_L1=1e-4;
    control.max_direct_rate_discrepancy=control.max_direct_debit_discrepancy=2e-3;
    control.max_knots=512;control.max_evaluations=4096;control.max_depth=16;
    const double lower=double(static_cast<long double>(in.ion_energy_J_m3)/(1.5L*in.thermal_number[2]));
    std::array<fusion_beam_birth_table_v1*,3> tables{};
    struct Cleanup {std::array<fusion_beam_birth_table_v1*,3>& tables;
        ~Cleanup(){for(auto*t:tables)fusion_c_beam_birth_table_destroy(t);}} cleanup{tables};
    std::vector<fusion_beam_table_entry_v1> entries;
    for(int j=cell-1;j<=cell+1;++j){
        double energy=double((static_cast<long double>(in.grid.edges[j])+in.grid.edges[j+1])/2);
        int status=fusion_c_beam_birth_table_create(3,0,energy,lower,1.05*lower,&options,&control,
            n,in.grid.edges.data(),&tables[j-cell+1]);
        require(status==0,"coupled beam table construction status="+std::to_string(status));
        entries.push_back({3,0,j,tables[j-cell+1]});
    }
    fusion_handoff_diagnostics_v1 diagnostic{};
    fusion_beam_table_usage_v1 usage{};
    auto evaluate=[&](const std::vector<fusion_beam_table_entry_v1>& list,Trial& out,int effective=0){
        return fusion_c_coupled_sources_trial(1e-4,&in.options,&in.fast_options,nullptr,
            int(list.size()),list.empty()?nullptr:list.data(),effective,n,in.grid.edges.data(),
            in.thermal_number.data(),in.electron_energy_J_m3,in.ion_energy_J_m3,in.electron_density_m3,
            in.thermal_charge_squared.data(),0,nullptr,in.coulomb_logs.data(),in.old_s.data(),in.old_t.data(),
            in.external_birth.data(),in.escape.data(),out.thermal_number.data(),out.s.data(),out.t.data(),
            &out.result,&diagnostic,&usage);
    };
    auto evaluate_covered=[&](const std::vector<fusion_beam_table_entry_v1>& list,Trial& out,
                              int policy,uint64_t* outside,
                              const fusion_coupled_floor_limits_v1* limits=nullptr,
                              fusion_birth_floor_ledger_v1* floor=nullptr){
        return fusion_c_coupled_sources_covered_trial(1e-4,&in.options,&in.fast_options,nullptr,
            int(list.size()),list.empty()?nullptr:list.data(),0,n,in.grid.edges.data(),
            in.thermal_number.data(),in.electron_energy_J_m3,in.ion_energy_J_m3,in.electron_density_m3,
            in.thermal_charge_squared.data(),0,nullptr,in.coulomb_logs.data(),in.old_s.data(),in.old_t.data(),
            in.external_birth.data(),in.escape.data(),out.thermal_number.data(),out.s.data(),out.t.data(),
            &out.result,&diagnostic,&usage,limits,floor,policy,outside);
    };
    Trial direct=make_trial(in),empty=make_trial(in),tabulated=make_trial(in);
    require(call_fast(in,1e-4,direct)==0,"source-table direct reference");
    require(evaluate({},empty)==0&&same_result(direct.result,empty.result)&&
        direct.s==empty.s&&direct.t==empty.t&&direct.thermal_number==empty.thermal_number,
        "empty beam set exactly preserves direct arithmetic");
    require(usage.direct_evaluations==3&&usage.table_evaluations==0,"direct calls counted");
    require(evaluate(entries,tabulated)==0,"beam-table endpoint coupled trial");
    require(tabulated.s==direct.s&&tabulated.t==direct.t&&tabulated.thermal_number==direct.thermal_number&&
        same_ledger(tabulated.result.ledger,direct.result.ledger)&&
        tabulated.result.electron_energy_J_m3==direct.result.electron_energy_J_m3&&
        tabulated.result.ion_energy_J_m3==direct.result.ion_energy_J_m3,"endpoint complete physical-state equality");
    require(usage.table_evaluations==3&&usage.direct_evaluations==0,"table calls counted");
    require(usage.max_validated_rate_error<=control.max_rate_error&&
        usage.max_validated_debit_error<=control.max_debit_error&&
        usage.max_validated_number_L1<=control.max_number_L1&&
        usage.max_validated_energy_L1<=control.max_energy_L1,"table envelopes exposed separately");
    auto subset=entries;subset.pop_back();
    require(evaluate(subset,tabulated)==0&&usage.table_evaluations==2&&usage.direct_evaluations==1&&
        tabulated.s==direct.s&&tabulated.t==direct.t,"explicit sparse table/direct mixture");
    for(double fraction:{.013,.029,.043}){
        in.ion_energy_J_m3=original.ion_energy_J_m3*(1+fraction);
        require(call_fast(in,1e-4,direct)==0&&evaluate(entries,tabulated)==0,"independent coupled temperature");
        require_finite_trial(tabulated,"table coupled state");
        auto l1=[](const double*a,const double*b,int count){
            long double difference=0,norm=0;
            for(int j=0;j<count;++j){difference+=std::abs(static_cast<long double>(a[j])-b[j]);norm+=std::abs(static_cast<long double>(a[j]));}
            return norm==0?(difference==0?0.:1.):double(difference/norm);
        };
        for(int id=0;id<6;++id){
            require(l1(direct.s.data()+id*n,tabulated.s.data()+id*n,n)<1e-3,"per-species S full-shape convergence");
            require(l1(direct.t.data()+id*n,tabulated.t.data()+id*n,n)<1e-3,"per-species T full-shape convergence");
            require(l1(direct.result.ledger.heat_to_bath_J_m3+7*id,tabulated.result.ledger.heat_to_bath_J_m3+7*id,7)<1e-3,"per-projectile heat partition convergence");
        }
        require(close_scaled(direct.result.ledger.events_m3[3],tabulated.result.ledger.events_m3[3],1e-4),"coupled event convergence");
        require(close_scaled(direct.result.ledger.neutron_energy_J_m3,tabulated.result.ledger.neutron_energy_J_m3,1e-4),"neutron energy convergence");
        require(close_scaled(direct.result.ledger.thermal_consumed_energy_J_m3[2],tabulated.result.ledger.thermal_consumed_energy_J_m3[2],1e-4),"selected target debit convergence");
        Trial retry=make_trial(in);require(evaluate(entries,retry)==0&&retry.s==tabulated.s&&retry.t==tabulated.t&&same_result(retry.result,tabulated.result),"immutable table retry determinism");
    }
    // Thermal burning changes the target Ti before the old-fast source is sampled.
    in=original;in.options.channels[3]=1;in.thermal_number[1]=5e19;
    in.ion_energy_J_m3=double(1.5L*(static_cast<long double>(in.thermal_number[1])+in.thermal_number[2])*lower*1.02L);
    require(call_fast(in,1e-4,direct)==0&&evaluate(entries,tabulated)==0,"thermal then beam-table split");
    require(close_scaled(direct.result.ledger.events_m3[3],tabulated.result.ledger.events_m3[3],1e-4),"combined thermal/fast events");
    require(usage.table_evaluations==3&&usage.direct_evaluations==0,"post-thermal source uses selected beam tables");
    in=original;

    auto same_floor=[](const fusion_birth_floor_ledger_v1&a,const fusion_birth_floor_ledger_v1&b){
        for(int i=0;i<6;++i){
            if(a.born_number_m3[i]!=b.born_number_m3[i]||a.born_energy_J_m3[i]!=b.born_energy_J_m3[i]||
               a.mapped_number_m3[i]!=b.mapped_number_m3[i]||a.mapped_energy_J_m3[i]!=b.mapped_energy_J_m3[i]||
               a.ion_energy_correction_J_m3[i]!=b.ion_energy_correction_J_m3[i]||
               a.energy_residual_J_m3[i]!=b.energy_residual_J_m3[i])return false;
        }
        return a.remaining_ion_energy_J_m3==b.remaining_ion_energy_J_m3;
    };
    auto covered_rejected=[&](const std::vector<fusion_beam_table_entry_v1>&list,const std::string&label,
                             int policy=FUSION_BEAM_TABLE_STRICT,
                             const fusion_coupled_floor_limits_v1*limits=nullptr,
                             fusion_birth_floor_ledger_v1*floor=nullptr,bool null_counter=false){
        Trial out=make_trial(in);usage={7,7,7,7,7,7};diagnostic={};diagnostic.candidate_number_m3[0]=7;
        uint64_t counter=77;
        if(floor){*floor={};floor->born_number_m3[0]=7;floor->remaining_ion_energy_J_m3=7;}
        int status=evaluate_covered(list,out,policy,null_counter?nullptr:&counter,limits,floor);
        require(status!=PB11_STATUS_OK,label+" rejects");require_cleared(out,label);
        require(usage.direct_evaluations==0&&usage.table_evaluations==0&&usage.max_validated_rate_error==0&&
            usage.max_validated_debit_error==0&&usage.max_validated_number_L1==0&&usage.max_validated_energy_L1==0&&
            same_diagnostics(diagnostic,fusion_handoff_diagnostics_v1{}),label+" clears diagnostics/usage");
        if(floor)require(same_floor(*floor,fusion_birth_floor_ledger_v1{}),label+" clears floor ledger");
        if(!null_counter)require(counter==0,label+" clears outside counter");
    };

    // The old entry point still rejects a used table outside its domain.  The
    // covered entry point preserves that strict behavior and adds an opt-in
    // direct-source policy with complete output parity.
    in=original;
    Trial legacy_inside=make_trial(in),strict_inside=make_trial(in),direct_inside=make_trial(in);
    uint64_t outside_counter=99;
    require(evaluate(entries,legacy_inside)==0,"legacy in-domain covered fixture");
    const auto legacy_diagnostic=diagnostic;const auto legacy_usage=usage;
    require(evaluate_covered(entries,strict_inside,FUSION_BEAM_TABLE_STRICT,&outside_counter)==0,
        "covered strict in-domain fixture");
    const auto strict_diagnostic=diagnostic;const auto strict_usage=usage;
    require(outside_counter==0&&same_result(legacy_inside.result,strict_inside.result)&&
        legacy_inside.thermal_number==strict_inside.thermal_number&&legacy_inside.s==strict_inside.s&&
        legacy_inside.t==strict_inside.t&&same_diagnostics(legacy_diagnostic,strict_diagnostic)&&
        same_usage(legacy_usage,strict_usage),"covered strict is legacy-exact inside the table");
    outside_counter=99;
    require(evaluate_covered(entries,direct_inside,FUSION_BEAM_TABLE_DIRECT_OUTSIDE,&outside_counter)==0,
        "covered direct policy in-domain fixture");
    require(outside_counter==0&&same_result(legacy_inside.result,direct_inside.result)&&
        legacy_inside.thermal_number==direct_inside.thermal_number&&legacy_inside.s==direct_inside.s&&
        legacy_inside.t==direct_inside.t&&same_diagnostics(legacy_diagnostic,diagnostic)&&
        same_usage(legacy_usage,usage),"covered direct policy is legacy-exact inside the table");

    in=original;in.ion_energy_J_m3*=1.1;
    Trial legacy_outside=make_trial(in),strict_outside=make_trial(in);
    require(evaluate(entries,legacy_outside)==PB11_STATUS_OUT_OF_RANGE,
        "legacy used out-of-domain table remains an error");
    require_cleared(legacy_outside,"legacy out-of-domain table");
    outside_counter=77;
    require(evaluate_covered(entries,strict_outside,FUSION_BEAM_TABLE_STRICT,&outside_counter)==PB11_STATUS_OUT_OF_RANGE,
        "covered strict used out-of-domain table remains an error");
    require_cleared(strict_outside,"covered strict out-of-domain table");
    require(outside_counter==0&&usage.direct_evaluations==0&&usage.table_evaluations==0&&
        same_diagnostics(diagnostic,fusion_handoff_diagnostics_v1{}),"strict failure clears coverage diagnostics");


    const auto original_outside_options=in.options;
    in.options.max_source_rate_error=0;in.options.max_source_debit_error=0;
    Trial gate_direct=make_trial(in),gate_covered=make_trial(in);
    const int gate_status=evaluate({},gate_direct);
    require(gate_status==PB11_STATUS_NUMERICAL_FAILURE,"zero source gate direct rejection");
    outside_counter=99;
    require(evaluate_covered(entries,gate_covered,FUSION_BEAM_TABLE_DIRECT_OUTSIDE,&outside_counter)==gate_status,
        "coverage policy must not mask direct numerical gate errors");
    require_cleared(gate_covered,"covered direct gate rejection");
    require(outside_counter==0&&usage.direct_evaluations==0&&usage.table_evaluations==0,
        "direct gate failure clears provisional coverage counts");
    in.options=original_outside_options;
    Trial direct_reference=make_trial(in),direct_fallback=make_trial(in);
    const int direct_reference_status=evaluate({},direct_reference);
    require(direct_reference_status==0,"direct source reference outside table status="+
        std::to_string(direct_reference_status));
    const auto reference_diagnostic=diagnostic;const auto reference_usage=usage;
    outside_counter=0;
    require(evaluate_covered(entries,direct_fallback,FUSION_BEAM_TABLE_DIRECT_OUTSIDE,&outside_counter)==0,
        "covered direct fallback outside table");
    require(outside_counter>0&&outside_counter<=static_cast<uint64_t>(usage.direct_evaluations)&&
        usage.table_evaluations==0&&same_result(direct_reference.result,direct_fallback.result)&&
        direct_reference.thermal_number==direct_fallback.thermal_number&&direct_reference.s==direct_fallback.s&&
        direct_reference.t==direct_fallback.t&&same_diagnostics(reference_diagnostic,diagnostic)&&
        same_usage(reference_usage,usage),"direct fallback has full direct-source parity and bounded counter");

    in=original;
    Trial empty_covered=make_trial(in),empty_legacy=make_trial(in);outside_counter=99;
    require(evaluate_covered({},empty_covered,FUSION_BEAM_TABLE_DIRECT_OUTSIDE,&outside_counter)==0,
        "covered empty entries direct source");
    const auto empty_diagnostic=diagnostic;const auto empty_usage=usage;
    require(evaluate({},empty_legacy)==0,"legacy empty entries direct source");
    require(outside_counter==0&&same_result(empty_covered.result,empty_legacy.result)&&
        empty_covered.thermal_number==empty_legacy.thermal_number&&empty_covered.s==empty_legacy.s&&
        empty_covered.t==empty_legacy.t&&same_diagnostics(empty_diagnostic,diagnostic)&&
        same_usage(empty_usage,usage),"empty covered entries preserve direct parity");

    fusion_coupled_floor_limits_v1 limits{.5,.5};fusion_birth_floor_ledger_v1 floor{};
    Trial floor_reference=make_trial(in),floor_covered=make_trial(in);
    fusion_birth_floor_ledger_v1 reference_floor{};
    require(fusion_c_coupled_sources_floor_trial(1e-4,&in.options,&in.fast_options,nullptr,
        0,nullptr,0,n,in.grid.edges.data(),in.thermal_number.data(),in.electron_energy_J_m3,
        in.ion_energy_J_m3,in.electron_density_m3,in.thermal_charge_squared.data(),0,nullptr,
        in.coulomb_logs.data(),in.old_s.data(),in.old_t.data(),in.external_birth.data(),in.escape.data(),
        floor_reference.thermal_number.data(),floor_reference.s.data(),floor_reference.t.data(),
        &floor_reference.result,&diagnostic,&usage,&limits,&reference_floor)==0,"legacy valid floor reference");
    const auto floor_diagnostic=diagnostic;const auto floor_usage=usage;
    require(evaluate_covered({},floor_covered,FUSION_BEAM_TABLE_DIRECT_OUTSIDE,&outside_counter,&limits,&floor)==0,
        "covered valid floor call");
    require(outside_counter==0&&same_result(floor_reference.result,floor_covered.result)&&
        floor_reference.thermal_number==floor_covered.thermal_number&&floor_reference.s==floor_covered.s&&
        floor_reference.t==floor_covered.t&&same_floor(reference_floor,floor)&&
        same_diagnostics(floor_diagnostic,diagnostic)&&same_usage(floor_usage,usage),"valid floor legacy parity");
    Trial both_null_floor=make_trial(in);outside_counter=99;
    require(evaluate_covered({},both_null_floor,FUSION_BEAM_TABLE_STRICT,&outside_counter,nullptr,nullptr)==0,
        "both null floor pointers preserve direct source");
    covered_rejected(entries,"floor limits without floor ledger",FUSION_BEAM_TABLE_STRICT,&limits,nullptr);
    covered_rejected(entries,"floor ledger without floor limits",FUSION_BEAM_TABLE_STRICT,nullptr,&floor);
    covered_rejected(entries,"invalid coverage policy",77);
    covered_rejected(entries,"missing mandatory outside counter",FUSION_BEAM_TABLE_STRICT,nullptr,nullptr,true);

    auto covered_bad=entries;covered_bad.push_back(entries[0]);covered_rejected(covered_bad,"covered duplicate entry");
    in.grid.edges[n-1]=std::nextafter(in.grid.edges[n-1],in.grid.edges[n]);
    covered_rejected(entries,"covered bad grid");in=original;
    in.options.birth.cutoff_J*=1.01;covered_rejected(entries,"covered bad source options");in=original;

    // A thermal burn can move the actual fast-source temperature outside a
    // narrow table after the table's initial temperature was in-domain.
    in.options.channels[3]=1;in.thermal_number[FUSION_DEUTERON]=5e19;
    in.ion_energy_J_m3=double(1.5L*(static_cast<long double>(in.thermal_number[FUSION_DEUTERON])+in.thermal_number[FUSION_TRITON])*lower*(1.0L+1.e-12L));
    const long double postburn_initial_ti=static_cast<long double>(in.ion_energy_J_m3)/
        (1.5L*(static_cast<long double>(in.thermal_number[FUSION_DEUTERON])+in.thermal_number[FUSION_TRITON]));
    require(postburn_initial_ti>=static_cast<long double>(lower)&&
        postburn_initial_ti<=1.05L*static_cast<long double>(lower),
        "post-burn fixture starts inside the narrow table domain");
    Trial postburn_strict=make_trial(in),postburn_reference=make_trial(in),postburn_fallback=make_trial(in);
    require(call_fast(in,1e-4,postburn_reference)==0,"post-burn direct source reference");
    outside_counter=77;
    int postburn_status=evaluate_covered(entries,postburn_strict,FUSION_BEAM_TABLE_STRICT,&outside_counter);
    require(postburn_status==PB11_STATUS_OUT_OF_RANGE,"post-burn strict table domain check");
    require_cleared(postburn_strict,"post-burn strict table");
    require(evaluate_covered(entries,postburn_fallback,FUSION_BEAM_TABLE_DIRECT_OUTSIDE,&outside_counter)==0,
        "post-burn direct fallback");
    require(outside_counter>0&&same_result(postburn_reference.result,postburn_fallback.result)&&
        postburn_reference.thermal_number==postburn_fallback.thermal_number&&postburn_reference.s==postburn_fallback.s&&
        postburn_reference.t==postburn_fallback.t,"post-burn fallback matches direct source");
    in=original;

    auto rejected=[&](const std::vector<fusion_beam_table_entry_v1>& list,const std::string& label,int effective=0){
        Trial out=make_trial(in);usage={7,7,7,7,7,7};diagnostic.candidate_number_m3[0]=7;
        require(evaluate(list,out,effective)!=0,label+" rejects");require_cleared(out,label);
        require(usage.direct_evaluations==0&&usage.table_evaluations==0&&usage.max_validated_rate_error==0&&
            usage.max_validated_debit_error==0&&usage.max_validated_number_L1==0&&usage.max_validated_energy_L1==0&&
            diagnostic.candidate_number_m3[0]==0,label+" clears diagnostics/usage");
    };
    auto bad=entries;bad.push_back(entries[0]);rejected(bad,"duplicate entry");
    bad=entries;bad[0].table=nullptr;rejected(bad,"null table");
    bad=entries;bad[0].energy_cell++;rejected(bad,"wrong projectile energy");
    bad=entries;bad[0].projectile_slot=1;rejected(bad,"wrong projectile slot");
    bad=entries;bad[0].channel=1;rejected(bad,"disabled channel");
    bad=entries;bad[0].energy_cell=-1;rejected(bad,"negative cell");
    rejected(entries,"invalid effective charge selector",2);
    in.ion_energy_J_m3*=1.1;rejected(entries,"used table temperature outside domain");in=original;
    // Every source option is part of the exact model key, even for inactive populations.
    for(int selector=0;selector<16;++selector){
        in=original;
        switch(selector){
        case 0:in.options.birth.relative_max_J*=1.01;break;
        case 1:in.fast_options.angular_max_exponent+=1;break;
        case 2:in.options.birth.ground_state_q_J*=1.01;break;
        case 3:in.options.birth.cutoff_J*=1.01;break;
        case 4:in.options.birth.l1_fraction*=.99;break;
        case 5:in.options.birth.relative_phase+=.01;break;
        case 6:in.options.birth.narrow_peak_fraction*=.99;break;
        case 7:in.options.birth.continuum_peak_scale*=.99;break;
        case 8:in.options.birth.remainder_policy=1;break;
        case 9:in.options.birth.broad_mode=1;break;
        case 10:in.options.birth.fsci_policy=1;break;
        case 11:in.options.birth.relative_order+=1;break;
        case 12:in.fast_options.angular_order+=1;break;
        case 13:in.options.birth.nq+=1;break;
        case 14:in.options.birth.ncos+=1;break;
        case 15:in.options.birth.continuation=FUSION_HIGH_FLAT;break;
        }
        std::fill(in.old_s.begin(),in.old_s.end(),0);std::fill(in.old_t.begin(),in.old_t.end(),0);
        rejected(entries,"complete model key "+std::to_string(selector));
    }
    in=original;in.grid.edges[n-1]=std::nextafter(in.grid.edges[n-1],in.grid.edges[n]);
    rejected(entries,"full output grid key");in=original;
    require(in.old_s==original.old_s&&in.old_t==original.old_t,"table trials retain caller inventory");
}

}  // namespace

int main() {
    try {
        test_explicit_full_source_tables();
        test_diagnosed_fast_interfaces();
        test_fast_dt_accounting_and_shared_components();
        test_all_fast_disabled_exact_legacy_parity();
        test_invalid_fast_options_clear_outputs();
        test_deterministic_retry_and_input_immutability();
        test_source_state_atomic_acceptance();
        test_late_charged_spill_rejects_and_retains_inputs();
        test_subnormal_fast_dt_spill_bound();
        std::cout << "PASS: coupled fast DT accounting, parity, rollback, "
                     "source-state acceptance and late spill tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
