#include "fusion_coupled_thermal.h"
#include "fusion_handoff.h"
#include "fusion_network.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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
constexpr double kAtomicMassUnitKg = 1.66053906660e-27;

using Species = std::array<double, kSpecies>;
using Channels = std::array<int, kChannels>;

void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

bool close_scaled(double actual, double expected, double relative,
                  double absolute = 0.0) {
    if (actual == expected)
        return true;
    if (!std::isfinite(actual) || !std::isfinite(expected))
        return false;
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return std::abs(actual - expected) <=
           std::max(absolute, relative * scale);
}

struct Grid {
    std::vector<double> edges;

    int cells() const { return static_cast<int>(edges.size()) - 1; }
};

Grid make_grid() {
    // The first center is positive but tiny.  This leaves a real zero lower
    // boundary for the FP operator while avoiding an E=0 cell center.  The
    // linear tail reaches 25 MeV, covering charged DT births throughout the
    // 1--20 keV bath range.  The grid is deliberately bounded at 257 cells so
    // the coupled test remains a small integration test.
    constexpr int log_cells = 96;
    constexpr int linear_cells = 160;
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
        grid.edges[1 + i] =
            std::exp(log_first + fraction * (log_split - log_first));
    }
    for (int i = 1; i <= linear_cells; ++i) {
        const double fraction = static_cast<double>(i) / linear_cells;
        grid.edges[1 + log_cells + i] =
            split + fraction * (maximum - split);
    }
    require(grid.cells() >= 200, "coupled grid has at least 200 cells");
    require(grid.edges.front() == 0.0 &&
                grid.edges.back() == 25.0 * kMeVJ,
            "coupled grid has requested zero and 25 MeV bounds");
    for (std::size_t i = 1; i < grid.edges.size(); ++i)
        require(grid.edges[i] > grid.edges[i - 1],
                "coupled grid edges are strictly increasing");
    return grid;
}

int nearest_cell(const Grid& grid, double energy) {
    int result = 0;
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < grid.cells(); ++i) {
        const double center =
            0.5 * (grid.edges[static_cast<std::size_t>(i)] +
                   grid.edges[static_cast<std::size_t>(i + 1)]);
        if (std::abs(center - energy) < best) {
            best = std::abs(center - energy);
            result = i;
        }
    }
    return result;
}

double cell_center(const Grid& grid, int cell) {
    return 0.5 * (grid.edges[static_cast<std::size_t>(cell)] +
                  grid.edges[static_cast<std::size_t>(cell + 1)]);
}

fusion_coupled_thermal_options_v1 make_options(const Channels& channels) {
    fusion_coupled_thermal_options_v1 options{};
    options.birth.relative_max_J = 5.0 * kMeVJ;
    options.birth.cm_max_kT = 40.0;
    options.birth.ground_state_q_J = 91.84 * kKeVJ;
    // Spell the lower bound exactly as the birth API's accepted expression;
    // 1 keV formed as (1e3 eV) is one roundoff below 0.001 MeV.
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
    for (int channel = 0; channel < kChannels; ++channel)
        options.channels[channel] = channels[static_cast<std::size_t>(channel)];
    options.handoff_enabled = 0;
    return options;
}

struct Inputs {
    Grid grid = make_grid();
    Species thermal_number{{0.0, 5.0e19, 5.0e19, 0.0, 0.0, 0.0}};
    double electron_energy_J_m3 = 0.0;
    double ion_energy_J_m3 = 0.0;
    double electron_density_m3 = 1.0e20;
    Species thermal_charge_squared{{1.0, 1.0, 1.0, 4.0, 4.0, 25.0}};
    int inert_count = 1;
    std::array<fusion_inert_ion_v1, 1> inert{{
        fusion_inert_ion_v1{2.0e19, 12.0 * kAtomicMassUnitKg, 36.0}}};
    std::vector<double> coulomb_logs;
    std::vector<double> old_s;
    std::vector<double> old_t;
    std::vector<double> external_birth;
    std::vector<double> escape;
    fusion_coupled_thermal_options_v1 options{};

    Inputs() {
        const double ti = 10.0 * kKeVJ;
        const double te = 5.0 * kKeVJ;
        const double total_ion_number = thermal_number[FUSION_DEUTERON] +
                                        thermal_number[FUSION_TRITON] +
                                        inert[0].density_m3;
        electron_energy_J_m3 = 1.5 * electron_density_m3 * te;
        ion_energy_J_m3 = 1.5 * total_ion_number * ti;
        const std::size_t fast_size =
            static_cast<std::size_t>(kSpecies) * grid.cells();
        coulomb_logs.assign(static_cast<std::size_t>(kSpecies) *
                                (kBaths + inert_count),
                            15.0);
        old_s.assign(fast_size, 0.0);
        old_t.assign(fast_size, 0.0);
        external_birth.assign(fast_size, 0.0);
        escape.assign(fast_size, 0.0);
        options = make_options(Channels{{0, 0, 0, 0, 0}});
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
    trial.s.assign(fast_size, -7.0);
    trial.t.assign(fast_size, -7.0);
    trial.thermal_number.fill(-7.0);
    trial.result = {};
    return trial;
}

int call_trial(const Inputs& inputs, double dt_s, Trial& trial) {
    return fusion_c_coupled_thermal_trial(
        dt_s, &inputs.options, inputs.grid.cells(), inputs.grid.edges.data(),
        inputs.thermal_number.data(), inputs.electron_energy_J_m3,
        inputs.ion_energy_J_m3, inputs.electron_density_m3,
        inputs.thermal_charge_squared.data(), inputs.inert_count,
        inputs.inert_count > 0 ? inputs.inert.data() : nullptr,
        inputs.coulomb_logs.data(), inputs.old_s.data(), inputs.old_t.data(),
        inputs.external_birth.data(), inputs.escape.data(),
        trial.thermal_number.data(), trial.s.data(), trial.t.data(),
        &trial.result);
}

bool zero_ledger(const fusion_source_ledger_v1& ledger) {
    for (double value : ledger.events_m3)
        if (value != 0.0) return false;
    for (double value : ledger.nuclear_born_number_m3)
        if (value != 0.0) return false;
    for (double value : ledger.nuclear_born_energy_J_m3)
        if (value != 0.0) return false;
    for (double value : ledger.external_born_number_m3)
        if (value != 0.0) return false;
    for (double value : ledger.external_born_energy_J_m3)
        if (value != 0.0) return false;
    for (double value : ledger.thermal_consumed_number_m3)
        if (value != 0.0) return false;
    for (double value : ledger.thermal_consumed_energy_J_m3)
        if (value != 0.0) return false;
    for (double value : ledger.fast_consumed_number_m3)
        if (value != 0.0) return false;
    for (double value : ledger.fast_consumed_energy_J_m3)
        if (value != 0.0) return false;
    for (double value : ledger.escaped_number_m3)
        if (value != 0.0) return false;
    for (double value : ledger.escaped_energy_J_m3)
        if (value != 0.0) return false;
    for (double value : ledger.handed_off_number_m3)
        if (value != 0.0) return false;
    for (double value : ledger.handed_off_energy_J_m3)
        if (value != 0.0) return false;
    for (double value : ledger.heat_to_bath_J_m3)
        if (value != 0.0) return false;
    return ledger.neutron_number_m3 == 0.0 &&
           ledger.neutron_energy_J_m3 == 0.0;
}

bool zero_result(const fusion_coupled_thermal_v1& result) {
    if (!zero_ledger(result.ledger)) return false;
    for (double value : result.inert_ion_heat_J_m3)
        if (value != 0.0) return false;
    if (result.electron_energy_J_m3 != 0.0 ||
        result.ion_energy_J_m3 != 0.0 ||
        result.energy_residual_J_m3 != 0.0 ||
        result.max_source_rate_discrepancy != 0.0 ||
        result.max_source_debit_discrepancy != 0.0) {
        return false;
    }
    for (double value : result.particle_residual_m3)
        if (value != 0.0) return false;
    for (double value : result.handoff_L1)
        if (value != 0.0) return false;
    for (double value : result.handoff_mean_error)
        if (value != 0.0) return false;
    for (int value : result.handoff_projected)
        if (value != 0) return false;
    return true;
}

void require_finite_result(const Trial& trial, const std::string& label) {
    for (double value : trial.thermal_number)
        require(std::isfinite(value) && value >= 0.0,
                label + " thermal number is finite and nonnegative");
    for (double value : trial.s)
        require(std::isfinite(value) && value >= 0.0,
                label + " S population is finite and nonnegative");
    for (double value : trial.t)
        require(std::isfinite(value) && value >= 0.0,
                label + " T population is finite and nonnegative");
    const auto& result = trial.result;
    for (double value : result.ledger.events_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " event ledger is finite and nonnegative");
    for (double value : result.ledger.nuclear_born_number_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " nuclear source number is finite");
    for (double value : result.ledger.nuclear_born_energy_J_m3)
        require(std::isfinite(value) && value >= 0.0,
                label + " nuclear source energy is finite");
    require(std::isfinite(result.ledger.neutron_number_m3) &&
                result.ledger.neutron_number_m3 >= 0.0,
            label + " neutron number is finite");
    require(std::isfinite(result.ledger.neutron_energy_J_m3) &&
                result.ledger.neutron_energy_J_m3 >= 0.0,
            label + " neutron energy is finite");
    for (double value : result.ledger.heat_to_bath_J_m3)
        require(std::isfinite(value), label + " network bath heat is finite");
    for (double value : result.inert_ion_heat_J_m3)
        require(std::isfinite(value), label + " inert bath heat is finite");
    require(std::isfinite(result.electron_energy_J_m3) &&
                std::isfinite(result.ion_energy_J_m3) &&
                std::isfinite(result.energy_residual_J_m3),
            label + " thermal energies and residual are finite");
    require(std::isfinite(result.max_source_rate_discrepancy) &&
                std::isfinite(result.max_source_debit_discrepancy),
            label + " source diagnostics are finite");
}

bool same_result(const Trial& left, const Trial& right) {
    if (left.thermal_number != right.thermal_number || left.s != right.s ||
        left.t != right.t) {
        return false;
    }
    const auto& a = left.result;
    const auto& b = right.result;
    const auto same_array = [](const double* x, const double* y,
                               std::size_t count) {
        for (std::size_t i = 0; i < count; ++i)
            if (x[i] != y[i]) return false;
        return true;
    };
    if (!same_array(a.ledger.events_m3, b.ledger.events_m3, kChannels) ||
        !same_array(a.ledger.nuclear_born_number_m3,
                    b.ledger.nuclear_born_number_m3, kSpecies) ||
        !same_array(a.ledger.nuclear_born_energy_J_m3,
                    b.ledger.nuclear_born_energy_J_m3, kSpecies) ||
        !same_array(a.ledger.external_born_number_m3,
                    b.ledger.external_born_number_m3, kSpecies) ||
        !same_array(a.ledger.external_born_energy_J_m3,
                    b.ledger.external_born_energy_J_m3, kSpecies) ||
        !same_array(a.ledger.thermal_consumed_number_m3,
                    b.ledger.thermal_consumed_number_m3, kSpecies) ||
        !same_array(a.ledger.thermal_consumed_energy_J_m3,
                    b.ledger.thermal_consumed_energy_J_m3, kSpecies) ||
        !same_array(a.ledger.fast_consumed_number_m3,
                    b.ledger.fast_consumed_number_m3, kSpecies) ||
        !same_array(a.ledger.fast_consumed_energy_J_m3,
                    b.ledger.fast_consumed_energy_J_m3, kSpecies) ||
        !same_array(a.ledger.escaped_number_m3,
                    b.ledger.escaped_number_m3, kSpecies) ||
        !same_array(a.ledger.escaped_energy_J_m3,
                    b.ledger.escaped_energy_J_m3, kSpecies) ||
        !same_array(a.ledger.handed_off_number_m3,
                    b.ledger.handed_off_number_m3, kSpecies) ||
        !same_array(a.ledger.handed_off_energy_J_m3,
                    b.ledger.handed_off_energy_J_m3, kSpecies) ||
        !same_array(a.ledger.heat_to_bath_J_m3,
                    b.ledger.heat_to_bath_J_m3, kSpecies * kBaths) ||
        a.ledger.neutron_number_m3 != b.ledger.neutron_number_m3 ||
        a.ledger.neutron_energy_J_m3 != b.ledger.neutron_energy_J_m3 ||
        !same_array(a.inert_ion_heat_J_m3, b.inert_ion_heat_J_m3, kSpecies) ||
        a.electron_energy_J_m3 != b.electron_energy_J_m3 ||
        a.ion_energy_J_m3 != b.ion_energy_J_m3 ||
        !same_array(a.particle_residual_m3, b.particle_residual_m3,
                    kSpecies) ||
        a.energy_residual_J_m3 != b.energy_residual_J_m3 ||
        a.max_source_rate_discrepancy != b.max_source_rate_discrepancy ||
        a.max_source_debit_discrepancy != b.max_source_debit_discrepancy ||
        !same_array(a.handoff_L1, b.handoff_L1, kSpecies) ||
        !same_array(a.handoff_mean_error, b.handoff_mean_error, kSpecies)) {
        return false;
    }
    for (int i = 0; i < kSpecies; ++i)
        if (a.handoff_projected[i] != b.handoff_projected[i]) return false;
    return true;
}

double fast_number(const Trial& trial, int species) {
    const int cells = static_cast<int>(trial.s.size() / kSpecies);
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

double sum_heat(const fusion_source_ledger_v1& ledger, int species,
                int first_bath, int last_bath) {
    double sum = 0.0;
    for (int bath = first_bath; bath < last_bath; ++bath)
        sum += ledger.heat_to_bath_J_m3[species * kBaths + bath];
    return sum;
}

double total_array(const double* values, int count) {
    double sum = 0.0;
    for (int i = 0; i < count; ++i) sum += values[i];
    return sum;
}

void require_outputs_cleared(const Trial& trial, const std::string& label) {
    for (double value : trial.thermal_number)
        require(value == 0.0, label + " clears thermal-number output");
    for (double value : trial.s)
        require(value == 0.0, label + " clears S output");
    for (double value : trial.t)
        require(value == 0.0, label + " clears T output");
    require(zero_result(trial.result), label + " clears result output");
}

void test_no_reaction_no_fast_is_identity() {
    Inputs inputs;
    const Species old_thermal = inputs.thermal_number;
    const std::vector<double> old_s = inputs.old_s;
    const std::vector<double> old_t = inputs.old_t;
    const double old_electron = inputs.electron_energy_J_m3;
    const double old_ion = inputs.ion_energy_J_m3;

    Trial trial = make_trial(inputs);
    const int status = call_trial(inputs, 1.0e-4, trial);
    require(status == PB11_STATUS_OK,
            "no-reaction/no-fast trial returns OK (status=" +
                std::to_string(status) + ")");
    require(trial.thermal_number == old_thermal,
            "no-reaction trial preserves thermal numbers");
    require(trial.s == old_s && trial.t == old_t,
            "no-reaction trial preserves empty fast populations");
    require(trial.result.electron_energy_J_m3 == old_electron &&
                trial.result.ion_energy_J_m3 == old_ion,
            "no-reaction trial preserves electron and ion energy");
    require(zero_ledger(trial.result.ledger),
            "no-reaction trial has zero source ledger");
    for (double value : trial.result.inert_ion_heat_J_m3)
        require(value == 0.0, "no-fast trial has zero inert-carbon heat");
    require(trial.result.energy_residual_J_m3 == 0.0,
            "no-reaction trial has zero energy residual");
    for (double value : trial.result.particle_residual_m3)
        require(value == 0.0, "no-reaction trial has zero particle residual");
    require(inputs.thermal_number == old_thermal && inputs.old_s == old_s &&
                inputs.old_t == old_t &&
                inputs.electron_energy_J_m3 == old_electron &&
                inputs.ion_energy_J_m3 == old_ion,
            "no-reaction trial leaves inputs unchanged");
}

void test_dt_burn_slowing_accounting_and_repeatability() {
    Inputs inputs;
    inputs.options = make_options(Channels{{0, 0, 0, 1, 0}});
    const int alpha_cell = nearest_cell(inputs.grid, 3.5 * kMeVJ);
    const std::size_t alpha_index =
        static_cast<std::size_t>(FUSION_HELIUM4) * inputs.grid.cells() +
        static_cast<std::size_t>(alpha_cell);
    inputs.old_s[alpha_index] = 1.0e12;

    const Species old_thermal = inputs.thermal_number;
    const std::vector<double> old_s = inputs.old_s;
    const std::vector<double> old_t = inputs.old_t;
    const std::vector<double> old_external = inputs.external_birth;
    const std::vector<double> old_escape = inputs.escape;
    const double old_electron = inputs.electron_energy_J_m3;
    const double old_ion = inputs.ion_energy_J_m3;

    Trial first = make_trial(inputs);
    const double dt = 1.0e-4;
    require(call_trial(inputs, dt, first) == PB11_STATUS_OK,
            "DT burn and alpha slowing trial returns OK");
    require_finite_result(first, "DT burn and alpha slowing");
    require(first.result.ledger.events_m3[FUSION_DT_ALPHAN] > 0.0,
            "DT event source is nonzero");
    for (int channel = 0; channel < kChannels; ++channel)
        if (channel != FUSION_DT_ALPHAN)
            require(first.result.ledger.events_m3[channel] == 0.0,
                    "disabled reaction channels stay zero");

    const auto& ledger = first.result.ledger;
    require(ledger.nuclear_born_number_m3[FUSION_HELIUM4] > 0.0 &&
                ledger.nuclear_born_energy_J_m3[FUSION_HELIUM4] > 0.0,
            "DT produces charged alpha source number and energy");
    require(ledger.neutron_number_m3 > 0.0 &&
                ledger.neutron_energy_J_m3 > 0.0,
            "DT produces explicit neutron number and energy");
    require(ledger.thermal_consumed_number_m3[FUSION_DEUTERON] > 0.0 &&
                ledger.thermal_consumed_number_m3[FUSION_TRITON] > 0.0 &&
                ledger.thermal_consumed_energy_J_m3[FUSION_DEUTERON] > 0.0 &&
                ledger.thermal_consumed_energy_J_m3[FUSION_TRITON] > 0.0,
            "DT thermal fuel source ledger is nonzero");
    require(first.thermal_number[FUSION_DEUTERON] <
                old_thermal[FUSION_DEUTERON] &&
                first.thermal_number[FUSION_TRITON] <
                    old_thermal[FUSION_TRITON],
            "DT trial depletes both thermal reactants");
    for (int species : {FUSION_PROTON, FUSION_HELIUM3, FUSION_HELIUM4,
                        FUSION_BORON11})
        require(first.thermal_number[species] == old_thermal[species],
                "DT trial leaves passive thermal species unchanged");

    const double initial_fast_energy =
        fast_energy(inputs.grid, inputs.old_s, inputs.old_t, FUSION_HELIUM4);
    const double final_fast_energy =
        fast_energy(inputs.grid, first.s, first.t, FUSION_HELIUM4);
    require(fast_number(first, FUSION_HELIUM4) >
                inputs.old_s[alpha_index],
            "DT trial retains seeded alpha and adds the birth population");
    require(final_fast_energy > initial_fast_energy,
            "DT alpha birth leaves positive fast alpha energy");
    require(first.result.inert_ion_heat_J_m3[FUSION_HELIUM4] > 0.0,
            "DT alpha slowing deposits distinct inert-carbon heat");

    // Explicit charged/neutron source stoichiometry and source energy budget.
    const double events = ledger.events_m3[FUSION_DT_ALPHAN];
    require(close_scaled(ledger.thermal_consumed_number_m3[FUSION_DEUTERON],
                         events, 3.0e-10, 1.0e-30),
            "DT consumes one deuteron per event");
    require(close_scaled(ledger.thermal_consumed_number_m3[FUSION_TRITON],
                         events, 3.0e-10, 1.0e-30),
            "DT consumes one triton per event");
    require(close_scaled(ledger.nuclear_born_number_m3[FUSION_HELIUM4],
                         events, 3.0e-10, 1.0e-30),
            "DT births one alpha per event");
    require(close_scaled(ledger.neutron_number_m3, events, 3.0e-10,
                         1.0e-30),
            "DT births one explicit neutron per event");
    fusion_nuclear_channel_v1 dt_channel{};
    require(fusion_c_nuclear_channel(FUSION_DT_ALPHAN, &dt_channel) ==
                PB11_STATUS_OK,
            "DT channel metadata is available for accounting");
    const double q_energy = events * dt_channel.q_J;
    const double source_born_energy =
        total_array(ledger.nuclear_born_energy_J_m3, kSpecies) +
        ledger.neutron_energy_J_m3;
    const double source_removed_energy =
        total_array(ledger.thermal_consumed_energy_J_m3, kSpecies) +
        total_array(ledger.fast_consumed_energy_J_m3, kSpecies);
    require(close_scaled(source_born_energy,
                         source_removed_energy + q_energy, 3.0e-9, 1.0e-20),
            "DT charged plus neutron source energy closes Q budget");

    // Every kinetic species has an independent inventory identity.  The
    // inert-ion account is deliberately included separately from the seven
    // source-state bath columns.
    for (int species = 0; species < kSpecies; ++species) {
        const double old_number =
            [&]() {
                double value = 0.0;
                const int cells = inputs.grid.cells();
                for (int cell = 0; cell < cells; ++cell) {
                    const std::size_t index =
                        static_cast<std::size_t>(species) * cells + cell;
                    value += inputs.old_s[index] + inputs.old_t[index];
                }
                return value;
            }();
        const double new_number = fast_number(first, species);
        const double expected_number =
            old_number + ledger.nuclear_born_number_m3[species] +
            ledger.external_born_number_m3[species] -
            ledger.fast_consumed_number_m3[species] -
            ledger.escaped_number_m3[species] -
            ledger.handed_off_number_m3[species];
        require(close_scaled(new_number, expected_number, 3.0e-9, 1.0e-20),
                "fast particle source/consumption/escape identity closes");

        const double old_energy =
            fast_energy(inputs.grid, inputs.old_s, inputs.old_t, species);
        const double new_energy =
            fast_energy(inputs.grid, first.s, first.t, species);
        const double network_heat = sum_heat(ledger, species, 0, kBaths);
        const double expected_energy =
            old_energy + ledger.nuclear_born_energy_J_m3[species] +
            ledger.external_born_energy_J_m3[species] -
            ledger.fast_consumed_energy_J_m3[species] -
            ledger.escaped_energy_J_m3[species] -
            ledger.handed_off_energy_J_m3[species] - network_heat -
            first.result.inert_ion_heat_J_m3[species];
        require(close_scaled(new_energy, expected_energy, 3.0e-8, 1.0e-20),
                "fast energy and separate inert heat identity closes");
    }

    const double electron_heat =
        [&]() {
            double value = 0.0;
            for (int species = 0; species < kSpecies; ++species)
                value += ledger.heat_to_bath_J_m3[species * kBaths];
            return value;
        }();
    double network_ion_heat = 0.0;
    for (int species = 0; species < kSpecies; ++species)
        network_ion_heat += sum_heat(ledger, species, 1, kBaths);
    const double inert_heat =
        total_array(first.result.inert_ion_heat_J_m3, kSpecies);
    const double thermal_debit =
        total_array(ledger.thermal_consumed_energy_J_m3, kSpecies);
    require(close_scaled(first.result.electron_energy_J_m3 - old_electron,
                         electron_heat, 3.0e-8, 1.0e-20),
            "electron energy change equals electron heat ledger");
    require(close_scaled(first.result.ion_energy_J_m3 - old_ion,
                         -thermal_debit + network_ion_heat + inert_heat,
                         3.0e-8, 1.0e-20),
            "ion energy change includes nuclear debit and both ion heat ledgers");

    double escaped_energy = 0.0;
    for (double value : ledger.escaped_energy_J_m3) escaped_energy += value;
    const double total_energy_change =
        first.result.electron_energy_J_m3 - old_electron +
        first.result.ion_energy_J_m3 - old_ion +
        (final_fast_energy - initial_fast_energy) +
        ledger.neutron_energy_J_m3 + escaped_energy;
    require(close_scaled(total_energy_change, q_energy, 3.0e-8, 1.0e-18),
            "DT total thermal/fast/neutron energy budget closes");
    require(std::abs(first.result.energy_residual_J_m3) <=
                1.0e-8 * std::max(1.0, std::abs(q_energy)),
            "DT coupled energy residual is small");

    Trial second = make_trial(inputs);
    require(call_trial(inputs, dt, second) == PB11_STATUS_OK,
            "repeated DT burn and slowing trial returns OK");
    require(same_result(first, second),
            "repeated DT trial from unchanged inputs is identical");
    require(inputs.thermal_number == old_thermal && inputs.old_s == old_s &&
                inputs.old_t == old_t && inputs.external_birth == old_external &&
                inputs.escape == old_escape &&
                inputs.electron_energy_J_m3 == old_electron &&
                inputs.ion_energy_J_m3 == old_ion,
            "DT trial leaves all original inputs unchanged");
}

void test_inert_heat_is_separate_and_closes() {
    Inputs inputs;
    inputs.options = make_options(Channels{{0, 0, 0, 0, 0}});
    const int alpha_cell = nearest_cell(inputs.grid, 3.5 * kMeVJ);
    inputs.old_s[static_cast<std::size_t>(FUSION_HELIUM4) * inputs.grid.cells() +
                 static_cast<std::size_t>(alpha_cell)] =
        1.0e12;
    const std::vector<double> old_s = inputs.old_s;
    const double old_electron = inputs.electron_energy_J_m3;
    const double old_ion = inputs.ion_energy_J_m3;

    Trial trial = make_trial(inputs);
    const int status = call_trial(inputs, 1.0e-5, trial);
    require(status == PB11_STATUS_OK,
            "inert-carbon alpha slowing trial returns OK (status=" +
                std::to_string(status) + ")");
    require_finite_result(trial, "inert-carbon alpha slowing");
    require(zero_ledger(trial.result.ledger) == false,
            "collision-only trial has a nonzero seven-bath heat ledger");
    require(trial.result.inert_ion_heat_J_m3[FUSION_HELIUM4] > 0.0,
            "collision-only trial records nonzero alpha-to-carbon heat");
    const double inert_heat =
        total_array(trial.result.inert_ion_heat_J_m3, kSpecies);
    require(inert_heat > 0.0, "inert heat total is positive");
    double network_ion_heat = 0.0;
    for (int species = 0; species < kSpecies; ++species)
        network_ion_heat += sum_heat(trial.result.ledger, species, 1, kBaths);
    require(close_scaled(trial.result.ion_energy_J_m3 - old_ion -
                             network_ion_heat,
                         inert_heat, 3.0e-8, 1.0e-10),
            "ion energy identifies inert heat outside source-state columns");
    double electron_heat = 0.0;
    for (int species = 0; species < kSpecies; ++species)
        electron_heat += sum_heat(trial.result.ledger, species, 0, 1);
    require(close_scaled(trial.result.electron_energy_J_m3 - old_electron,
                         electron_heat, 3.0e-8,
                         1.0e-10),
            "collision-only electron energy matches electron heat ledger");
    const double initial_fast =
        fast_energy(inputs.grid, inputs.old_s, inputs.old_t, FUSION_HELIUM4);
    const double final_fast =
        fast_energy(inputs.grid, trial.s, trial.t, FUSION_HELIUM4);
    const double total_change = trial.result.electron_energy_J_m3 - old_electron +
                                trial.result.ion_energy_J_m3 - old_ion +
                                final_fast - initial_fast;
    require(close_scaled(total_change, 0.0, 3.0e-8, 1.0e-10),
            "collision-only thermal and fast energy budget closes");
    require(inputs.old_s == old_s &&
                inputs.electron_energy_J_m3 == old_electron &&
                inputs.ion_energy_J_m3 == old_ion,
            "collision-only trial leaves original inputs unchanged");
}

void test_mixed_pool_handoff_uses_self_consistent_target() {
    Inputs inputs;
    inputs.options = make_options(Channels{{0, 0, 0, 0, 0}});
    inputs.options.handoff_enabled = 1;
    inputs.options.handoff_max_L1 = 0.02;
    inputs.options.handoff_max_mean_error = 0.02;
    const double ti = 10.0 * kKeVJ;
    inputs.electron_energy_J_m3 = 1.5 * inputs.electron_density_m3 * ti;

    std::vector<double> probability(static_cast<std::size_t>(inputs.grid.cells()));
    fusion_maxwellian_grid_v1 grid_result{};
    require(fusion_c_maxwellian_energy_grid(
                inputs.grid.cells(), 1.01 * ti, inputs.grid.edges.data(),
                probability.data(), &grid_result) == PB11_STATUS_OK,
            "handoff source Maxwellian grid returns OK");
    const double thermal_pool =
        inputs.thermal_number[FUSION_DEUTERON] +
        inputs.thermal_number[FUSION_TRITON] + inputs.inert[0].density_m3;
    const double fast_alpha_number = 0.01 * thermal_pool;
    for (int cell = 0; cell < inputs.grid.cells(); ++cell) {
        require(std::isfinite(probability[static_cast<std::size_t>(cell)]) &&
                    probability[static_cast<std::size_t>(cell)] >= 0.0,
                "handoff source Maxwellian probabilities are finite");
        inputs.old_t[static_cast<std::size_t>(FUSION_HELIUM4) *
                         inputs.grid.cells() + static_cast<std::size_t>(cell)] =
            fast_alpha_number * probability[static_cast<std::size_t>(cell)];
    }

    Trial trial = make_trial(inputs);
    require(call_trial(inputs, 1.0e-12, trial) == PB11_STATUS_OK,
            "mixed-pool handoff trial returns OK");
    require_finite_result(trial, "mixed-pool handoff");
    require(trial.result.handoff_projected[FUSION_HELIUM4] == 1,
            "mixed-pool alpha candidate is projected");
    for (int species = 0; species < kSpecies; ++species)
        if (species != FUSION_HELIUM4)
            require(trial.result.handoff_projected[species] == 0,
                    "only seeded alpha candidate is projected");
    const double handed_number =
        trial.result.ledger.handed_off_number_m3[FUSION_HELIUM4];
    const double handed_energy =
        trial.result.ledger.handed_off_energy_J_m3[FUSION_HELIUM4];
    require(handed_number > 0.0 && handed_energy > 0.0,
            "mixed-pool handoff returns positive alpha fluid amount");

    fusion_thermal_increment_v1 increment{};
    require(fusion_c_coupled_thermal_increment(&trial.result.ledger,
                trial.result.inert_ion_heat_J_m3, &increment) == PB11_STATUS_OK,
            "actual projected handoff increment");
    require(increment.thermal_number_m3[FUSION_HELIUM4] == handed_number,
            "handoff fluid number counted once");
    require(close_scaled(inputs.ion_energy_J_m3 + increment.ion_energy_J_m3,
                         trial.result.ion_energy_J_m3, 4e-15, 0.0),
            "handoff and signed correction reconstruct ion energy once");

    double final_thermal_number = 0.0;
    for (double value : trial.thermal_number) final_thermal_number += value;
    const double final_common_ti =
        trial.result.ion_energy_J_m3 /
        (1.5 * (final_thermal_number + inputs.inert[0].density_m3));
    const double handed_ti = handed_energy / (1.5 * handed_number);
    require(close_scaled(final_common_ti, handed_ti, 4.0e-10, 1.0e-30),
            "projected alpha fluid energy uses the final common ion temperature");
    require(final_thermal_number >
                inputs.thermal_number[FUSION_DEUTERON] +
                    inputs.thermal_number[FUSION_TRITON],
            "mixed-pool handoff adds alpha to the thermal population");
    require(trial.result.ledger.events_m3[FUSION_DT_ALPHAN] == 0.0 &&
                zero_ledger(trial.result.ledger) == false,
            "handoff trial has transfer ledger without nuclear reactions");
}

void test_negative_thermal_energy_and_invalid_inputs_reject_and_clear() {
    Inputs extreme;
    extreme.options = make_options(Channels{{0, 0, 0, 0, 0}});
    extreme.thermal_number = Species{{0.0, 1.0e20, 1.0e20, 0.0, 0.0, 0.0}};
    extreme.electron_density_m3 = 2.0e20;
    extreme.electron_energy_J_m3 = 1.5 * extreme.electron_density_m3 *
                                   (1.0e-3 * kElectronVoltJ);
    extreme.ion_energy_J_m3 = 1.5 *
                              (extreme.thermal_number[FUSION_DEUTERON] +
                               extreme.thermal_number[FUSION_TRITON] +
                               extreme.inert[0].density_m3) *
                              (1.0e-3 * kElectronVoltJ);
    const int alpha_cell = nearest_cell(extreme.grid, 3.5 * kMeVJ);
    extreme.old_s[static_cast<std::size_t>(FUSION_HELIUM4) *
                      extreme.grid.cells() +
                  static_cast<std::size_t>(alpha_cell)] =
        1.0e20;
    Trial negative_energy = make_trial(extreme);
    require(call_trial(extreme, 1.0e3, negative_energy) != PB11_STATUS_OK,
            "extreme collision step rejects negative thermal energy");
    require_outputs_cleared(negative_energy,
                            "negative-thermal-energy rejection");

    Inputs invalid = Inputs{};
    invalid.options.birth.relative_order = 3;
    Trial bad_options = make_trial(invalid);
    require(call_trial(invalid, 1.0e-4, bad_options) != PB11_STATUS_OK,
            "invalid birth quadrature option is rejected");
    require_outputs_cleared(bad_options, "invalid options rejection");

    invalid = Inputs{};
    invalid.coulomb_logs[0] = std::numeric_limits<double>::quiet_NaN();
    Trial bad_log = make_trial(invalid);
    require(call_trial(invalid, 1.0e-4, bad_log) != PB11_STATUS_OK,
            "nonfinite Coulomb log is rejected");
    require_outputs_cleared(bad_log, "invalid Coulomb-log rejection");

    invalid = Inputs{};
    Trial bad_array = make_trial(invalid);
    require(fusion_c_coupled_thermal_trial(
                1.0e-4, &invalid.options, invalid.grid.cells(),
                invalid.grid.edges.data(), invalid.thermal_number.data(),
                invalid.electron_energy_J_m3, invalid.ion_energy_J_m3,
                invalid.electron_density_m3,
                invalid.thermal_charge_squared.data(), invalid.inert_count,
                invalid.inert.data(), invalid.coulomb_logs.data(),
                invalid.old_s.data(), invalid.old_t.data(), nullptr,
                invalid.escape.data(), bad_array.thermal_number.data(),
                bad_array.s.data(), bad_array.t.data(), &bad_array.result) !=
                PB11_STATUS_OK,
            "null external source array is rejected");
    require_outputs_cleared(bad_array, "invalid array rejection");
}

void test_clipped_charged_birth_grid_rejects_and_clears() {
    Inputs inputs;
    inputs.options = make_options(Channels{{0, 0, 0, 1, 0}});
    const int alpha_cell = nearest_cell(inputs.grid, 3.5 * kMeVJ);
    inputs.old_s[static_cast<std::size_t>(FUSION_HELIUM4) * inputs.grid.cells() +
                 static_cast<std::size_t>(alpha_cell)] =
        1.0e12;
    // The coupled API uses the caller grid for the physical birth source.  A
    // deliberately clipped upper bound must reject the whole trial instead
    // of silently classifying charged source spill as escape.
    const int split_edge = 1 + 96;
    const int last_edge = inputs.grid.cells();
    for (int edge = split_edge + 1; edge <= last_edge; ++edge) {
        const double fraction = static_cast<double>(edge - split_edge) /
                                (last_edge - split_edge);
        inputs.grid.edges[static_cast<std::size_t>(edge)] =
            0.5 * kMeVJ + fraction * (2.0 * kMeVJ - 0.5 * kMeVJ);
    }
    Trial clipped = make_trial(inputs);
    require(call_trial(inputs, 1.0e-4, clipped) != PB11_STATUS_OK,
            "charged DT birth outside grid rejects whole trial");
    require_outputs_cleared(clipped, "clipped charged-birth rejection");
}

}  // namespace

int main() {
    try {
        test_no_reaction_no_fast_is_identity();
        test_dt_burn_slowing_accounting_and_repeatability();
        test_inert_heat_is_separate_and_closes();
        test_mixed_pool_handoff_uses_self_consistent_target();
        test_negative_thermal_energy_and_invalid_inputs_reject_and_clear();
        test_clipped_charged_birth_grid_rejects_and_clears();
        std::cout << "PASS: coupled thermal identity, DT source/energy closure, "
                     "inert-ion heat, rejection and repeatability tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
