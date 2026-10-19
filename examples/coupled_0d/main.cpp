#include "fusion_coupled_thermal.h"
#include "fusion_beam_birth_table.h"
#include "fusion_handoff.h"
#include "fusion_network.h"
#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"

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


// This tutorial explicitly assumes every network ion and carbon is fully
// stripped. A charge mixture would require caller-supplied mean Z, which
// cannot generally be recovered from its mean Z squared collision factor.
long double positive_charge(const Inputs& inputs, const Species& thermal,
                            const std::vector<double>& s,
                            const std::vector<double>& t) {
    constexpr double carbon_mean_charge = 6.0;
    require(inputs.inert_count == 1, "tutorial has one explicit carbon bath");
    long double charge = carbon_mean_charge * inputs.inert[0].density_m3;
    for (int species = 0; species < kSpecies; ++species) {
        fusion_nuclear_mass_v1 mass{};
        require(fusion_c_nuclear_mass(species, &mass) == PB11_STATUS_OK,
                "public nuclear charge status");
        long double number = thermal[species];
        for (int cell = 0; cell < inputs.grid.cells(); ++cell) {
            const auto i = static_cast<std::size_t>(species) * inputs.grid.cells() + cell;
            number += static_cast<long double>(s[i]) + t[i];
        }
        charge += mass.nuclear_charge * number;
    }
    require(std::isfinite(charge) && charge > 0, "positive finite charge inventory");
    return charge;
}

// A trial owns every candidate buffer. Nothing below is accepted until commit.
void run(const std::string& fuel) {
    Inputs accepted;
    const bool pb = fuel == "pb";
    const int channel = pb ? FUSION_PB11_3ALPHA : FUSION_DT_ALPHAN;
    accepted.options.channels[channel] = 1;
    accepted.options.handoff_enabled = 1;
    if (pb) {
        // Use the existing study_coupled_thermal pB cutoff and lower bound.
        // The public pB source requires parent energy at most 12 MeV.
        accepted.options.birth.relative_max_J = 2.5 * kMeVJ;
        const double first_edge = 1.0e-18 * kKeVJ;
        for (int j = 1; j <= 97; ++j)
            accepted.grid.edges[j] = std::exp(std::log(first_edge) +
                double(j - 1) / 96 * (std::log(0.5 * kMeVJ) - std::log(first_edge)));
        // Same explicit pB initial-temperature choice as study_coupled_thermal.
        accepted.thermal_number = Species{{5.0e19, 0, 0, 0, 0, 5.0e19}};
        accepted.ion_energy_J_m3 = 1.5 * 1.2e20 * 100.0 * kKeVJ;
    }
    // A seeded hot alpha makes the collision heat exchange visible immediately.
    accepted.old_s[static_cast<std::size_t>(FUSION_HELIUM4) *
        accepted.grid.cells() + nearest_cell(accepted.grid, 3.5 * kMeVJ)] = 1.0e12;
    accepted.electron_density_m3 = static_cast<double>(positive_charge(
        accepted, accepted.thermal_number, accepted.old_s, accepted.old_t));
    accepted.electron_energy_J_m3 = 1.5 * accepted.electron_density_m3 * 5.0 * kKeVJ;
    constexpr double dt = 1.0e-4; // Existing DT integration fixture step.
    fusion_nuclear_channel_v1 metadata{};
    require(fusion_c_nuclear_channel(channel, &metadata) == PB11_STATUS_OK,
            "channel metadata status");
    for (int step = 1; step <= 2; ++step) {
        const long double old_charge = positive_charge(
            accepted, accepted.thermal_number, accepted.old_s, accepted.old_t);
        // Caller updates charge balance between steps; Ue remains accepted
        // energy. Only electron density is frozen inside the local trial.
        accepted.electron_density_m3 = static_cast<double>(old_charge);
        const Inputs before = accepted;
        Trial discarded = make_trial(accepted);
        int status = call_trial(accepted, dt, discarded);
        require(status == PB11_STATUS_OK,
                fuel + " discarded trial status=" + std::to_string(status));
        require(accepted.old_s == before.old_s && accepted.old_t == before.old_t &&
                    accepted.thermal_number == before.thermal_number &&
                    accepted.electron_energy_J_m3 == before.electron_energy_J_m3 &&
                    accepted.ion_energy_J_m3 == before.ion_energy_J_m3,
                "trial must leave accepted state untouched");
        Trial candidate = make_trial(accepted);
        status = call_trial(accepted, dt, candidate);
        require(status == PB11_STATUS_OK,
                fuel + " retry status=" + std::to_string(status));
        require(same_result(discarded, candidate), "discard/retry is identical");
        require_finite_result(candidate, fuel);
        require(candidate.result.electron_energy_J_m3 >= 0 &&
                    candidate.result.ion_energy_J_m3 >= 0,
                "nonnegative thermal energies");
        const long double new_charge = positive_charge(
            accepted, candidate.thermal_number, candidate.s, candidate.t);
        const long double charge_relative_residual = std::abs(new_charge - old_charge) /
            std::max(std::abs(new_charge), std::abs(old_charge));
        require(charge_relative_residual <= 1.0e-12L,
                "total positive charge conserved across trial");
        const auto& r = candidate.result;
        const auto& l = r.ledger;
        require(l.events_m3[channel] > 0 &&
                    l.nuclear_born_number_m3[FUSION_HELIUM4] > 0,
                "actual fusion birth required");
        require(r.max_source_rate_discrepancy <= accepted.options.max_source_rate_error &&
                    r.max_source_debit_discrepancy <= accepted.options.max_source_debit_error,
                "source discrepancy gates");
        fusion_thermal_increment_v1 increment{};
        require(fusion_c_coupled_thermal_increment(&l, r.inert_ion_heat_J_m3,
                    &increment) == PB11_STATUS_OK, "thermal increment status");
        double old_fast = 0, new_fast = 0, max_particle_residual = 0;
        for (int species = 0; species < kSpecies; ++species) {
            const double old_u = fast_energy(accepted.grid, accepted.old_s,
                                             accepted.old_t, species);
            const double new_u = fast_energy(accepted.grid, candidate.s,
                                             candidate.t, species);
            old_fast += old_u; new_fast += new_u;
            double old_n = 0;
            for (int cell = 0; cell < accepted.grid.cells(); ++cell) {
                const auto i = static_cast<std::size_t>(species) * accepted.grid.cells() + cell;
                old_n += accepted.old_s[i] + accepted.old_t[i];
            }
            const double expected_n = old_n + l.nuclear_born_number_m3[species] +
                l.external_born_number_m3[species] - l.fast_consumed_number_m3[species] -
                l.escaped_number_m3[species] - l.handed_off_number_m3[species];
            require(close_scaled(fast_number(candidate, species), expected_n,
                                 3.0e-9, 1.0e-20), "fast number identity");
            const double expected_u = old_u + l.nuclear_born_energy_J_m3[species] +
                l.external_born_energy_J_m3[species] - l.fast_consumed_energy_J_m3[species] -
                l.escaped_energy_J_m3[species] - l.handed_off_energy_J_m3[species] -
                sum_heat(l, species, 0, kBaths) - r.inert_ion_heat_J_m3[species];
            require(close_scaled(new_u, expected_u, 3.0e-8, 1.0e-20),
                    "fast energy identity including separate inert heat");
            require(std::isfinite(r.particle_residual_m3[species]), "finite particle residual");
            max_particle_residual = std::max(max_particle_residual,
                                             std::abs(r.particle_residual_m3[species]));
        }
        const double q = l.events_m3[channel] * metadata.q_J;
        const double born = total_array(l.nuclear_born_energy_J_m3, kSpecies) + l.neutron_energy_J_m3;
        const double removed = total_array(l.thermal_consumed_energy_J_m3, kSpecies) +
                               total_array(l.fast_consumed_energy_J_m3, kSpecies);
        require(close_scaled(born, removed + q, 3.0e-9, 1.0e-20), "source Q identity");
        const double delta = increment.electron_energy_J_m3 + increment.ion_energy_J_m3 +
            new_fast - old_fast + l.neutron_energy_J_m3 +
            total_array(l.escaped_energy_J_m3, kSpecies);
        require(close_scaled(delta, q, 3.0e-8, 1.0e-18), "total energy identity");
        require(std::abs(r.energy_residual_J_m3) <= 1.0e-8 * std::max(1.0, std::abs(q)),
                "coupled energy residual");
        // Commit all thermal and both kinetic inventories together. Returned
        // Ue/Ui are authoritative; explicit increments preserve weak source
        // diagnostics even when subtracting two large inventories loses digits.
        accepted.thermal_number = candidate.thermal_number;
        accepted.old_s.swap(candidate.s); accepted.old_t.swap(candidate.t);
        accepted.electron_energy_J_m3 = r.electron_energy_J_m3;
        accepted.ion_energy_J_m3 = r.ion_energy_J_m3;
        std::cout << fuel << ',' << step << ',' << step * dt;
        for (double n : accepted.thermal_number) std::cout << ',' << n;
        std::cout << ',' << accepted.electron_energy_J_m3 << ',' << accepted.ion_energy_J_m3
                  << ',' << new_fast << ',' << l.events_m3[channel]
                  << ',' << l.handed_off_number_m3[FUSION_HELIUM4]
                  << ',' << increment.electron_energy_J_m3 << ',' << increment.ion_energy_J_m3
                  << ',' << r.energy_residual_J_m3 << ',' << max_particle_residual
                  << ',' << accepted.electron_density_m3 << ',' << charge_relative_residual << '\n';
    }
}
} // namespace
int main() {
    try {
        std::cout.precision(17);
        std::cerr << "kernel=" << fusion_c_beam_birth_table_kernel_identity() << "\n";
        std::cout << "fuel,step,time_s,nH,nD,nT,nHe3,nHe4,nB11,Ue_J_m3,Ui_J_m3,fast_U_J_m3,events_m3,ash_handoff_m3,dUe_J_m3,dUi_J_m3,energy_residual_J_m3,max_particle_residual_m3,ne_frozen_m3,charge_relative_residual\n";
        run("dt"); run("pb");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; return 1;
    }
}
