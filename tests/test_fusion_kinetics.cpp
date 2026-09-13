#include "fusion_kinetics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr double kJoulesPerKeV = 1.602176634e-16;
constexpr double kDiffusion = 1.0e-30;

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool close_relative(double actual, double expected, double tolerance) {
    if (actual == expected) {
        return true;
    }
    if (!std::isfinite(actual) || !std::isfinite(expected)) {
        return false;
    }
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return scale > 0.0 && std::abs(actual - expected) <= tolerance * scale;
}

bool all_zero(const fusion_kinetic_ledger_v1& ledger) {
    return ledger.initial_number_m3 == 0.0 &&
           ledger.final_number_m3 == 0.0 &&
           ledger.initial_energy_J_m3 == 0.0 &&
           ledger.final_energy_J_m3 == 0.0 &&
           ledger.born_number_m3 == 0.0 &&
           ledger.born_energy_J_m3 == 0.0 &&
           ledger.escaped_number_m3 == 0.0 &&
           ledger.escaped_energy_J_m3 == 0.0 &&
           ledger.thermalized_number_m3 == 0.0 &&
           ledger.thermalized_energy_J_m3 == 0.0 &&
           ledger.particle_balance_error_m3 == 0.0 &&
           ledger.energy_balance_error_J_m3 == 0.0;
}

void fill_ledger(fusion_kinetic_ledger_v1* ledger, double value) {
    ledger->initial_number_m3 = value;
    ledger->final_number_m3 = value;
    ledger->initial_energy_J_m3 = value;
    ledger->final_energy_J_m3 = value;
    ledger->born_number_m3 = value;
    ledger->born_energy_J_m3 = value;
    ledger->escaped_number_m3 = value;
    ledger->escaped_energy_J_m3 = value;
    ledger->thermalized_number_m3 = value;
    ledger->thermalized_energy_J_m3 = value;
    ledger->particle_balance_error_m3 = value;
    ledger->energy_balance_error_J_m3 = value;
}

void check_finite_nonnegative(const double* values, int count,
                              const std::string& label) {
    for (int i = 0; i < count; ++i) {
        check(std::isfinite(values[i]) && values[i] >= 0.0,
              label + " is finite and nonnegative");
    }
}

void check_ledger_balance(const fusion_kinetic_ledger_v1& ledger,
                          const std::string& label,
                          double net_heat_to_baths = 0.0) {
    const double particle_rhs = ledger.initial_number_m3 +
                                 ledger.born_number_m3 -
                                 ledger.escaped_number_m3 -
                                 ledger.thermalized_number_m3;
    const double energy_rhs = ledger.initial_energy_J_m3 +
                              ledger.born_energy_J_m3 -
                              ledger.escaped_energy_J_m3 -
                              ledger.thermalized_energy_J_m3 -
                              net_heat_to_baths;
    check(close_relative(ledger.final_number_m3, particle_rhs, 1.0e-10),
          label + " particle ledger closes");
    check(close_relative(ledger.final_energy_J_m3, energy_rhs, 1.0e-10),
          label + " energy ledger closes");
    check(std::isfinite(ledger.particle_balance_error_m3) &&
              std::isfinite(ledger.energy_balance_error_J_m3),
          label + " ledger errors are finite");
}

void test_one_cell_analytic_step() {
    const double edges[2]{1.0e-16, 3.0e-16};
    const double old[1]{4.0e20};
    const double birth[1]{3.0e20};
    const double escape[1]{0.2};
    constexpr double thermalization = 0.4;
    constexpr double dt = 0.5;
    double trial[1]{-1.0};
    double heat[1]{-1.0};
    fusion_kinetic_ledger_v1 ledger;
    fill_ledger(&ledger, -1.0);

    const int status = fusion_c_energy_fp_trial(
        1, 0, dt, edges, old, nullptr, nullptr, birth, escape,
        thermalization, trial, heat, &ledger);
    check(status == PB11_STATUS_OK, "one-cell analytic step returns OK");

    const double cell_energy = 2.0e-16;
    const double expected_final =
        (old[0] + dt * birth[0]) /
        (1.0 + dt * (escape[0] + thermalization));
    check(close_relative(trial[0], expected_final, 1.0e-14),
          "one-cell step matches implicit birth/loss solution");
    check(old[0] == 4.0e20,
          "one-cell trial step leaves the accepted old state unchanged");
    check(heat[0] == -1.0,
          "zero-bath call does not write a nonexistent bath slot");

    const double expected_born = dt * birth[0];
    const double expected_escaped = dt * escape[0] * expected_final;
    const double expected_thermalized =
        dt * thermalization * expected_final;
    check(close_relative(ledger.initial_number_m3, old[0], 1.0e-14),
          "one-cell ledger records initial particles");
    check(close_relative(ledger.final_number_m3, expected_final, 1.0e-14),
          "one-cell ledger records returned particles");
    check(close_relative(ledger.born_number_m3, expected_born, 1.0e-14),
          "one-cell ledger records birth particles");
    check(close_relative(ledger.escaped_number_m3, expected_escaped, 1.0e-14),
          "one-cell ledger records escape particles once");
    check(close_relative(ledger.thermalized_number_m3, expected_thermalized,
                         1.0e-14),
          "one-cell ledger records thermalized particles once");

    check(close_relative(ledger.initial_energy_J_m3,
                         cell_energy * old[0], 1.0e-14),
          "one-cell ledger records initial energy");
    check(close_relative(ledger.final_energy_J_m3,
                         cell_energy * expected_final, 1.0e-14),
          "one-cell ledger records final energy");
    check(close_relative(ledger.born_energy_J_m3,
                         cell_energy * expected_born, 1.0e-14),
          "one-cell ledger records birth energy");
    check(close_relative(ledger.escaped_energy_J_m3,
                         cell_energy * expected_escaped, 1.0e-14),
          "one-cell ledger records escaped energy once");
    check(close_relative(ledger.thermalized_energy_J_m3,
                         cell_energy * expected_thermalized, 1.0e-14),
          "one-cell ledger records thermalized energy once");
    check(std::abs(ledger.particle_balance_error_m3) < 1.0e10,
          "one-cell particle balance error is small");
    check(std::abs(ledger.energy_balance_error_J_m3) < 1.0e-5,
          "one-cell energy balance error is small");
    check_ledger_balance(ledger, "one-cell");
}

std::array<double, 3> discrete_maxwellian(const std::array<double, 4>& edges,
                                          double temperature,
                                          double scale) {
    std::array<double, 3> values{};
    for (int i = 0; i < 3; ++i) {
        const double width = edges[i + 1] - edges[i];
        const double center = 0.5 * (edges[i + 1] + edges[i]);
        values[i] = scale * width * std::sqrt(center) *
                    std::exp(-center / temperature);
    }
    return values;
}

void test_common_bath_equilibrium_and_reflection() {
    const std::array<double, 4> edges{{1.0e-16, 2.0e-16, 4.0e-16,
                                       8.0e-16}};
    const double temperature = 4.0e-16;
    const std::array<double, 3> equilibrium =
        discrete_maxwellian(edges, temperature, 1.0e40);
    const double bath_temperature[1]{temperature};
    const double diffusion[2]{kDiffusion, kDiffusion};
    const double birth[3]{0.0, 0.0, 0.0};
    const double escape[3]{0.0, 0.0, 0.0};
    double trial_first[3]{-1.0, -1.0, -1.0};
    double trial_second[3]{-1.0, -1.0, -1.0};
    double heat_first[1]{-1.0};
    double heat_second[1]{-1.0};
    fusion_kinetic_ledger_v1 ledger_first;
    fusion_kinetic_ledger_v1 ledger_second;
    fill_ledger(&ledger_first, -1.0);
    fill_ledger(&ledger_second, -1.0);

    const int status_first = fusion_c_energy_fp_trial(
        3, 1, 0.1, edges.data(), equilibrium.data(), bath_temperature,
        diffusion, birth, escape, 0.0, trial_first, heat_first,
        &ledger_first);
    const int status_second = fusion_c_energy_fp_trial(
        3, 1, 0.1, edges.data(), equilibrium.data(), bath_temperature,
        diffusion, birth, escape, 0.0, trial_second, heat_second,
        &ledger_second);
    check(status_first == PB11_STATUS_OK && status_second == PB11_STATUS_OK,
          "common-bath Maxwellian trials return OK");
    for (int i = 0; i < 3; ++i) {
        check(close_relative(trial_first[i], equilibrium[i], 1.0e-11),
              "discrete Maxwellian is invariant under a common bath");
        check(trial_first[i] == trial_second[i],
              "repeated trials return identical cell populations");
    }
    check(std::abs(heat_first[0]) < 1.0e-8,
          "common-bath Maxwellian has zero bath heat to roundoff");
    check(heat_first[0] == heat_second[0],
          "repeated trials return identical bath heat");
    check(std::memcmp(&ledger_first, &ledger_second, sizeof(ledger_first)) ==
              0,
          "repeated trials return an identical ledger");
    check_ledger_balance(ledger_first, "common-bath equilibrium");
    for (int i = 0; i < 3; ++i) {
        check(equilibrium[i] == discrete_maxwellian(edges, temperature, 1.0e40)[i],
              "equilibrium input remains unchanged after trial");
    }

    const double arbitrary_old[3]{1.0e20, 2.0e20, 4.0e20};
    double arbitrary_trial[3]{-1.0, -1.0, -1.0};
    double arbitrary_heat[1]{-1.0};
    fusion_kinetic_ledger_v1 arbitrary_ledger;
    const int arbitrary_status = fusion_c_energy_fp_trial(
        3, 1, 0.1, edges.data(), arbitrary_old, bath_temperature, diffusion,
        birth, escape, 0.0, arbitrary_trial, arbitrary_heat,
        &arbitrary_ledger);
    check(arbitrary_status == PB11_STATUS_OK,
          "reflecting arbitrary distribution returns OK");
    check_finite_nonnegative(arbitrary_trial, 3,
                             "reflected arbitrary trial populations");
    const double old_total = arbitrary_old[0] + arbitrary_old[1] + arbitrary_old[2];
    const double trial_total = arbitrary_trial[0] + arbitrary_trial[1] +
                               arbitrary_trial[2];
    check(close_relative(trial_total, old_total, 1.0e-12),
          "reflecting energy grid conserves particles");
    check(close_relative(arbitrary_ledger.initial_number_m3, old_total,
                         1.0e-12),
          "reflection ledger records initial particle total");
    check(close_relative(arbitrary_ledger.final_number_m3, trial_total,
                         1.0e-12),
          "reflection ledger records final particle total");
}

void test_bath_exchange_and_energy_ledger() {
    const double edges[3]{1.0e-16, 2.0e-16, 4.0e-16};
    const double old[2]{0.0, 1.0e20};
    const double bath_temperature[1]{0.8e-16};
    const double diffusion[1]{kDiffusion};
    const double birth[2]{0.0, 0.0};
    const double escape[2]{0.0, 0.0};
    double trial[2]{-1.0, -1.0};
    double heat[1]{-1.0};
    fusion_kinetic_ledger_v1 ledger;

    const int status = fusion_c_energy_fp_trial(
        2, 1, 0.1, edges, old, bath_temperature, diffusion, birth, escape,
        0.0, trial, heat, &ledger);
    check(status == PB11_STATUS_OK,
          "hot fast group coupled to a cold bath returns OK");
    check_finite_nonnegative(trial, 2, "cold-bath trial populations");
    check(heat[0] > 0.0,
          "hot fast group deposits positive heat into a colder bath");
    check(ledger.final_energy_J_m3 < ledger.initial_energy_J_m3,
          "cold-bath cooling lowers fast-particle energy");
    check_ledger_balance(ledger, "cold-bath exchange", heat[0]);
    const double cold_energy_scale = std::abs(ledger.initial_energy_J_m3) +
                                     std::abs(ledger.final_energy_J_m3) +
                                     std::abs(heat[0]) + 1.0;
    check(std::abs(ledger.final_energy_J_m3 -
                   ledger.initial_energy_J_m3 + heat[0]) <=
              1.0e-10 * cold_energy_scale,
          "cold-bath heat closes the fast-particle energy change");

    const std::array<double, 4> multi_edges{{1.0e-16, 2.0e-16, 4.0e-16,
                                             8.0e-16}};
    const std::array<double, 3> middle_state =
        discrete_maxwellian(multi_edges, 4.0e-16, 1.0e40);
    const double bath_temperatures[2]{1.0e-16, 8.0e-16};
    const double multi_diffusion[4]{kDiffusion, kDiffusion, kDiffusion,
                                    kDiffusion};
    const double multi_birth[3]{0.0, 0.0, 0.0};
    const double multi_escape[3]{0.0, 0.0, 0.0};
    double multi_trial[3]{-1.0, -1.0, -1.0};
    double multi_heat[2]{-1.0, -1.0};
    fusion_kinetic_ledger_v1 multi_ledger;
    const int multi_status = fusion_c_energy_fp_trial(
        3, 2, 0.05, multi_edges.data(), middle_state.data(),
        bath_temperatures, multi_diffusion, multi_birth, multi_escape, 0.0,
        multi_trial, multi_heat, &multi_ledger);
    check(multi_status == PB11_STATUS_OK,
          "two-bath hot/cold exchange returns OK");
    check(multi_heat[0] > 0.0,
          "colder bath receives heat from the intermediate distribution");
    check(multi_heat[1] < 0.0,
          "hotter bath supplies heat to the intermediate distribution");
    check_finite_nonnegative(multi_trial, 3,
                             "two-bath trial populations");
    check_ledger_balance(multi_ledger, "two-bath exchange",
                         multi_heat[0] + multi_heat[1]);
    const double multi_energy_scale =
        std::abs(multi_ledger.initial_energy_J_m3) +
        std::abs(multi_ledger.final_energy_J_m3) + std::abs(multi_heat[0]) +
        std::abs(multi_heat[1]) + 1.0;
    check(std::abs(multi_ledger.final_energy_J_m3 -
                   multi_ledger.initial_energy_J_m3 + multi_heat[0] +
                   multi_heat[1]) <=
              1.0e-10 * multi_energy_scale,
          "two-bath heat sum closes the total energy ledger");
}

double integrate_one_cell(double dt, int steps) {
    const double edges[2]{1.0e-16, 3.0e-16};
    const double birth[1]{1.0e20};
    const double escape[1]{0.5};
    constexpr double thermalization = 0.2;
    double old[1]{3.0e20};
    double trial[1]{0.0};
    double heat[1]{0.0};
    fusion_kinetic_ledger_v1 ledger;
    for (int step = 0; step < steps; ++step) {
        const int status = fusion_c_energy_fp_trial(
            1, 0, dt, edges, old, nullptr, nullptr, birth, escape,
            thermalization, trial, heat, &ledger);
        check(status == PB11_STATUS_OK,
              "time-refinement implicit Euler step returns OK");
        old[0] = trial[0];
    }
    return old[0];
}

void test_time_refinement() {
    constexpr double initial = 3.0e20;
    constexpr double birth = 1.0e20;
    constexpr double lambda = 0.7;
    constexpr double final_time = 2.0;
    const double exact = birth / lambda +
                         (initial - birth / lambda) *
                             std::exp(-lambda * final_time);
    const double coarse = integrate_one_cell(0.2, 10);
    const double fine = integrate_one_cell(0.1, 20);
    const double coarse_error = std::abs(coarse - exact);
    const double fine_error = std::abs(fine - exact);
    check(fine_error < coarse_error,
          "halving the implicit Euler step reduces the error");
    check(coarse_error / fine_error > 1.7 && coarse_error / fine_error < 2.3,
          "implicit Euler convergence is first order");
}

void check_failure(const double* edges, const double* old,
                   const double* birth, const double* escape, double dt,
                   int expected_status, const std::string& label) {
    const double bath_temperature[1]{1.0e-16};
    const double diffusion[1]{kDiffusion};
    double trial[2]{-1.0, -1.0};
    double heat[1]{-1.0};
    fusion_kinetic_ledger_v1 ledger;
    fill_ledger(&ledger, -1.0);
    const int status = fusion_c_energy_fp_trial(
        2, 1, dt, edges, old, bath_temperature, diffusion, birth, escape,
        0.0, trial, heat, &ledger);
    check(status == expected_status, label + " returns expected status");
    check(trial[0] == 0.0 && trial[1] == 0.0 && heat[0] == 0.0 &&
              all_zero(ledger),
          label + " clears all outputs");
}

void test_zero_invalid_and_overflow() {
    const double edges_one[2]{1.0e-16, 3.0e-16};
    const double zero[1]{0.0};
    const double birth_one[1]{0.0};
    const double escape_one[1]{0.0};
    double trial_one[1]{-1.0};
    double heat_one[1]{-1.0};
    fusion_kinetic_ledger_v1 zero_ledger;
    fill_ledger(&zero_ledger, -1.0);
    int status = fusion_c_energy_fp_trial(
        1, 0, 0.1, edges_one, zero, nullptr, nullptr, birth_one, escape_one,
        0.0, trial_one, heat_one, &zero_ledger);
    check(status == PB11_STATUS_OK && trial_one[0] == 0.0 &&
              all_zero(zero_ledger),
          "zero one-cell data return zero state and ledger");

    const double edges_two[3]{1.0e-16, 2.0e-16, 4.0e-16};
    const double old_two[2]{1.0e20, 2.0e20};
    const double birth_two[2]{0.0, 0.0};
    const double escape_two[2]{0.0, 0.0};
    check_failure(edges_two, old_two, birth_two, escape_two, 0.0,
                  PB11_STATUS_OUT_OF_RANGE, "zero time step");

    const double nan_old[2]{std::numeric_limits<double>::quiet_NaN(), 2.0e20};
    check_failure(edges_two, nan_old, birth_two, escape_two, 0.1,
                  PB11_STATUS_INVALID_ARGUMENT, "NaN old population");

    const double bad_edges[3]{1.0e-16, 4.0e-16, 2.0e-16};
    check_failure(bad_edges, old_two, birth_two, escape_two, 0.1,
                  PB11_STATUS_OUT_OF_RANGE, "non-increasing edges");

    double trial_null[2]{-1.0, -1.0};
    double heat_null[1]{-1.0};
    fusion_kinetic_ledger_v1 ledger_null;
    fill_ledger(&ledger_null, -1.0);
    status = fusion_c_energy_fp_trial(
        2, 1, 0.1, edges_two, old_two, nullptr, nullptr, birth_two,
        escape_two, 0.0, trial_null, heat_null, &ledger_null);
    check(status == PB11_STATUS_INVALID_ARGUMENT && trial_null[0] == 0.0 &&
              trial_null[1] == 0.0 && heat_null[0] == 0.0 &&
              all_zero(ledger_null),
          "missing required bath arrays clear valid-dimension outputs");

    check(fusion_c_energy_fp_trial(
              2, 1, 0.1, edges_two, old_two, nullptr, nullptr, birth_two,
              escape_two, 0.0, nullptr, heat_null, &ledger_null) ==
              PB11_STATUS_NULL_OUTPUT,
          "null trial output is reported");
    check(fusion_c_energy_fp_trial(
              2, 1, 0.1, edges_two, old_two, nullptr, nullptr, birth_two,
              escape_two, 0.0, trial_null, nullptr, &ledger_null) ==
              PB11_STATUS_NULL_OUTPUT,
          "null bath-heat output is reported");
    check(fusion_c_energy_fp_trial(
              2, 1, 0.1, edges_two, old_two, nullptr, nullptr, birth_two,
              escape_two, 0.0, trial_null, heat_null, nullptr) ==
              PB11_STATUS_NULL_OUTPUT,
          "null ledger output is reported");

    const double huge_old[1]{std::numeric_limits<double>::max()};
    const double huge_birth[1]{std::numeric_limits<double>::max()};
    const double huge_escape[1]{0.0};
    double huge_trial[1]{-1.0};
    double huge_heat[1]{-1.0};
    fusion_kinetic_ledger_v1 huge_ledger;
    fill_ledger(&huge_ledger, -1.0);
    status = fusion_c_energy_fp_trial(
        1, 0, 1.0, edges_one, huge_old, nullptr, nullptr, huge_birth,
        huge_escape, 0.0, huge_trial, huge_heat, &huge_ledger);
    check(status == PB11_STATUS_NUMERICAL_FAILURE && huge_trial[0] == 0.0 &&
              all_zero(huge_ledger),
          "unrepresentable trial output is rejected and cleared");
}

}  // namespace

int main() {
    test_one_cell_analytic_step();
    test_common_bath_equilibrium_and_reflection();
    test_bath_exchange_and_energy_ledger();
    test_time_refinement();
    test_zero_invalid_and_overflow();

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All fusion-kinetics tests passed\n";
    return 0;
}
