#include "fusion_two_component.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr double kJoulesPerKeV = 1.602176634e-16;
constexpr double kAlphaMassKg = 6.644657230e-27;
constexpr double kProtonMassKg = 1.67262192369e-27;

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
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

void fill_ledger(fusion_two_component_ledger_v1& ledger, double value) {
    ledger.total.initial_number_m3 = value;
    ledger.total.final_number_m3 = value;
    ledger.total.initial_energy_J_m3 = value;
    ledger.total.final_energy_J_m3 = value;
    ledger.total.born_number_m3 = value;
    ledger.total.born_energy_J_m3 = value;
    ledger.total.escaped_number_m3 = value;
    ledger.total.escaped_energy_J_m3 = value;
    ledger.total.thermalized_number_m3 = value;
    ledger.total.thermalized_energy_J_m3 = value;
    ledger.total.particle_balance_error_m3 = value;
    ledger.total.energy_balance_error_J_m3 = value;
    ledger.transferred_number_m3 = value;
    ledger.transferred_energy_J_m3 = value;
}

bool zero_ledger(const fusion_two_component_ledger_v1& ledger) {
    return ledger.total.initial_number_m3 == 0.0 &&
           ledger.total.final_number_m3 == 0.0 &&
           ledger.total.initial_energy_J_m3 == 0.0 &&
           ledger.total.final_energy_J_m3 == 0.0 &&
           ledger.total.born_number_m3 == 0.0 &&
           ledger.total.born_energy_J_m3 == 0.0 &&
           ledger.total.escaped_number_m3 == 0.0 &&
           ledger.total.escaped_energy_J_m3 == 0.0 &&
           ledger.total.thermalized_number_m3 == 0.0 &&
           ledger.total.thermalized_energy_J_m3 == 0.0 &&
           ledger.total.particle_balance_error_m3 == 0.0 &&
           ledger.total.energy_balance_error_J_m3 == 0.0 &&
           ledger.transferred_number_m3 == 0.0 &&
           ledger.transferred_energy_J_m3 == 0.0;
}

void check_ledger_close(const fusion_two_component_ledger_v1& actual,
                        const fusion_kinetic_ledger_v1& expected,
                        const std::string& label) {
    check(close_relative(actual.total.initial_number_m3,
                         expected.initial_number_m3, 2.0e-12),
          label + " initial particle ledger");
    check(close_relative(actual.total.final_number_m3,
                         expected.final_number_m3, 2.0e-12),
          label + " final particle ledger");
    check(close_relative(actual.total.initial_energy_J_m3,
                         expected.initial_energy_J_m3, 2.0e-12),
          label + " initial energy ledger");
    check(close_relative(actual.total.final_energy_J_m3,
                         expected.final_energy_J_m3, 2.0e-12),
          label + " final energy ledger");
    check(close_relative(actual.total.born_number_m3,
                         expected.born_number_m3, 2.0e-12),
          label + " external birth number");
    check(close_relative(actual.total.born_energy_J_m3,
                         expected.born_energy_J_m3, 2.0e-12),
          label + " external birth energy");
    check(close_relative(actual.total.escaped_number_m3,
                         expected.escaped_number_m3, 2.0e-12),
          label + " physical escape number");
    check(close_relative(actual.total.escaped_energy_J_m3,
                         expected.escaped_energy_J_m3, 2.0e-12),
          label + " physical escape energy");
    check(actual.total.thermalized_number_m3 == 0.0 &&
              actual.total.thermalized_energy_J_m3 == 0.0,
          label + " has no thermalized ledger");
}

void test_one_cell_analytic() {
    const double edges[2]{1.0e-16, 3.0e-16};
    const double old_s[1]{4.0e20};
    const double old_t[1]{1.0e20};
    const double birth_s[1]{3.0e20};
    const double birth_t[1]{2.0e20};
    const double escape[1]{0.2};
    const double transfer[1]{0.4};
    constexpr double dt = 0.5;
    double trial_s[1]{-1.0};
    double trial_t[1]{-1.0};
    fusion_two_component_ledger_v1 ledger;
    fill_ledger(ledger, -1.0);

    const int status = fusion_c_two_component_trial(
        1, 0, dt, edges, old_s, old_t, nullptr, nullptr, birth_s, birth_t,
        escape, transfer, trial_s, trial_t, nullptr, &ledger);
    check(status == PB11_STATUS_OK, "one-cell two-component step returns OK");

    const double expected_s = (old_s[0] + dt * birth_s[0]) /
                              (1.0 + dt * (escape[0] + transfer[0]));
    const double expected_t =
        (old_t[0] + dt * (birth_t[0] + transfer[0] * expected_s)) /
        (1.0 + dt * escape[0]);
    check(close_relative(trial_s[0], expected_s, 2.0e-14),
          "one-cell S matches analytic implicit solution");
    check(close_relative(trial_t[0], expected_t, 2.0e-14),
          "one-cell T matches analytic implicit solution");

    const double cell_energy = 2.0e-16;
    const double expected_born = dt * (birth_s[0] + birth_t[0]);
    const double expected_escape =
        dt * escape[0] * (expected_s + expected_t);
    const double expected_transfer = dt * transfer[0] * expected_s;
    check(close_relative(ledger.total.born_number_m3, expected_born, 2.0e-14),
          "one-cell total birth excludes internal transfer");
    check(close_relative(ledger.total.escaped_number_m3, expected_escape,
                         2.0e-14),
          "one-cell total escape contains physical escape only");
    check(close_relative(ledger.transferred_number_m3, expected_transfer,
                         2.0e-14),
          "one-cell internal transfer number is separate");
    check(close_relative(ledger.total.born_energy_J_m3,
                         cell_energy * expected_born, 2.0e-14),
          "one-cell total birth energy is external only");
    check(close_relative(ledger.total.escaped_energy_J_m3,
                         cell_energy * expected_escape, 2.0e-14),
          "one-cell total escape energy is physical only");
    check(close_relative(ledger.transferred_energy_J_m3,
                         cell_energy * expected_transfer, 2.0e-14),
          "one-cell carried transfer energy is separate");
    check(ledger.total.thermalized_number_m3 == 0.0 &&
              ledger.total.thermalized_energy_J_m3 == 0.0,
          "one-cell total thermalized fields are exactly zero");
    check(std::abs(ledger.total.particle_balance_error_m3) < 1.0e10,
          "one-cell particle ledger closes");
    check(std::abs(ledger.total.energy_balance_error_J_m3) < 1.0e-5,
          "one-cell energy ledger closes");
}

bool coulomb_diffusion(const std::array<double, 5>& edges,
                       const double bath_kT[2], const double density[2],
                       const double coulomb_log[2], double diffusion[6]) {
    for (int b = 0; b < 2; ++b) {
        fusion_maxwellian_bath_v1 bath{
            density[b], kProtonMassKg, 1.0, bath_kT[b], coulomb_log[b]};
        for (int f = 0; f < 3; ++f) {
            const double energy =
                0.5 * (edges[static_cast<std::size_t>(f + 1)] +
                       edges[static_cast<std::size_t>(f + 2)]);
            fusion_coulomb_energy_v1 coefficient{};
            const int status = fusion_c_coulomb_energy(
                energy, kAlphaMassKg, 2.0, &bath, &coefficient);
            check(status == PB11_STATUS_OK,
                  "actual Coulomb coefficient returns OK");
            if (status != PB11_STATUS_OK ||
                !std::isfinite(coefficient.diffusion_J2_s) ||
                coefficient.diffusion_J2_s < 0.0) {
                return false;
            }
            diffusion[static_cast<std::size_t>(b * 3 + f)] =
                coefficient.diffusion_J2_s;
        }
    }
    return true;
}

void test_multicell_coulomb_full_parity() {
    constexpr int cells = 4;
    constexpr int baths = 2;
    const std::array<double, 5> edges{{1.0e-16, 1.5e-16, 3.0e-16,
                                       6.0e-16, 1.2e-15}};
    std::array<double, cells> old_s{{1.0e19, 2.0e19, 3.0e19, 4.0e19}};
    std::array<double, cells> old_t{{2.0e19, 1.0e19, 0.5e19, 3.0e19}};

    for (int step = 0; step < 4; ++step) {
        const double dt = 0.02 + 0.005 * static_cast<double>(step);
        const double bath_kT[baths]{
            (0.8 + 0.1 * step) * kJoulesPerKeV,
            (4.0 - 0.2 * step) * kJoulesPerKeV};
        const double density[baths]{
            2.0e20 * (1.0 + 0.1 * step),
            5.0e19 * (1.0 + 0.05 * step)};
        const double coulomb_log[baths]{12.0 + step, 14.0 - 0.25 * step};
        double diffusion[baths * (cells - 1)]{};
        check(coulomb_diffusion(edges, bath_kT, density, coulomb_log,
                                diffusion),
              "time-varying actual Coulomb diffusion is available");

        std::array<double, cells> birth_s{};
        std::array<double, cells> birth_t{};
        std::array<double, cells> escape{};
        std::array<double, cells> transfer{};
        for (int i = 0; i < cells; ++i) {
            birth_s[static_cast<std::size_t>(i)] =
                (0.5 + 0.1 * i + 0.05 * step) * 1.0e18;
            birth_t[static_cast<std::size_t>(i)] =
                (0.2 + 0.07 * i + 0.03 * step) * 1.0e18;
            escape[static_cast<std::size_t>(i)] =
                0.01 + 0.002 * i + 0.001 * step;
            transfer[static_cast<std::size_t>(i)] =
                0.03 + 0.01 * i + 0.002 * step;
        }

        std::array<double, cells> two_s{};
        std::array<double, cells> two_t{};
        std::array<double, baths> two_heat{};
        fusion_two_component_ledger_v1 two_ledger;
        fill_ledger(two_ledger, -1.0);
        const int two_status = fusion_c_two_component_trial(
            cells, baths, dt, edges.data(), old_s.data(), old_t.data(),
            bath_kT, diffusion, birth_s.data(), birth_t.data(), escape.data(),
            transfer.data(), two_s.data(), two_t.data(), two_heat.data(),
            &two_ledger);
        check(two_status == PB11_STATUS_OK,
              "multi-cell two-component trial returns OK");

        std::array<double, cells> full_old{};
        std::array<double, cells> full_birth{};
        for (int i = 0; i < cells; ++i) {
            full_old[static_cast<std::size_t>(i)] =
                old_s[static_cast<std::size_t>(i)] +
                old_t[static_cast<std::size_t>(i)];
            full_birth[static_cast<std::size_t>(i)] =
                birth_s[static_cast<std::size_t>(i)] +
                birth_t[static_cast<std::size_t>(i)];
        }
        std::array<double, cells> full_trial{};
        std::array<double, baths> full_heat{};
        fusion_kinetic_ledger_v1 full_ledger{};
        const int full_status = fusion_c_energy_fp_trial(
            cells, baths, dt, edges.data(), full_old.data(), bath_kT, diffusion,
            full_birth.data(), escape.data(), 0.0, full_trial.data(),
            full_heat.data(), &full_ledger);
        check(full_status == PB11_STATUS_OK,
              "single-component FULL comparison trial returns OK");

        double l1 = 0.0;
        double full_scale = 0.0;
        for (int i = 0; i < cells; ++i) {
            const double sum = two_s[static_cast<std::size_t>(i)] +
                               two_t[static_cast<std::size_t>(i)];
            const double difference = sum - full_trial[static_cast<std::size_t>(i)];
            l1 += std::abs(difference);
            full_scale += std::abs(full_trial[static_cast<std::size_t>(i)]);
            check(close_relative(sum, full_trial[static_cast<std::size_t>(i)],
                                 2.0e-11),
                  "S+T equals FULL cell value at step " +
                      std::to_string(step));
        }
        /* Complete L1 average: no legacy 0.5 factor. */
        l1 /= static_cast<double>(cells);
        full_scale /= static_cast<double>(cells);
        check(l1 <= 2.0e-11 * std::max(1.0, full_scale),
              "complete L1(S+T,FULL) parity at step " +
                  std::to_string(step));
        for (int b = 0; b < baths; ++b) {
            check(close_relative(two_heat[static_cast<std::size_t>(b)],
                                 full_heat[static_cast<std::size_t>(b)],
                                 2.0e-11, 1.0e-30),
                  "S+T heat equals FULL heat at step " +
                      std::to_string(step));
        }
        check_ledger_close(two_ledger, full_ledger,
                           "two-component versus FULL step " +
                               std::to_string(step));

        check(std::isfinite(two_ledger.transferred_number_m3) &&
                  two_ledger.transferred_number_m3 > 0.0,
              "multi-cell transfer number is positive and finite");
        check(std::isfinite(two_ledger.transferred_energy_J_m3) &&
                  two_ledger.transferred_energy_J_m3 > 0.0,
              "multi-cell transfer energy is positive and finite");

        old_s = two_s;
        old_t = two_t;
    }
}

void test_rollback_repeat() {
    constexpr int cells = 3;
    constexpr int baths = 1;
    const double edges[cells + 1]{1.0e-16, 2.0e-16, 4.0e-16, 8.0e-16};
    const double old_s[cells]{1.0e20, 2.0e20, 3.0e20};
    const double old_t[cells]{3.0e20, 2.0e20, 1.0e20};
    const double bath_kT[baths]{4.0e-16};
    const double diffusion[baths * (cells - 1)]{1.0e-30, 2.0e-30};
    const double birth_s[cells]{1.0e18, 2.0e18, 3.0e18};
    const double birth_t[cells]{4.0e18, 5.0e18, 6.0e18};
    const double escape[cells]{0.02, 0.03, 0.04};
    const double transfer[cells]{0.1, 0.2, 0.3};
    double first_s[cells]{}, first_t[cells]{}, first_heat[baths]{};
    double second_s[cells]{}, second_t[cells]{}, second_heat[baths]{};
    fusion_two_component_ledger_v1 first_ledger{}, second_ledger{};

    const int first_status = fusion_c_two_component_trial(
        cells, baths, 0.1, edges, old_s, old_t, bath_kT, diffusion, birth_s,
        birth_t, escape, transfer, first_s, first_t, first_heat, &first_ledger);
    const int second_status = fusion_c_two_component_trial(
        cells, baths, 0.1, edges, old_s, old_t, bath_kT, diffusion, birth_s,
        birth_t, escape, transfer, second_s, second_t, second_heat,
        &second_ledger);
    check(first_status == PB11_STATUS_OK && second_status == PB11_STATUS_OK,
          "repeated trial calls return OK");
    check(std::memcmp(first_s, second_s, sizeof(first_s)) == 0 &&
              std::memcmp(first_t, second_t, sizeof(first_t)) == 0 &&
              std::memcmp(first_heat, second_heat, sizeof(first_heat)) == 0 &&
              std::memcmp(&first_ledger, &second_ledger,
                          sizeof(first_ledger)) == 0,
          "repeated trial calls are bitwise identical for rollback");
    check(old_s[0] == 1.0e20 && old_s[1] == 2.0e20 && old_s[2] == 3.0e20 &&
              old_t[0] == 3.0e20 && old_t[1] == 2.0e20 && old_t[2] == 1.0e20,
          "rollback trials leave accepted old states unchanged");
}

void test_invalid_and_clear() {
    const double edges[2]{1.0e-16, 3.0e-16};
    const double old_s[1]{1.0e20};
    const double old_t[1]{2.0e20};
    const double birth_s[1]{1.0e18};
    const double birth_t[1]{2.0e18};
    const double escape[1]{0.1};
    const double transfer[1]{0.2};
    double trial_s[1]{-1.0};
    double trial_t[1]{-1.0};
    fusion_two_component_ledger_v1 ledger;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    fill_ledger(ledger, -1.0);
    int status = fusion_c_two_component_trial(
        1, 0, 0.1, edges, old_s, old_t, nullptr, nullptr, birth_s, birth_t,
        escape, &nan, trial_s, trial_t, nullptr, &ledger);
    check(status == PB11_STATUS_INVALID_ARGUMENT && trial_s[0] == 0.0 &&
              trial_t[0] == 0.0 && zero_ledger(ledger),
          "nonfinite transfer rejects and clears outputs");

    fill_ledger(ledger, -1.0);
    trial_s[0] = trial_t[0] = -1.0;
    const double negative_old[1]{-1.0};
    status = fusion_c_two_component_trial(
        1, 0, 0.1, edges, negative_old, old_t, nullptr, nullptr, birth_s,
        birth_t, escape, transfer, trial_s, trial_t, nullptr, &ledger);
    check(status == PB11_STATUS_INVALID_ARGUMENT && trial_s[0] == 0.0 &&
              trial_t[0] == 0.0 && zero_ledger(ledger),
          "negative input rejects without clipping and clears outputs");

    fill_ledger(ledger, -1.0);
    trial_s[0] = trial_t[0] = -1.0;
    status = fusion_c_two_component_trial(
        1, 0, 0.1, edges, old_s, old_t, nullptr, nullptr, birth_s, birth_t,
        escape, transfer, nullptr, trial_t, nullptr, &ledger);
    check(status == PB11_STATUS_NULL_OUTPUT && trial_t[0] == 0.0 &&
              zero_ledger(ledger),
          "null S output is reported and remaining outputs clear");

    fill_ledger(ledger, -1.0);
    trial_s[0] = trial_t[0] = -1.0;
    status = fusion_c_two_component_trial(
        1, 1, 0.1, edges, old_s, old_t, &escape[0], nullptr, birth_s,
        birth_t, escape, transfer, trial_s, trial_t, nullptr, &ledger);
    check(status == PB11_STATUS_NULL_OUTPUT && trial_s[0] == 0.0 &&
              trial_t[0] == 0.0 && zero_ledger(ledger),
          "missing bath heat output is reported and arrays clear");

    const double one_bath_kT[1]{4.0e-16};
    double one_bath_heat[1]{-1.0};
    trial_s[0] = trial_t[0] = -1.0;
    fill_ledger(ledger, -1.0);
    status = fusion_c_two_component_trial(
        1, 1, 0.1, edges, old_s, old_t, one_bath_kT, nullptr, birth_s,
        birth_t, escape, transfer, trial_s, trial_t, one_bath_heat, &ledger);
    check(status == PB11_STATUS_OK && std::isfinite(trial_s[0]) &&
              std::isfinite(trial_t[0]) && one_bath_heat[0] == 0.0,
          "one-cell bath accepts null diffusion and has zero collision heat");

    trial_s[0] = trial_t[0] = -1.0;
    one_bath_heat[0] = -1.0;
    fill_ledger(ledger, -1.0);
    status = fusion_c_two_component_trial(
        1, 1, 0.1, edges, old_s, old_t, nullptr, nullptr, birth_s, birth_t,
        escape, transfer, trial_s, trial_t, one_bath_heat, &ledger);
    check(status == PB11_STATUS_INVALID_ARGUMENT && trial_s[0] == 0.0 &&
              trial_t[0] == 0.0 && one_bath_heat[0] == 0.0 &&
              zero_ledger(ledger),
          "missing bath temperature rejects and clears outputs");

    const double huge = std::numeric_limits<double>::max();
    const double huge_old[1]{0.0};
    const double huge_birth[1]{huge};
    const double huge_transfer[1]{huge};
    const double zero_escape[1]{huge};
    trial_s[0] = trial_t[0] = -1.0;
    fill_ledger(ledger, -1.0);
    status = fusion_c_two_component_trial(
        1, 0, 1.0, edges, huge_old, huge_old, nullptr, nullptr, huge_birth,
        huge_old, zero_escape, huge_transfer, trial_s, trial_t, nullptr,
        &ledger);
    check(status == PB11_STATUS_NUMERICAL_FAILURE && trial_s[0] == 0.0 &&
              trial_t[0] == 0.0 && zero_ledger(ledger),
          "derived escape overflow rejects and clears outputs");

    trial_s[0] = trial_t[0] = -1.0;
    fill_ledger(ledger, -1.0);
    status = fusion_c_two_component_trial(
        0, 0, 0.1, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, trial_s, trial_t, nullptr, &ledger);
    check(status == PB11_STATUS_INVALID_ARGUMENT && zero_ledger(ledger),
          "zero cell count is rejected");
}

void test_transfer_tail_underflow() {
    const double edges[3]={0.,1.,2.}, old_s[2]={4.,1e-200}, old_t[2]={1.,0.};
    const double zero[2]={0.,0.}, transfer[2]={.4,1e-200};
    double sn[2],tn[2]; fusion_two_component_ledger_v1 ledger{};
    const int status=fusion_c_two_component_trial(2,0,.5,edges,old_s,old_t,
        nullptr,nullptr,zero,zero,zero,transfer,sn,tn,nullptr,&ledger);
    check(status==PB11_STATUS_OK,"negligible per-cell transfer underflow accepted");
    check(close_relative(sn[0],4./1.2,1e-14)&&
          close_relative(tn[0],1.+.2*sn[0],1e-14),"resolved transfer preserved");
    check(sn[1]==old_s[1] && tn[1]==0.,"positive tail rounded without state clipping");
    check(std::abs(ledger.total.particle_balance_error_m3)<1e-14 &&
          std::abs(ledger.total.energy_balance_error_J_m3)<1e-14,
          "combined ledger bounds transfer rounding loss");
}

}  // namespace

int main() {
    test_one_cell_analytic();
    test_multicell_coulomb_full_parity();
    test_rollback_repeat();
    test_invalid_and_clear();
    test_transfer_tail_underflow();

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All two-component fusion tests passed\n";
    return 0;
}
