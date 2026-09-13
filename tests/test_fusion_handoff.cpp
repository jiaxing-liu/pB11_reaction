#include "fusion_handoff.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kJoulesPerKeV = 1.602176634e-16;

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

/* Gamma(3/2) regularized CDF, independently evaluated through erf.  The
 * tests use x away from zero so this comparison does not test cancellation in
 * the independent reference itself. */
double gamma32_cdf(double energy_J, double kT_J) {
    if (energy_J <= 0.0) {
        return 0.0;
    }
    const double x = energy_J / kT_J;
    const double y = std::sqrt(x);
    return std::erf(y) - (2.0 / std::sqrt(kPi)) * y * std::exp(-x);
}

void fill_grid(fusion_maxwellian_grid_v1& out, double value) {
    out.below_probability = value;
    out.above_probability = value;
    out.represented_mean_energy_J = value;
    out.probability_balance_error = value;
}

bool zero_grid(const fusion_maxwellian_grid_v1& out) {
    return out.below_probability == 0.0 && out.above_probability == 0.0 &&
           out.represented_mean_energy_J == 0.0 &&
           out.probability_balance_error == 0.0;
}

void fill_handoff(fusion_handoff_ledger_v1& out, double value) {
    out.initial_number_m3 = value;
    out.initial_energy_J_m3 = value;
    out.remaining_number_m3 = value;
    out.remaining_energy_J_m3 = value;
    out.fluid_number_m3 = value;
    out.fluid_energy_J_m3 = value;
    out.bath_energy_correction_J_m3 = value;
    out.distribution_L1 = value;
    out.relative_mean_energy_error = value;
    out.outside_grid_probability = value;
    out.represented_Maxwellian_mean_energy_J = value;
    out.particle_balance_error_m3 = value;
    out.energy_balance_error_J_m3 = value;
}

bool zero_handoff(const fusion_handoff_ledger_v1& out) {
    return out.initial_number_m3 == 0.0 &&
           out.initial_energy_J_m3 == 0.0 &&
           out.remaining_number_m3 == 0.0 &&
           out.remaining_energy_J_m3 == 0.0 && out.fluid_number_m3 == 0.0 &&
           out.fluid_energy_J_m3 == 0.0 &&
           out.bath_energy_correction_J_m3 == 0.0 &&
           out.distribution_L1 == 0.0 &&
           out.relative_mean_energy_error == 0.0 &&
           out.outside_grid_probability == 0.0 &&
           out.represented_Maxwellian_mean_energy_J == 0.0 &&
           out.particle_balance_error_m3 == 0.0 &&
           out.energy_balance_error_J_m3 == 0.0;
}

bool same_handoff(const fusion_handoff_ledger_v1& lhs,
                  const fusion_handoff_ledger_v1& rhs) {
    return lhs.initial_number_m3 == rhs.initial_number_m3 &&
           lhs.initial_energy_J_m3 == rhs.initial_energy_J_m3 &&
           lhs.remaining_number_m3 == rhs.remaining_number_m3 &&
           lhs.remaining_energy_J_m3 == rhs.remaining_energy_J_m3 &&
           lhs.fluid_number_m3 == rhs.fluid_number_m3 &&
           lhs.fluid_energy_J_m3 == rhs.fluid_energy_J_m3 &&
           lhs.bath_energy_correction_J_m3 == rhs.bath_energy_correction_J_m3 &&
           lhs.distribution_L1 == rhs.distribution_L1 &&
           lhs.relative_mean_energy_error == rhs.relative_mean_energy_error &&
           lhs.outside_grid_probability == rhs.outside_grid_probability &&
           lhs.represented_Maxwellian_mean_energy_J ==
               rhs.represented_Maxwellian_mean_energy_J &&
           lhs.particle_balance_error_m3 == rhs.particle_balance_error_m3 &&
           lhs.energy_balance_error_J_m3 == rhs.energy_balance_error_J_m3;
}

void check_handoff_finite(const fusion_handoff_ledger_v1& out,
                          const std::string& label) {
    check(std::isfinite(out.initial_number_m3), label + " initial N finite");
    check(std::isfinite(out.initial_energy_J_m3), label + " initial U finite");
    check(std::isfinite(out.remaining_number_m3), label + " remaining N finite");
    check(std::isfinite(out.remaining_energy_J_m3), label + " remaining U finite");
    check(std::isfinite(out.fluid_number_m3), label + " fluid N finite");
    check(std::isfinite(out.fluid_energy_J_m3), label + " fluid U finite");
    check(std::isfinite(out.bath_energy_correction_J_m3),
          label + " bath correction finite");
    check(std::isfinite(out.distribution_L1), label + " L1 finite");
    check(std::isfinite(out.relative_mean_energy_error),
          label + " mean error finite");
    check(std::isfinite(out.outside_grid_probability),
          label + " outside probability finite");
    check(std::isfinite(out.represented_Maxwellian_mean_energy_J),
          label + " represented mean finite");
    check(std::isfinite(out.particle_balance_error_m3),
          label + " particle closure finite");
    check(std::isfinite(out.energy_balance_error_J_m3),
          label + " energy closure finite");
}

void test_gamma32_grid_and_truncation() {
    const double kT = 2.0 * kJoulesPerKeV;
    const std::array<double, 5> edges{{0.35 * kT, 0.80 * kT, 1.40 * kT,
                                       2.20 * kT, 3.70 * kT}};
    std::array<double, 4> probability{{-1.0, -1.0, -1.0, -1.0}};
    fusion_maxwellian_grid_v1 out;
    fill_grid(out, -1.0);

    const int status = fusion_c_maxwellian_energy_grid(
        static_cast<int>(probability.size()), kT, edges.data(),
        probability.data(), &out);
    check(status == PB11_STATUS_OK, "Gamma(3/2) grid returns OK");

    double actual_sum = 0.0;
    double expected_mean = 0.0;
    for (std::size_t i = 0; i < probability.size(); ++i) {
        const double expected = gamma32_cdf(edges[i + 1], kT) -
                                gamma32_cdf(edges[i], kT);
        actual_sum += probability[i];
        expected_mean += expected *
                         (edges[i] + edges[i + 1]) / 2.0;
        check(close_relative(probability[i], expected, 2.0e-12),
              "Gamma(3/2) bin probability agrees with independent erf CDF");
    }
    const double expected_below = gamma32_cdf(edges.front(), kT);
    const double expected_above = 1.0 - gamma32_cdf(edges.back(), kT);
    check(close_relative(out.below_probability, expected_below, 2.0e-12),
          "Gamma(3/2) below tail agrees with independent CDF");
    check(close_relative(out.above_probability, expected_above, 2.0e-12),
          "Gamma(3/2) above tail agrees with independent CDF");
    check(close_relative(actual_sum + out.below_probability +
                             out.above_probability,
                         1.0, 2.0e-12),
          "bin probabilities plus both tails normalize to one");
    check(close_relative(out.represented_mean_energy_J, expected_mean,
                         2.0e-12),
          "represented mean uses arithmetic cell centers");
    check(close_relative(out.probability_balance_error,
                         actual_sum + out.below_probability +
                             out.above_probability - 1.0,
                         2.0e-12, 1.0e-18),
          "grid probability balance error reports the full normalization");

    const std::array<double, 2> truncated_edges{{0.40 * kT, 2.00 * kT}};
    double truncated_probability[1]{-1.0};
    fusion_maxwellian_grid_v1 truncated_out;
    fill_grid(truncated_out, -1.0);
    const int truncated_status = fusion_c_maxwellian_energy_grid(
        1, kT, truncated_edges.data(), truncated_probability, &truncated_out);
    check(truncated_status == PB11_STATUS_OK,
          "truncated Maxwellian grid returns OK");
    const double truncated_sum = truncated_probability[0];
    check(truncated_sum < 1.0 && truncated_sum > 0.0,
          "truncated grid retains less than one full-distribution particle");
    check(close_relative(truncated_sum + truncated_out.below_probability +
                             truncated_out.above_probability,
                         1.0, 2.0e-12),
          "truncated grid tails complete normalization without renormalizing q");
    check(close_relative(truncated_sum,
                         gamma32_cdf(truncated_edges[1], kT) -
                             gamma32_cdf(truncated_edges[0], kT),
                         2.0e-12),
          "truncated q is the original probability mass, not a renormalized mass");
}

struct FineGrid {
    static constexpr int cells = 1000;
    double kT = 2.0 * kJoulesPerKeV;
    std::vector<double> edges;
    std::vector<double> probability;

    FineGrid() : edges(cells + 1), probability(cells) {
        for (int i = 0; i <= cells; ++i) {
            edges[static_cast<std::size_t>(i)] =
                kT * (20.0 * static_cast<double>(i) / cells);
        }
    }
};

bool make_fine_grid(FineGrid& grid) {
    fusion_maxwellian_grid_v1 diagnostics{};
    const int status = fusion_c_maxwellian_energy_grid(
        FineGrid::cells, grid.kT, grid.edges.data(), grid.probability.data(),
        &diagnostics);
    check(status == PB11_STATUS_OK, "fine Maxwellian grid returns OK");
    check(std::isfinite(diagnostics.represented_mean_energy_J),
          "fine Maxwellian represented mean is finite");
    check(diagnostics.above_probability < 1.0e-6,
          "fine grid has a small omitted high-energy tail");
    return status == PB11_STATUS_OK;
}

void calculate_old_moments(const FineGrid& grid, const std::vector<double>& old,
                           long double& number, long double& energy) {
    number = 0.0L;
    energy = 0.0L;
    for (int i = 0; i < FineGrid::cells; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const long double center =
            (static_cast<long double>(grid.edges[index]) +
             static_cast<long double>(grid.edges[index + 1])) /
            2.0L;
        number += static_cast<long double>(old[index]);
        energy += center * static_cast<long double>(old[index]);
    }
}

void check_rejected_handoff(const fusion_handoff_ledger_v1& out,
                            const std::vector<double>& old,
                            const std::vector<double>& trial, int projected,
                            const std::string& label) {
    check(projected == 0, label + " is not projected");
    check(trial == old, label + " returns the old candidate exactly");
    check(out.remaining_number_m3 == out.initial_number_m3,
          label + " retains candidate particle inventory");
    check(out.remaining_energy_J_m3 == out.initial_energy_J_m3,
          label + " retains candidate energy");
    check(out.fluid_number_m3 == 0.0 && out.fluid_energy_J_m3 == 0.0 &&
              out.bath_energy_correction_J_m3 == 0.0,
          label + " has no fluid or bath source on rejection");
    check(out.particle_balance_error_m3 == 0.0 &&
              out.energy_balance_error_J_m3 == 0.0,
          label + " rejection ledger has zero closure sources");
}

void test_fine_projection_and_signed_correction(const FineGrid& grid) {
    constexpr double full_number = 1.0e20;
    std::vector<double> old(grid.probability.size());
    for (std::size_t i = 0; i < old.size(); ++i) {
        old[i] = full_number * grid.probability[i];
    }
    const std::vector<double> old_copy = old;
    std::vector<double> trial(old.size(), -1.0);
    int projected = -1;
    fusion_handoff_ledger_v1 out;
    fill_handoff(out, -1.0);

    const double max_l1 = 1.0e-3;
    const double max_mean_error = 1.0e-4;
    const int status = fusion_c_maxwellian_handoff_trial(
        FineGrid::cells, grid.kT, max_l1, max_mean_error, grid.edges.data(),
        old.data(), trial.data(), &projected, &out);
    check(status == PB11_STATUS_OK, "fine Maxwellian handoff returns OK");
    check(projected == 1, "fine Maxwellian candidate is projected");
    check(std::all_of(trial.begin(), trial.end(),
                      [](double value) { return value == 0.0; }),
          "projected trial distribution is zero");
    check(old == old_copy, "projection trial leaves old state unchanged");
    check(out.remaining_number_m3 == 0.0 &&
              out.remaining_energy_J_m3 == 0.0,
          "projected ledger has no remaining kinetic candidate");

    long double old_number = 0.0L;
    long double old_energy = 0.0L;
    calculate_old_moments(grid, old, old_number, old_energy);
    const double expected_initial_number = static_cast<double>(old_number);
    const double expected_initial_energy = static_cast<double>(old_energy);
    const double expected_fluid_energy = static_cast<double>(
        1.5L * static_cast<long double>(grid.kT) *
        static_cast<long double>(out.fluid_number_m3));
    const double expected_correction =
        static_cast<double>(old_energy -
                            static_cast<long double>(expected_fluid_energy));

    check(close_relative(out.initial_number_m3, expected_initial_number,
                         2.0e-14),
          "projected ledger records initial particle number");
    check(close_relative(out.initial_energy_J_m3, expected_initial_energy,
                         2.0e-14),
          "projected ledger records actual candidate energy");
    check(close_relative(out.fluid_number_m3, expected_initial_number,
                         2.0e-14),
          "projected fluid number equals candidate number");
    check(close_relative(out.fluid_energy_J_m3, expected_fluid_energy,
                         2.0e-14),
          "projected fluid energy is exactly 1.5 kT N");
    check(close_relative(out.bath_energy_correction_J_m3, expected_correction,
                         2.0e-14, 1.0e-35),
          "signed bath correction equals actual candidate energy minus fluid energy");
    check(out.bath_energy_correction_J_m3 != 0.0,
          "fine midpoint candidate retains a signed correction");
    check(out.distribution_L1 <= max_l1 &&
              out.relative_mean_energy_error <= max_mean_error,
          "fine candidate diagnostics satisfy both projection tolerances");
    check(close_relative(
              out.relative_mean_energy_error,
              std::abs(old_energy /
                            (1.5L * static_cast<long double>(grid.kT) *
                             old_number) -
                        1.0L),
              2.0e-12, 1.0e-18),
          "projected mean-error diagnostic uses actual candidate energy");
    check(close_relative(out.fluid_number_m3 + out.remaining_number_m3 -
                             out.initial_number_m3,
                         0.0, 2.0e-12, 1.0e-20),
          "projected particle N closure holds");
    check(close_relative(out.fluid_energy_J_m3 +
                             out.bath_energy_correction_J_m3 +
                             out.remaining_energy_J_m3 -
                             out.initial_energy_J_m3,
                         0.0, 2.0e-12, 1.0e-30),
          "projected energy closure includes the signed bath correction once");
    check_handoff_finite(out, "projected handoff");
}

void test_hot_cold_rejection(const FineGrid& grid) {
    constexpr double full_number = 1.0e20;
    std::vector<double> old(grid.probability.size());
    for (std::size_t i = 0; i < old.size(); ++i) {
        old[i] = full_number * grid.probability[i];
    }
    const std::array<double, 2> target_temperatures{{0.7 * grid.kT,
                                                      1.4 * grid.kT}};
    for (double target : target_temperatures) {
        std::vector<double> trial(old.size(), -1.0);
        int projected = -1;
        fusion_handoff_ledger_v1 out;
        fill_handoff(out, -1.0);
        const int status = fusion_c_maxwellian_handoff_trial(
            FineGrid::cells, target, 1.0e-3, 1.0e-4, grid.edges.data(),
            old.data(), trial.data(), &projected, &out);
        check(status == PB11_STATUS_OK,
              "hot/cold nonqualifying handoff returns OK");
        check_rejected_handoff(out, old, trial, projected,
                               target < grid.kT ? "cold candidate" :
                                                  "hot candidate");
        check(out.distribution_L1 > 1.0e-3 ||
                  out.relative_mean_energy_error > 1.0e-4,
              "hot/cold candidate fails at least one explicit tolerance");
    }
}

void test_sparse_high_energy_tail(const FineGrid& grid) {
    constexpr double full_number = 1.0e20;
    constexpr double epsilon = 1.0e-3;
    const int removed = 100;
    const int added = FineGrid::cells - 1;
    std::vector<double> candidate_probability = grid.probability;
    check(candidate_probability[static_cast<std::size_t>(removed)] > epsilon,
          "tail test has enough central probability to remove");
    candidate_probability[static_cast<std::size_t>(removed)] -= epsilon;
    candidate_probability[static_cast<std::size_t>(added)] += epsilon;
    std::vector<double> old(candidate_probability.size());
    for (std::size_t i = 0; i < old.size(); ++i) {
        old[i] = full_number * candidate_probability[i];
    }
    std::vector<double> trial(old.size(), -1.0);
    int projected = -1;
    fusion_handoff_ledger_v1 out;
    fill_handoff(out, -1.0);
    const double max_l1 = 3.0e-3;
    const double max_mean_error = 1.0e-3;
    const int status = fusion_c_maxwellian_handoff_trial(
        FineGrid::cells, grid.kT, max_l1, max_mean_error, grid.edges.data(),
        old.data(), trial.data(), &projected, &out);
    check(status == PB11_STATUS_OK,
          "sparse high-energy-tail handoff returns OK");
    check(out.distribution_L1 < max_l1,
          "sparse high-energy tail stays below the L1 tolerance");
    check(out.relative_mean_energy_error > max_mean_error,
          "sparse high-energy tail exceeds the independent mean-energy tolerance");
    check_rejected_handoff(out, old, trial, projected,
                           "sparse high-energy-tail candidate");
}

void test_zero_candidate_and_rollback(const FineGrid& grid) {
    std::vector<double> zero_old(grid.probability.size(), 0.0);
    std::vector<double> zero_trial(zero_old.size(), -1.0);
    int projected = -1;
    fusion_handoff_ledger_v1 zero_out;
    fill_handoff(zero_out, -1.0);
    const int zero_status = fusion_c_maxwellian_handoff_trial(
        FineGrid::cells, grid.kT, 0.0, 0.0, grid.edges.data(), zero_old.data(),
        zero_trial.data(), &projected, &zero_out);
    check(zero_status == PB11_STATUS_OK && projected == 0,
          "zero-N candidate is accepted as a nonprojection");
    check(std::all_of(zero_trial.begin(), zero_trial.end(),
                      [](double value) { return value == 0.0; }),
          "zero-N candidate trial is zero");
    check(zero_out.initial_number_m3 == 0.0 &&
              zero_out.initial_energy_J_m3 == 0.0 &&
              zero_out.remaining_number_m3 == 0.0 &&
              zero_out.remaining_energy_J_m3 == 0.0 &&
              zero_out.fluid_number_m3 == 0.0 &&
              zero_out.fluid_energy_J_m3 == 0.0 &&
              zero_out.bath_energy_correction_J_m3 == 0.0 &&
              zero_out.distribution_L1 == 0.0 &&
              zero_out.relative_mean_energy_error == 0.0,
          "zero-N candidate has zero inventory and source fields");

    constexpr double full_number = 1.0e20;
    std::vector<double> old(grid.probability.size());
    for (std::size_t i = 0; i < old.size(); ++i) {
        old[i] = full_number * grid.probability[i];
    }
    const std::vector<double> old_copy = old;
    std::vector<double> first_trial(old.size(), -1.0);
    std::vector<double> second_trial(old.size(), -1.0);
    int first_projected = -1;
    int second_projected = -1;
    fusion_handoff_ledger_v1 first_out;
    fusion_handoff_ledger_v1 second_out;
    fill_handoff(first_out, -1.0);
    fill_handoff(second_out, -1.0);
    const int first_status = fusion_c_maxwellian_handoff_trial(
        FineGrid::cells, grid.kT, 1.0e-3, 1.0e-4, grid.edges.data(), old.data(),
        first_trial.data(), &first_projected, &first_out);
    const int second_status = fusion_c_maxwellian_handoff_trial(
        FineGrid::cells, grid.kT, 1.0e-3, 1.0e-4, grid.edges.data(), old.data(),
        second_trial.data(), &second_projected, &second_out);
    check(first_status == PB11_STATUS_OK && second_status == PB11_STATUS_OK &&
              first_projected == second_projected && first_trial == second_trial &&
              same_handoff(first_out, second_out),
          "repeated handoff trials are identical for rollback");
    check(old == old_copy, "repeated rollback trials leave old state unchanged");
}

void test_grid_errors() {
    const double kT = 2.0 * kJoulesPerKeV;
    const double valid_edges[3]{0.4 * kT, 1.0 * kT, 2.0 * kT};
    double probability[2]{-1.0, -1.0};
    fusion_maxwellian_grid_v1 out;

    const double bad_edges[3]{0.4 * kT, 1.0 * kT, 0.9 * kT};
    fill_grid(out, -1.0);
    int status = fusion_c_maxwellian_energy_grid(
        2, kT, bad_edges, probability, &out);
    check(status == PB11_STATUS_OUT_OF_RANGE && probability[0] == 0.0 &&
              probability[1] == 0.0 && zero_grid(out),
          "non-increasing grid edges reject and clear outputs");

    fill_grid(out, -1.0);
    probability[0] = probability[1] = -1.0;
    status = fusion_c_maxwellian_energy_grid(
        2, std::numeric_limits<double>::quiet_NaN(), valid_edges, probability,
        &out);
    check(status == PB11_STATUS_INVALID_ARGUMENT && probability[0] == 0.0 &&
              probability[1] == 0.0 && zero_grid(out),
          "nonfinite kT rejects and clears grid outputs");

    fill_grid(out, -1.0);
    probability[0] = probability[1] = -1.0;
    status = fusion_c_maxwellian_energy_grid(
        2, 0.0, valid_edges, probability, &out);
    check(status == PB11_STATUS_OUT_OF_RANGE && probability[0] == 0.0 &&
              probability[1] == 0.0 && zero_grid(out),
          "nonpositive kT rejects and clears grid outputs");

    fill_grid(out, -1.0);
    probability[0] = probability[1] = -1.0;
    status = fusion_c_maxwellian_energy_grid(2, kT, valid_edges, nullptr, &out);
    check(status == PB11_STATUS_NULL_OUTPUT && zero_grid(out),
          "null grid probability output is reported");

    fill_grid(out, -1.0);
    probability[0] = probability[1] = -1.0;
    status = fusion_c_maxwellian_energy_grid(2, kT, valid_edges, probability,
                                             nullptr);
    check(status == PB11_STATUS_NULL_OUTPUT && probability[0] == 0.0 &&
              probability[1] == 0.0,
          "null grid diagnostics output is reported and q clears");

    status = fusion_c_maxwellian_energy_grid(0, kT, nullptr, nullptr, nullptr);
    check(status == PB11_STATUS_INVALID_ARGUMENT,
          "zero grid extent is rejected");
}

void test_handoff_errors() {
    const double kT = 2.0 * kJoulesPerKeV;
    const double edges[2]{0.4 * kT, 2.0 * kT};
    const double old[1]{1.0e20};
    double trial[1]{-1.0};
    int projected = -1;
    fusion_handoff_ledger_v1 out;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    fill_handoff(out, -1.0);
    {
        const int status = fusion_c_maxwellian_handoff_trial(
            1, kT, 0.1, 0.1, edges, &nan, trial, &projected, &out);
        check(status == PB11_STATUS_INVALID_ARGUMENT && trial[0] == 0.0 &&
                  projected == 0 && zero_handoff(out),
              "nonfinite old candidate rejects and clears handoff outputs");
    }

    const double negative_old[1]{-1.0};
    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int negative_status = fusion_c_maxwellian_handoff_trial(
        1, kT, 0.1, 0.1, edges, negative_old, trial, &projected, &out);
    check(negative_status == PB11_STATUS_INVALID_ARGUMENT && trial[0] == 0.0 &&
              projected == 0 && zero_handoff(out),
          "negative old candidate rejects without clipping and clears outputs");

    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int bad_tolerance_status = fusion_c_maxwellian_handoff_trial(
        1, kT, -0.1, 0.1, edges, old, trial, &projected, &out);
    check(bad_tolerance_status == PB11_STATUS_OUT_OF_RANGE &&
              trial[0] == 0.0 && projected == 0 && zero_handoff(out),
          "negative L1 tolerance rejects and clears outputs");

    const double nonincreasing_edges[2]{2.0 * kT, 1.0 * kT};
    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int bad_extent_status = fusion_c_maxwellian_handoff_trial(
        1, kT, 0.1, 0.1, nonincreasing_edges, old, trial, &projected, &out);
    check(bad_extent_status == PB11_STATUS_OUT_OF_RANGE && trial[0] == 0.0 &&
              projected == 0 && zero_handoff(out),
          "non-increasing handoff edges reject and clear outputs");

    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int null_trial_status = fusion_c_maxwellian_handoff_trial(
        1, kT, 0.1, 0.1, edges, old, nullptr, &projected, &out);
    check(null_trial_status == PB11_STATUS_NULL_OUTPUT && projected == 0 &&
              zero_handoff(out),
          "null handoff trial output is reported");

    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int null_projected_status = fusion_c_maxwellian_handoff_trial(
        1, kT, 0.1, 0.1, edges, old, trial, nullptr, &out);
    check(null_projected_status == PB11_STATUS_NULL_OUTPUT && trial[0] == 0.0 &&
              zero_handoff(out),
          "null projected output is reported and trial clears");

    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int null_ledger_status = fusion_c_maxwellian_handoff_trial(
        1, kT, 0.1, 0.1, edges, old, trial, &projected, nullptr);
    check(null_ledger_status == PB11_STATUS_NULL_OUTPUT && trial[0] == 0.0 &&
              projected == 0,
          "null handoff ledger output is reported and trial clears");

    const double huge = std::numeric_limits<double>::max();
    const double huge_edges[2]{0.5 * huge, 0.75 * huge};
    const double huge_old[1]{huge};
    fill_handoff(out, -1.0);
    trial[0] = -1.0;
    projected = -1;
    const int overflow_status = fusion_c_maxwellian_handoff_trial(
        1, 1.0e308, 1.0, 1.0, huge_edges, huge_old, trial, &projected, &out);
    check(overflow_status == PB11_STATUS_NUMERICAL_FAILURE &&
              trial[0] == 0.0 && projected == 0 && zero_handoff(out),
          "unrepresentable candidate energy rejects and clears outputs");

    const int invalid_extent_status = fusion_c_maxwellian_handoff_trial(
        0, kT, 0.1, 0.1, nullptr, nullptr, nullptr, nullptr, nullptr);
    check(invalid_extent_status == PB11_STATUS_INVALID_ARGUMENT,
          "zero handoff extent is rejected");
}

void test_both_correction_signs() {
    const double T=1.602176634e-15, old[1]={2.e19};
    for(double shift : {-0.01,0.01}) {
        double edges[2]={(1.+shift)*T,(2.+shift)*T}, trial[1];
        fusion_handoff_ledger_v1 out{};int projected=0;
        const int status=fusion_c_maxwellian_handoff_trial(1,T,2.,.01,
            edges,old,trial,&projected,&out);
        check(status==PB11_STATUS_OK && projected==1,
              "explicit broad shape bound permits signed-correction fixture");
        check(shift*out.bath_energy_correction_J_m3>0,
              "bath correction preserves heating/cooling sign");
        check(close_relative(out.bath_energy_correction_J_m3,shift*T*old[0],1e-12),
              "signed correction equals the candidate's excess/deficit energy");
    }
}

}  // namespace

int main() {
    test_gamma32_grid_and_truncation();

    FineGrid grid;
    if (make_fine_grid(grid)) {
        test_fine_projection_and_signed_correction(grid);
        test_hot_cold_rejection(grid);
        test_sparse_high_energy_tail(grid);
        test_zero_candidate_and_rollback(grid);
    }
    test_both_correction_signs();
    test_grid_errors();
    test_handoff_errors();

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All fusion-handoff tests passed\n";
    return 0;
}
