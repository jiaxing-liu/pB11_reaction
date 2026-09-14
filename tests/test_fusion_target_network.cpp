#include "fusion_target_network.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TestRunner {
  int failures = 0;

  void check(bool condition, const std::string& what) {
    if (!condition) {
      ++failures;
      std::cerr << "FAIL: " << what << '\n';
    }
  }

  void close(double actual, double expected, const std::string& what,
             double relative = 5e-12, double absolute = 5e-13) {
    const double scale = std::max(std::abs(actual), std::abs(expected));
    const bool ok = std::isfinite(actual) && std::isfinite(expected) &&
                    std::abs(actual - expected) <=
                        absolute + relative * scale;
    if (!ok) {
      std::ostringstream message;
      message << what << " (actual=" << std::setprecision(17) << actual
              << ", expected=" << expected << ')';
      check(false, message.str());
    }
  }

  void zero(double value, const std::string& what) {
    check(value == 0.0, what);
  }
};

struct Inputs {
  double dt = 1.0;
  std::vector<double> fast_energy;
  std::vector<double> old_fast;
  std::vector<double> target_number;
  std::vector<double> target_energy;
  std::vector<fusion_target_network_edge_v1> edges;
};

struct Result {
  int rc = 0;
  std::vector<double> fast;
  std::vector<double> target;
  std::vector<double> target_energy;
  std::vector<double> events;
  fusion_target_network_v1 out{};
};

double* data_or_null(std::vector<double>& values) {
  return values.empty() ? nullptr : values.data();
}

const double* data_or_null(const std::vector<double>& values) {
  return values.empty() ? nullptr : values.data();
}

const fusion_target_network_edge_v1* data_or_null(
    const std::vector<fusion_target_network_edge_v1>& values) {
  return values.empty() ? nullptr : values.data();
}

Result run(const Inputs& input, double sentinel = 0.0) {
  Result result;
  result.fast.assign(input.old_fast.size(), sentinel);
  result.target.assign(input.target_number.size(), sentinel);
  result.target_energy.assign(input.target_energy.size(), sentinel);
  result.events.assign(input.edges.size(), sentinel);
  result.out.reactions_m3 = sentinel;
  result.out.removed_fast_energy_J_m3 = sentinel;
  result.out.removed_target_energy_J_m3 = sentinel;
  result.out.max_fast_number_relative_residual = sentinel;
  result.out.max_target_number_relative_residual = sentinel;
  result.out.max_target_energy_relative_residual = sentinel;
  result.out.iterations = static_cast<int>(sentinel);
  result.rc = fusion_c_target_network_trial(
      static_cast<int>(input.old_fast.size()),
      static_cast<int>(input.target_number.size()),
      static_cast<int>(input.edges.size()), input.dt,
      data_or_null(input.fast_energy), data_or_null(input.old_fast),
      data_or_null(input.target_number), data_or_null(input.target_energy),
      data_or_null(input.edges), data_or_null(result.fast),
      data_or_null(result.target), data_or_null(result.target_energy),
      data_or_null(result.events), &result.out);
  return result;
}

void expect_success(TestRunner& tests, const Result& result,
                    const std::string& name) {
  tests.check(result.rc == 0, name + " returns success");
}

void expect_diagnostics(TestRunner& tests, const fusion_target_network_v1& out,
                        const std::string& name) {
  tests.check(std::isfinite(out.reactions_m3) && out.reactions_m3 >= 0.0,
              name + " reaction total is finite and nonnegative");
  tests.check(std::isfinite(out.removed_fast_energy_J_m3) &&
                  out.removed_fast_energy_J_m3 >= 0.0,
              name + " fast energy ledger is finite and nonnegative");
  tests.check(std::isfinite(out.removed_target_energy_J_m3) &&
                  out.removed_target_energy_J_m3 >= 0.0,
              name + " target energy ledger is finite and nonnegative");
  tests.check(std::isfinite(out.max_fast_number_relative_residual) &&
                  out.max_fast_number_relative_residual <= 2.0e-12,
              name + " fast residual meets the numerical gate");
  tests.check(std::isfinite(out.max_target_number_relative_residual) &&
                  out.max_target_number_relative_residual <= 2.0e-12,
              name + " target residual meets the numerical gate");
  tests.check(std::isfinite(out.max_target_energy_relative_residual) &&
                  out.max_target_energy_relative_residual <= 2.0e-12,
              name + " target-energy residual meets the numerical gate");
  tests.check(out.iterations >= 0 && out.iterations <= 512,
              name + " iteration count is bounded");
}

void expect_all_zero(TestRunner& tests, const Result& result,
                     const std::string& name) {
  for (std::size_t i = 0; i < result.fast.size(); ++i)
    tests.zero(result.fast[i], name + " fast output[" + std::to_string(i) + "] cleared");
  for (std::size_t i = 0; i < result.target.size(); ++i)
    tests.zero(result.target[i], name + " target output[" + std::to_string(i) + "] cleared");
  for (std::size_t i = 0; i < result.target_energy.size(); ++i)
    tests.zero(result.target_energy[i],
               name + " target-energy output[" + std::to_string(i) + "] cleared");
  for (std::size_t i = 0; i < result.events.size(); ++i)
    tests.zero(result.events[i], name + " event output[" + std::to_string(i) + "] cleared");
  tests.zero(result.out.reactions_m3, name + " total reactions cleared");
  tests.zero(result.out.removed_fast_energy_J_m3,
             name + " removed fast energy cleared");
  tests.zero(result.out.removed_target_energy_J_m3,
             name + " removed target energy cleared");
  tests.zero(result.out.max_fast_number_relative_residual,
             name + " fast residual cleared");
  tests.zero(result.out.max_target_number_relative_residual,
             name + " target residual cleared");
  tests.zero(result.out.max_target_energy_relative_residual,
             name + " target-energy residual cleared");
  tests.check(result.out.iterations == 0, name + " iteration count cleared");
}

void test_golden_root(TestRunner& tests) {
  Inputs input;
  input.fast_energy = {5.0};
  input.old_fast = {1.0};
  input.target_number = {1.0};
  input.target_energy = {10.0};
  input.edges = {{0, 0, 1.0, 0.0}};

  const Result result = run(input);
  expect_success(tests, result, "one-fast/one-target golden root");
  if (result.rc != 0) return;

  const double root = (std::sqrt(5.0) - 1.0) / 2.0;
  const double reaction = 1.0 - root;
  tests.close(result.fast[0], root, "golden fast inventory");
  tests.close(result.target[0], root, "golden target inventory");
  tests.close(result.target_energy[0], 10.0, "zero-M target energy");
  tests.close(result.events[0], reaction, "golden edge reaction");
  tests.close(result.out.reactions_m3, reaction, "golden total reaction");
  tests.close(result.out.removed_fast_energy_J_m3, 5.0 * reaction,
              "golden removed fast energy");
  tests.close(result.out.removed_target_energy_J_m3, 0.0,
              "golden removed target energy");
  expect_diagnostics(tests, result.out, "golden root");
}

void test_extreme_stiff_scale(TestRunner& tests) {
  Inputs input;
  input.fast_energy = {3.0};
  input.old_fast = {1.0};
  input.target_number = {1.0};
  input.target_energy = {10.0};
  // M/K is one joule per event even though K is deliberately very large.
  input.edges = {{0, 0, 1e100, 1e100}};

  const Result result = run(input);
  expect_success(tests, result, "extreme stiff one-bin root");
  if (result.rc != 0) return;

  const double root = 2.0 / (1.0 + std::sqrt(1.0 + 4.0e100));
  const double reaction = 1.0 - root;
  tests.close(result.fast[0], root, "extreme stiff fast root", 1e-10, 1e-60);
  tests.close(result.target[0], root, "extreme stiff target root", 1e-10,
              1e-60);
  tests.close(result.events[0], reaction, "extreme stiff event", 1e-10,
              1e-12);
  tests.close(result.target_energy[0], 10.0 - reaction,
              "extreme stiff target energy", 1e-10, 1e-12);
  tests.close(result.out.removed_fast_energy_J_m3, 3.0 * reaction,
              "extreme stiff fast energy ledger", 1e-10, 1e-12);
  tests.close(result.out.removed_target_energy_J_m3, reaction,
              "extreme stiff M/K target energy ledger", 1e-10, 1e-12);
  expect_diagnostics(tests, result.out, "extreme stiff root");
}

void test_parallel_channels(TestRunner& tests) {
  Inputs input;
  input.fast_energy = {7.0};
  input.old_fast = {1.0};
  input.target_number = {1.0};
  input.target_energy = {8.0};
  // M/K is 0.75 J for both channels; using M directly would debit the
  // 3-rate edge three times too much.
  input.edges = {{0, 0, 1.0, 0.75}, {0, 0, 3.0, 2.25}};

  const Result result = run(input);
  expect_success(tests, result, "parallel 1:3 channels");
  if (result.rc != 0) return;

  const double root = (std::sqrt(17.0) - 1.0) / 8.0;
  const double total = 1.0 - root;
  tests.close(result.fast[0], root, "parallel shared fast inventory");
  tests.close(result.target[0], root, "parallel shared target inventory");
  tests.close(result.events[0], total / 4.0, "parallel 1-rate event");
  tests.close(result.events[1], 3.0 * total / 4.0, "parallel 3-rate event");
  tests.close(result.out.reactions_m3, total, "parallel total event");
  tests.close(result.out.removed_fast_energy_J_m3, 7.0 * total,
              "parallel removed fast energy");
  tests.close(result.out.removed_target_energy_J_m3, 0.75 * total,
              "parallel M/K target energy");
  tests.close(result.target_energy[0], 8.0 - 0.75 * total,
              "parallel target energy");
  expect_diagnostics(tests, result.out, "parallel channels");
}

void test_competing_targets(TestRunner& tests) {
  Inputs input;
  input.fast_energy = {2.0};
  input.old_fast = {2.0};
  input.target_number = {1.0, 1.0};
  input.target_energy = {3.0, 3.0};
  input.edges = {{0, 0, 1.0, 0.0}, {0, 1, 1.0, 0.0}};

  const Result result = run(input);
  expect_success(tests, result, "one-fast/two-target symmetry");
  if (result.rc != 0) return;

  // f = 2/(1+b0+b1), b_j = 1/(1+f), hence f=1 and b_j=1/2.
  tests.close(result.fast[0], 1.0, "competing-target fast root");
  tests.close(result.target[0], 0.5, "first competing target root");
  tests.close(result.target[1], 0.5, "second competing target root");
  tests.close(result.events[0], 0.5, "first competing event");
  tests.close(result.events[1], 0.5, "second competing event");
  tests.close(result.out.reactions_m3, 1.0, "competing total event");
  tests.close(result.out.removed_fast_energy_J_m3, 2.0,
              "competing removed fast energy");
  expect_diagnostics(tests, result.out, "competing targets");
}

struct ExactNetwork {
  Inputs input;
  std::vector<double> final_fast;
  std::vector<double> final_target;
  std::vector<double> final_target_energy;
  std::vector<double> expected_events;
};

ExactNetwork make_exact_network() {
  ExactNetwork network;
  network.input.dt = 0.75;
  network.input.fast_energy = {2.0, 5.0, 3.5};
  network.final_fast = {0.8, 1.3, 0.45};
  network.final_target = {1.1, 0.7, 0.9};
  network.final_target_energy = {5.0, 6.0, 7.0};
  network.input.edges = {
      {0, 0, 0.7, 0.7 * 0.4},
      {0, 2, 1.1, 1.1 * 0.9},
      {1, 0, 0.3, 0.3 * 0.2},
      {1, 1, 0.8, 0.8 * 0.8},
      {2, 1, 1.4, 1.4 * 0.6},
      {2, 2, 0.5, 0.5 * 1.1},
  };

  network.input.old_fast.assign(network.final_fast.size(), 0.0);
  network.input.target_number.assign(network.final_target.size(), 0.0);
  network.input.target_energy = network.final_target_energy;
  network.expected_events.assign(network.input.edges.size(), 0.0);

  for (std::size_t i = 0; i < network.final_fast.size(); ++i) {
    double factor = 1.0;
    for (const auto& edge : network.input.edges) {
      if (edge.fast_index == static_cast<int>(i))
        factor += network.input.dt * edge.reactivity_m3_s *
                  network.final_target[edge.target_index];
    }
    network.input.old_fast[i] = network.final_fast[i] * factor;
  }

  network.input.target_number = network.final_target;
  for (std::size_t e = 0; e < network.input.edges.size(); ++e) {
    const auto& edge = network.input.edges[e];
    const double loss = network.input.dt * edge.reactivity_m3_s *
                        network.final_fast[edge.fast_index] *
                        network.final_target[edge.target_index];
    const double target_energy_per_event =
        edge.target_energy_reactivity_J_m3_s / edge.reactivity_m3_s;
    network.expected_events[e] = loss;
    network.input.target_number[edge.target_index] += loss;
    network.input.target_energy[edge.target_index] +=
        loss * target_energy_per_event;
  }
  return network;
}

void test_constructed_exact_solution(TestRunner& tests,
                                     const ExactNetwork& network) {
  const Result result = run(network.input);
  expect_success(tests, result, "constructed multi-bin exact solution");
  if (result.rc != 0) return;

  for (std::size_t i = 0; i < network.final_fast.size(); ++i)
    tests.close(result.fast[i], network.final_fast[i],
                "constructed fast inventory[" + std::to_string(i) + "]",
                3e-11, 2e-12);
  for (std::size_t j = 0; j < network.final_target.size(); ++j) {
    tests.close(result.target[j], network.final_target[j],
                "constructed target inventory[" + std::to_string(j) + "]",
                3e-11, 2e-12);
    tests.close(result.target_energy[j], network.final_target_energy[j],
                "constructed target energy[" + std::to_string(j) + "]",
                3e-11, 2e-12);
  }

  double total = 0.0;
  double removed_fast = 0.0;
  double removed_target = 0.0;
  for (std::size_t e = 0; e < network.expected_events.size(); ++e) {
    const auto& edge = network.input.edges[e];
    const double event = network.expected_events[e];
    total += event;
    removed_fast += event * network.input.fast_energy[edge.fast_index];
    removed_target +=
        event * edge.target_energy_reactivity_J_m3_s /
        edge.reactivity_m3_s;
    tests.close(result.events[e], event,
                "constructed edge event[" + std::to_string(e) + "]",
                4e-11, 2e-12);
  }
  tests.close(result.out.reactions_m3, total, "constructed total events",
              4e-11, 2e-12);
  tests.close(result.out.removed_fast_energy_J_m3, removed_fast,
              "constructed fast energy ledger", 4e-11, 2e-12);
  tests.close(result.out.removed_target_energy_J_m3, removed_target,
              "constructed target M/K ledger", 4e-11, 2e-12);
  expect_diagnostics(tests, result.out, "constructed solution");
}

int find_index(const std::vector<int>& permutation, int old_index) {
  for (std::size_t i = 0; i < permutation.size(); ++i)
    if (permutation[i] == old_index) return static_cast<int>(i);
  return -1;
}

Inputs permute_inputs(const Inputs& source, const std::vector<int>& fast_perm,
                      const std::vector<int>& target_perm,
                      const std::vector<int>& edge_order) {
  Inputs result;
  result.dt = source.dt;
  result.fast_energy.resize(fast_perm.size());
  result.old_fast.resize(fast_perm.size());
  result.target_number.resize(target_perm.size());
  result.target_energy.resize(target_perm.size());
  for (std::size_t i = 0; i < fast_perm.size(); ++i) {
    result.fast_energy[i] = source.fast_energy[fast_perm[i]];
    result.old_fast[i] = source.old_fast[fast_perm[i]];
  }
  for (std::size_t j = 0; j < target_perm.size(); ++j) {
    result.target_number[j] = source.target_number[target_perm[j]];
    result.target_energy[j] = source.target_energy[target_perm[j]];
  }
  result.edges.reserve(edge_order.size());
  for (int old_edge_index : edge_order) {
    fusion_target_network_edge_v1 edge = source.edges[old_edge_index];
    edge.fast_index = find_index(fast_perm, edge.fast_index);
    edge.target_index = find_index(target_perm, edge.target_index);
    result.edges.push_back(edge);
  }
  return result;
}

void test_permutation_invariance(TestRunner& tests,
                                 const ExactNetwork& network) {
  const Result original = run(network.input);
  const std::vector<int> fast_perm = {2, 0, 1};
  const std::vector<int> target_perm = {1, 2, 0};
  const std::vector<int> edge_order = {5, 2, 0, 4, 1, 3};
  const Inputs permuted_inputs =
      permute_inputs(network.input, fast_perm, target_perm, edge_order);
  const Result permuted = run(permuted_inputs);

  expect_success(tests, original, "permutation reference call");
  expect_success(tests, permuted, "permutation call");
  if (original.rc != 0 || permuted.rc != 0) return;

  for (std::size_t i = 0; i < fast_perm.size(); ++i)
    tests.close(permuted.fast[i], original.fast[fast_perm[i]],
                "permuted fast inventory[" + std::to_string(i) + "]",
                5e-11, 2e-12);
  for (std::size_t j = 0; j < target_perm.size(); ++j) {
    tests.close(permuted.target[j], original.target[target_perm[j]],
                "permuted target inventory[" + std::to_string(j) + "]",
                5e-11, 2e-12);
    tests.close(permuted.target_energy[j],
                original.target_energy[target_perm[j]],
                "permuted target energy[" + std::to_string(j) + "]",
                5e-11, 2e-12);
  }
  for (std::size_t e = 0; e < edge_order.size(); ++e)
    tests.close(permuted.events[e], original.events[edge_order[e]],
                "permuted event[" + std::to_string(e) + "]", 5e-11,
                2e-12);

  tests.close(permuted.out.reactions_m3, original.out.reactions_m3,
              "permuted total event", 5e-11, 2e-12);
  tests.close(permuted.out.removed_fast_energy_J_m3,
              original.out.removed_fast_energy_J_m3,
              "permuted fast energy ledger", 5e-11, 2e-12);
  tests.close(permuted.out.removed_target_energy_J_m3,
              original.out.removed_target_energy_J_m3,
              "permuted target energy ledger", 5e-11, 2e-12);
}

void test_weak_event(TestRunner& tests) {
  Inputs input;
  // The relative depletion is about 1e-17, below one ulp at 1e150.  The
  // explicit event ledger must still retain a positive reaction amount.
  input.fast_energy = {1.25};
  input.old_fast = {1e150};
  input.target_number = {1e150};
  input.target_energy = {1e150};
  input.edges = {{0, 0, 1e-167, 0.0}};

  const Result result = run(input);
  expect_success(tests, result, "weak event");
  if (result.rc != 0) return;

  tests.check(result.fast[0] == input.old_fast[0],
              "weak event fast inventory difference rounds away");
  tests.check(result.target[0] == input.target_number[0],
              "weak event target inventory difference rounds away");
  tests.check(std::isfinite(result.events[0]) && result.events[0] > 0.0,
              "weak event keeps a positive per-edge event");
  tests.check(std::isfinite(result.out.reactions_m3) &&
                  result.out.reactions_m3 == result.events[0] &&
                  result.out.reactions_m3 > 0.0,
              "weak event total uses explicit edge ledger");
  tests.check(result.out.removed_fast_energy_J_m3 > 0.0,
              "weak event retains fast-energy debit");
}

void test_energy_overdraw_clears_outputs(TestRunner& tests) {
  Inputs input;
  input.fast_energy = {2.0};
  input.old_fast = {1.0};
  input.target_number = {1.0};
  input.target_energy = {0.1};
  input.edges = {{0, 0, 1.0, 2.0}};

  const Result result = run(input, 42.0);
  tests.check(result.rc != 0, "target-energy overdraw fails");
  expect_all_zero(tests, result, "target-energy overdraw");
}

void test_retry_is_deterministic_and_inputs_immutable(TestRunner& tests,
                                                        const Inputs& input) {
  const Inputs before = input;
  const Result first = run(input);
  const Result second = run(input);
  expect_success(tests, first, "retry first call");
  expect_success(tests, second, "retry second call");
  if (first.rc != 0 || second.rc != 0) return;

  tests.check(std::memcmp(first.fast.data(), second.fast.data(),
                          first.fast.size() * sizeof(double)) == 0,
              "retry fast output is deterministic");
  tests.check(std::memcmp(first.target.data(), second.target.data(),
                          first.target.size() * sizeof(double)) == 0,
              "retry target output is deterministic");
  tests.check(std::memcmp(first.target_energy.data(), second.target_energy.data(),
                          first.target_energy.size() * sizeof(double)) == 0,
              "retry target-energy output is deterministic");
  tests.check(std::memcmp(first.events.data(), second.events.data(),
                          first.events.size() * sizeof(double)) == 0,
              "retry event output is deterministic");
  tests.check(first.out.reactions_m3 == second.out.reactions_m3 &&
                  first.out.removed_fast_energy_J_m3 ==
                      second.out.removed_fast_energy_J_m3 &&
                  first.out.removed_target_energy_J_m3 ==
                      second.out.removed_target_energy_J_m3 &&
                  first.out.max_fast_number_relative_residual ==
                      second.out.max_fast_number_relative_residual &&
                  first.out.max_target_number_relative_residual ==
                      second.out.max_target_number_relative_residual &&
                  first.out.max_target_energy_relative_residual ==
                      second.out.max_target_energy_relative_residual &&
                  first.out.iterations == second.out.iterations,
              "retry diagnostics are deterministic");

  tests.check(std::memcmp(input.fast_energy.data(), before.fast_energy.data(),
                          input.fast_energy.size() * sizeof(double)) == 0,
              "retry leaves fast energies unchanged");
  tests.check(std::memcmp(input.old_fast.data(), before.old_fast.data(),
                          input.old_fast.size() * sizeof(double)) == 0,
              "retry leaves fast inventories unchanged");
  tests.check(std::memcmp(input.target_number.data(), before.target_number.data(),
                          input.target_number.size() * sizeof(double)) == 0,
              "retry leaves target inventories unchanged");
  tests.check(std::memcmp(input.target_energy.data(), before.target_energy.data(),
                          input.target_energy.size() * sizeof(double)) == 0,
              "retry leaves target energies unchanged");
  tests.check(std::memcmp(input.edges.data(), before.edges.data(),
                          input.edges.size() * sizeof(fusion_target_network_edge_v1)) == 0,
              "retry leaves edge coefficients unchanged");
}

void test_zero_edge_boundary(TestRunner& tests) {
  Inputs input;
  input.fast_energy = {2.0, 3.0};
  input.old_fast = {4.0, 5.0};
  input.target_number = {6.0, 7.0};
  input.target_energy = {8.0, 9.0};

  const Result result = run(input);
  expect_success(tests, result, "zero-edge boundary");
  if (result.rc != 0) return;
  tests.close(result.fast[0], 4.0, "zero-edge fast[0]");
  tests.close(result.fast[1], 5.0, "zero-edge fast[1]");
  tests.close(result.target[0], 6.0, "zero-edge target[0]");
  tests.close(result.target[1], 7.0, "zero-edge target[1]");
  tests.close(result.target_energy[0], 8.0, "zero-edge energy[0]");
  tests.close(result.target_energy[1], 9.0, "zero-edge energy[1]");
  tests.zero(result.out.reactions_m3, "zero-edge total reaction");
  tests.zero(result.out.removed_fast_energy_J_m3,
             "zero-edge removed fast energy");
  tests.zero(result.out.removed_target_energy_J_m3,
             "zero-edge removed target energy");

  Inputs zero_rate;
  zero_rate.fast_energy = {2.0};
  zero_rate.old_fast = {4.0};
  zero_rate.target_number = {6.0};
  zero_rate.target_energy = {8.0};
  zero_rate.edges = {{0, 0, 0.0, 0.0}};
  const Result zero_rate_result = run(zero_rate);
  expect_success(tests, zero_rate_result, "zero-reactivity edge");
  if (zero_rate_result.rc == 0) {
    tests.close(zero_rate_result.fast[0], 4.0,
                "zero-reactivity fast inventory");
    tests.close(zero_rate_result.target[0], 6.0,
                "zero-reactivity target inventory");
    tests.close(zero_rate_result.target_energy[0], 8.0,
                "zero-reactivity target energy");
    tests.zero(zero_rate_result.events[0], "zero-reactivity event");
    tests.zero(zero_rate_result.out.reactions_m3,
               "zero-reactivity total reaction");
  }
}

void test_invalid_arguments(TestRunner& tests) {
  const auto invalid_with_clear = [&](const std::string& name,
                                      double dt,
                                      const double* fast_energy,
                                      const double* old_fast,
                                      const double* target_number,
                                      const double* target_energy,
                                      const fusion_target_network_edge_v1* edges) {
    Inputs input;
    input.dt = dt;
    input.fast_energy = {1.0};
    input.old_fast = {1.0};
    input.target_number = {1.0};
    input.target_energy = {1.0};
    input.edges = {{0, 0, 1.0, 0.0}};
    Result result;
    result.fast.assign(1, 77.0);
    result.target.assign(1, 77.0);
    result.target_energy.assign(1, 77.0);
    result.events.assign(1, 77.0);
    result.out.reactions_m3 = 77.0;
    result.out.removed_fast_energy_J_m3 = 77.0;
    result.out.removed_target_energy_J_m3 = 77.0;
    result.out.max_fast_number_relative_residual = 77.0;
    result.out.max_target_number_relative_residual = 77.0;
    result.out.max_target_energy_relative_residual = 77.0;
    result.out.iterations = 77;
    result.rc = fusion_c_target_network_trial(
        1, 1, 1, dt, fast_energy, old_fast, target_number, target_energy,
        edges, result.fast.data(), result.target.data(),
        result.target_energy.data(), result.events.data(), &result.out);
    tests.check(result.rc != 0, name + " rejects invalid input");
    expect_all_zero(tests, result, name);
  };

  double energy = 1.0;
  double old = 1.0;
  double target = 1.0;
  double target_energy = 1.0;
  fusion_target_network_edge_v1 edge{0, 0, 1.0, 0.0};
  invalid_with_clear("zero dt", 0.0, &energy, &old, &target, &target_energy,
                    &edge);
  invalid_with_clear("NaN dt", std::numeric_limits<double>::quiet_NaN(),
                    &energy, &old, &target, &target_energy, &edge);

  fusion_target_network_edge_v1 bad_index{1, 0, 1.0, 0.0};
  invalid_with_clear("edge index", 1.0, &energy, &old, &target,
                    &target_energy, &bad_index);
  fusion_target_network_edge_v1 negative_rate{0, 0, -1.0, 0.0};
  invalid_with_clear("negative reactivity", 1.0, &energy, &old, &target,
                    &target_energy, &negative_rate);
  fusion_target_network_edge_v1 zero_with_m{0, 0, 0.0, 1.0};
  invalid_with_clear("M without K", 1.0, &energy, &old, &target,
                    &target_energy, &zero_with_m);

  double negative_old = -1.0;
  invalid_with_clear("negative fast inventory", 1.0, &energy, &negative_old,
                    &target, &target_energy, &edge);
  double negative_target = -1.0;
  invalid_with_clear("negative target inventory", 1.0, &energy, &old,
                    &negative_target, &target_energy, &edge);
  double negative_energy = -1.0;
  invalid_with_clear("negative target energy", 1.0, &energy, &old, &target,
                    &negative_energy, &edge);
  double nonzero_energy_without_target = 1.0;
  double zero_target = 0.0;
  invalid_with_clear("energy without target particles", 1.0, &energy, &old,
                    &zero_target, &nonzero_energy_without_target, &edge);
  invalid_with_clear("null fast energy", 1.0, nullptr, &old, &target,
                    &target_energy, &edge);
  invalid_with_clear("null edge array", 1.0, &energy, &old, &target,
                    &target_energy, nullptr);

  // Invalid dimensions are rejected as well.  The output arrays have valid
  // storage, but the contract only promises clearing for valid dimensions.
  double next_fast = 17.0;
  double next_target = 17.0;
  double next_energy = 17.0;
  double event = 17.0;
  fusion_target_network_v1 out{};
  out.reactions_m3 = 17.0;
  const int bad_dimensions = fusion_c_target_network_trial(
      0, 1, 1, 1.0, &energy, &old, &target, &target_energy, &edge,
      &next_fast, &next_target, &next_energy, &event, &out);
  tests.check(bad_dimensions != 0, "zero fast dimension is rejected");
}

}  // namespace

int main() {
  TestRunner tests;
  test_golden_root(tests);
  test_extreme_stiff_scale(tests);
  test_parallel_channels(tests);
  test_competing_targets(tests);

  const ExactNetwork exact = make_exact_network();
  test_constructed_exact_solution(tests, exact);
  test_permutation_invariance(tests, exact);
  test_weak_event(tests);
  test_energy_overdraw_clears_outputs(tests);
  test_retry_is_deterministic_and_inputs_immutable(tests, exact.input);
  test_zero_edge_boundary(tests);
  test_invalid_arguments(tests);

  if (tests.failures != 0) {
    std::cerr << tests.failures << " test assertion(s) failed\n";
    return 1;
  }
  std::cout << "fusion target network tests passed\n";
  return 0;
}
