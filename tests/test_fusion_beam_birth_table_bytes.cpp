#include "fusion_beam_birth_table.h"
#include "fusion_nuclear_data.h"
#include "pb11_c.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kev = 1.602176634e-16;
constexpr double mev = 1.602176634e-13;
constexpr uint64_t entry_cap = 12500000ULL;

void check(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

template <class T>
using owner =
    std::unique_ptr<T, decltype(&fusion_c_beam_birth_table_destroy)>;

fusion_beam_birth_options_v1 beam_options() {
  fusion_beam_birth_options_v1 o{};
  o.relative_max_J = 2.5 * mev;
  o.angular_max_exponent = 40;
  o.ground_state_q_J = 91.84 * kev;
  o.cutoff_J = .001 * mev;
  o.l1_fraction = .76;
  o.relative_phase = 0;
  o.narrow_peak_fraction = .051;
  o.continuum_peak_scale = 1;
  o.continuation = 1;
  o.pb_low = 0;
  o.remainder_policy = 0;
  o.broad_mode = 13;
  o.fsci_policy = 0;
  o.relative_order = 4;
  o.angular_order = 4;
  o.nq = 4;
  o.ncos = 4;
  return o;
}

fusion_birth_table_control_v1 table_control() {
  // Five direct samples and two stored endpoint knots are enough for the
  // deliberately broad, cheap DT fixture used by these persistence tests.
  return fusion_birth_table_control_v1{.5, .5, .5, .5, .1, .1, 8, 64, 0};
}

std::vector<double> ordinary_edges(int cells) {
  std::vector<double> edges(static_cast<size_t>(cells) + 1);
  for (int j = 0; j <= cells; ++j)
    edges[static_cast<size_t>(j)] = 25 * mev * j / cells;
  return edges;
}

struct Fixture {
  owner<fusion_beam_birth_table_v1> table;
  fusion_beam_birth_options_v1 options{};
  fusion_birth_table_control_v1 control{};
  std::vector<double> edges;
  int channel = 3;
  int slot = 0;
  double projectile_energy_J = 1.5 * mev;
  double lower_kT_J = .08 * kev;
  double upper_kT_J = .12 * kev;

  Fixture() : table(nullptr, fusion_c_beam_birth_table_destroy) {}
};

Fixture make_ordinary_fixture() {
  Fixture fixture;
  fixture.options = beam_options();
  fixture.control = table_control();
  fixture.edges = ordinary_edges(8);
  fusion_beam_birth_table_v1 *raw = nullptr;
  const int status = fusion_c_beam_birth_table_create(
      fixture.channel, fixture.slot, fixture.projectile_energy_J,
      fixture.lower_kT_J, fixture.upper_kT_J, &fixture.options,
      &fixture.control, static_cast<int>(fixture.edges.size() - 1),
      fixture.edges.data(), &raw);
  check(status == PB11_STATUS_OK && raw, "ordinary DT fixture creation");
  fixture.table = owner<fusion_beam_birth_table_v1>(
      raw, fusion_c_beam_birth_table_destroy);
  fusion_beam_birth_table_info_v1 info{};
  check(fusion_c_beam_birth_table_info(fixture.table.get(), &info) ==
            PB11_STATUS_OK,
        "ordinary DT fixture info");
  check(info.knots >= 2 && info.direct_evaluations >= 5 && info.cells == 8,
        "ordinary DT fixture metadata");
  return fixture;
}

std::vector<double> subnormal_edges(int cells) {
  std::vector<double> edges(static_cast<size_t>(cells) + 1);
  for (int j = 1; j <= cells / 3; ++j) {
    const double fraction = double(j - 1) / double(cells / 3 - 1);
    edges[static_cast<size_t>(j)] = 1e-10 * kev * std::pow(1e12, fraction);
  }
  for (int j = cells / 3 + 1; j <= cells; ++j)
    edges[static_cast<size_t>(j)] =
        (100 + 24900. * (j - cells / 3) / (cells - cells / 3)) * kev;
  return edges;
}

Fixture make_subnormal_fixture() {
  Fixture fixture;
  constexpr int cells = 800;
  fixture.options = beam_options();
  fixture.options.relative_order = 16;
  fixture.options.angular_order = 16;
  fixture.options.nq = fixture.options.ncos = 8;
  fixture.control = fusion_birth_table_control_v1{
      .01, .01, .01, .01, 1e-5, 1e-5, 2, 5, 0};
  fixture.edges = subnormal_edges(cells);
  fixture.slot = 1;
  fixture.projectile_energy_J = 2.1399859254691007e-13;
  fixture.lower_kT_J = 1.4401072590598843e-17;
  fixture.upper_kT_J = fixture.lower_kT_J * 1.001;
  fusion_beam_birth_table_v1 *raw = nullptr;
  const int status = fusion_c_beam_birth_table_create(
      fixture.channel, fixture.slot, fixture.projectile_energy_J,
      fixture.lower_kT_J, fixture.upper_kT_J, &fixture.options,
      &fixture.control, cells, fixture.edges.data(), &raw);
  check(status == PB11_STATUS_OK && raw,
        "subnormal DT fixture creation status=" + std::to_string(status));
  fixture.table = owner<fusion_beam_birth_table_v1>(
      raw, fusion_c_beam_birth_table_destroy);
  return fixture;
}

uint64_t load_u64(const std::vector<unsigned char> &bytes, size_t offset) {
  check(offset <= bytes.size() && bytes.size() - offset >= 8,
        "wire u64 offset");
  uint64_t value = 0;
  for (int j = 0; j < 8; ++j)
    value |= uint64_t(bytes[offset + static_cast<size_t>(j)]) << (8 * j);
  return value;
}

void store_u64(std::vector<unsigned char> &bytes, size_t offset,
               uint64_t value) {
  check(offset <= bytes.size() && bytes.size() - offset >= 8,
        "wire u64 store offset");
  for (int j = 0; j < 8; ++j)
    bytes[offset + static_cast<size_t>(j)] =
        static_cast<unsigned char>(value >> (8 * j));
}

double load_double(const std::vector<unsigned char> &bytes, size_t offset) {
  const uint64_t bits = load_u64(bytes, offset);
  double value = 0;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

void store_double(std::vector<unsigned char> &bytes, size_t offset,
                  double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(value));
  store_u64(bytes, offset, bits);
}

uint64_t checksum(const unsigned char *data, size_t length) {
  uint64_t value = 14695981039346656037ULL;
  for (size_t j = 0; j < length; ++j) {
    value ^= data[j];
    value *= 1099511628211ULL;
  }
  return value;
}

void rechecksum(std::vector<unsigned char> &bytes) {
  check(bytes.size() >= 8, "checksum fixture size");
  store_u64(bytes, bytes.size() - 8,
            checksum(bytes.data(), bytes.size() - 8));
}

struct SparseEntryWire {
  size_t index_offset = 0;
  size_t value_offset = 0;
};

struct KnotWire {
  size_t temperature_offset = 0;
  std::array<size_t, 31> coefficient_offsets{};
  size_t count_offset = 0;
  std::vector<SparseEntryWire> entries;
};

struct WireLayout {
  std::array<size_t, 42> metadata_offsets{};
  std::vector<std::string> metadata_names;
  size_t identity_offset = 24;
  size_t edges_offset = 0;
  std::vector<size_t> edge_offsets;
  std::vector<KnotWire> knots;
  size_t checksum_offset = 0;

  size_t metadata(const std::string &name) const {
    for (size_t j = 0; j < metadata_names.size(); ++j)
      if (metadata_names[j] == name) return metadata_offsets[j];
    throw std::runtime_error("missing wire metadata field: " + name);
  }
};

const std::vector<std::string> &metadata_names() {
  static const std::vector<std::string> names = {
      "projectile_energy_J",
      "lower_kT_J",
      "upper_kT_J",
      "max_validated_rate_error",
      "max_validated_debit_error",
      "max_validated_number_L1",
      "max_validated_energy_L1",
      "max_sampled_direct_rate_discrepancy",
      "max_sampled_direct_debit_discrepancy",
      "channel",
      "projectile_slot",
      "cells",
      "knots",
      "direct_evaluations",
      "spectral_entries_evaluated",
      "stored_spectral_entries",
      "source.relative_max_J",
      "source.angular_max_exponent",
      "source.ground_state_q_J",
      "source.cutoff_J",
      "source.l1_fraction",
      "source.relative_phase",
      "source.narrow_peak_fraction",
      "source.continuum_peak_scale",
      "source.continuation",
      "source.pb_low",
      "source.remainder_policy",
      "source.broad_mode",
      "source.fsci_policy",
      "source.relative_order",
      "source.angular_order",
      "source.nq",
      "source.ncos",
      "control.max_rate_error",
      "control.max_debit_error",
      "control.max_number_L1",
      "control.max_energy_L1",
      "control.max_direct_rate_discrepancy",
      "control.max_direct_debit_discrepancy",
      "control.max_knots",
      "control.max_evaluations",
      "control.max_depth"};
  return names;
}

WireLayout locate_wire(const std::vector<unsigned char> &bytes,
                       const fusion_beam_birth_table_info_v1 &info) {
  check(bytes.size() >= 96, "packed fixture minimum size");
  WireLayout layout;
  layout.metadata_names = metadata_names();
  check(layout.metadata_names.size() == layout.metadata_offsets.size(),
        "metadata field inventory");
  size_t position = 24 + 64;
  for (size_t j = 0; j < layout.metadata_offsets.size(); ++j) {
    layout.metadata_offsets[j] = position;
    position += 8;
  }
  layout.edges_offset = position;
  for (int j = 0; j <= info.cells; ++j) {
    layout.edge_offsets.push_back(position);
    position += 8;
  }
  const uint64_t knot_count = load_u64(bytes, layout.metadata("knots"));
  check(knot_count <= 100000, "packed knot count inventory");
  for (uint64_t k = 0; k < knot_count; ++k) {
    KnotWire knot;
    knot.temperature_offset = position;
    position += 8;
    for (size_t j = 0; j < knot.coefficient_offsets.size(); ++j) {
      knot.coefficient_offsets[j] = position;
      position += 8;
    }
    knot.count_offset = position;
    const uint64_t count = load_u64(bytes, position);
    position += 8;
    check(count <= uint64_t(7 * info.cells), "packed sparse count inventory");
    for (uint64_t j = 0; j < count; ++j) {
      knot.entries.push_back({position, position + 8});
      position += 16;
    }
    layout.knots.push_back(std::move(knot));
  }
  layout.checksum_offset = position;
  check(layout.checksum_offset + 8 == bytes.size(),
        "canonical wire layout reaches footer");
  return layout;
}

std::vector<unsigned char> pack_table(const fusion_beam_birth_table_v1 *table,
                                      size_t *required_out = nullptr) {
  size_t required = 0;
  check(fusion_c_beam_birth_table_pack_size(table, &required) ==
            PB11_STATUS_OK,
        "pack size");
  std::vector<unsigned char> bytes(required, 0);
  size_t written = 0;
  check(fusion_c_beam_birth_table_pack(table, bytes.data(), bytes.size(),
                                       &written) == PB11_STATUS_OK,
        "pack");
  check(written == required, "pack written size");
  if (required_out) *required_out = required;
  return bytes;
}

bool same_info(const fusion_beam_birth_table_info_v1 &a,
               const fusion_beam_birth_table_info_v1 &b) {
  if (a.projectile_energy_J != b.projectile_energy_J ||
      a.lower_kT_J != b.lower_kT_J || a.upper_kT_J != b.upper_kT_J ||
      a.max_validated_rate_error != b.max_validated_rate_error ||
      a.max_validated_debit_error != b.max_validated_debit_error ||
      a.max_validated_number_L1 != b.max_validated_number_L1 ||
      a.max_validated_energy_L1 != b.max_validated_energy_L1 ||
      a.max_sampled_direct_rate_discrepancy !=
          b.max_sampled_direct_rate_discrepancy ||
      a.max_sampled_direct_debit_discrepancy !=
          b.max_sampled_direct_debit_discrepancy ||
      a.channel != b.channel || a.projectile_slot != b.projectile_slot ||
      a.cells != b.cells || a.knots != b.knots ||
      a.direct_evaluations != b.direct_evaluations ||
      a.spectral_entries_evaluated != b.spectral_entries_evaluated ||
      a.stored_spectral_entries != b.stored_spectral_entries)
    return false;
  const auto &x = a.source;
  const auto &y = b.source;
  if (x.relative_max_J != y.relative_max_J ||
      x.angular_max_exponent != y.angular_max_exponent ||
      x.ground_state_q_J != y.ground_state_q_J || x.cutoff_J != y.cutoff_J ||
      x.l1_fraction != y.l1_fraction || x.relative_phase != y.relative_phase ||
      x.narrow_peak_fraction != y.narrow_peak_fraction ||
      x.continuum_peak_scale != y.continuum_peak_scale ||
      x.continuation != y.continuation || x.pb_low != y.pb_low ||
      x.remainder_policy != y.remainder_policy ||
      x.broad_mode != y.broad_mode || x.fsci_policy != y.fsci_policy ||
      x.relative_order != y.relative_order ||
      x.angular_order != y.angular_order || x.nq != y.nq ||
      x.ncos != y.ncos)
    return false;
  const auto &p = a.control;
  const auto &q = b.control;
  return p.max_rate_error == q.max_rate_error &&
         p.max_debit_error == q.max_debit_error &&
         p.max_number_L1 == q.max_number_L1 &&
         p.max_energy_L1 == q.max_energy_L1 &&
         p.max_direct_rate_discrepancy == q.max_direct_rate_discrepancy &&
         p.max_direct_debit_discrepancy == q.max_direct_debit_discrepancy &&
         p.max_knots == q.max_knots && p.max_evaluations == q.max_evaluations &&
         p.max_depth == q.max_depth;
}

struct Evaluation {
  std::vector<double> grid;
  fusion_birth_coefficients_v1 coefficients{};
};

Evaluation evaluate(const fusion_beam_birth_table_v1 *table, int cells,
                    double temperature) {
  Evaluation result;
  result.grid.assign(static_cast<size_t>(7 * cells), 7.0);
  std::fill(reinterpret_cast<unsigned char *>(&result.coefficients),
            reinterpret_cast<unsigned char *>(&result.coefficients) +
                sizeof(result.coefficients),
            0xA5);
  check(fusion_c_beam_birth_table_evaluate(
            table, temperature, cells, result.grid.data(),
            &result.coefficients) == PB11_STATUS_OK,
        "table evaluation");
  return result;
}

void exact_evaluation(const Evaluation &a, const Evaluation &b,
                      const std::string &where) {
  check(a.grid == b.grid, where + " grid parity");
  check(std::memcmp(&a.coefficients, &b.coefficients,
                    sizeof(a.coefficients)) == 0,
        where + " coefficient parity");
}

fusion_birth_coefficients_v1 direct_coefficients(
    const fusion_beam_birth_v1 &b) {
  fusion_birth_coefficients_v1 result{};
  const auto &s = b.spectrum;
  result.reactivity_m3_s = s.reactivity_m3_s;
  for (int j = 0; j < 2; ++j)
    result.reactant_energy_moment_J_m3_s[j] =
        s.reactant_energy_moment_J_m3_s[j];
  for (int j = 0; j < 7; ++j) {
    result.below_number_m3_s[j] = s.below_number_m3_s[j];
    result.below_energy_J_m3_s[j] = s.below_energy_J_m3_s[j];
    result.above_number_m3_s[j] = s.above_number_m3_s[j];
    result.above_energy_J_m3_s[j] = s.above_energy_J_m3_s[j];
  }
  return result;
}

void test_identity_and_roundtrip(
    const Fixture &fixture, const std::vector<unsigned char> &bytes,
    const WireLayout &layout) {
  const char *identity = fusion_c_beam_birth_table_kernel_identity();
  check(identity && std::strlen(identity) == 64, "kernel identity length");
  for (size_t j = 0; j < 64; ++j)
    check(std::isxdigit(static_cast<unsigned char>(identity[j])) != 0,
          "kernel identity hexadecimal");
  check(std::memcmp(bytes.data() + layout.identity_offset, identity, 64) == 0,
        "wire identity");

  fusion_beam_birth_table_v1 *raw = nullptr;
  check(fusion_c_beam_birth_table_unpack(bytes.data(), bytes.size(), &raw) ==
            PB11_STATUS_OK &&
            raw,
        "unpack roundtrip");
  owner<fusion_beam_birth_table_v1> restored(raw,
                                              fusion_c_beam_birth_table_destroy);
  fusion_beam_birth_table_info_v1 before{}, after{};
  check(fusion_c_beam_birth_table_info(fixture.table.get(), &before) ==
            PB11_STATUS_OK,
        "original info");
  check(fusion_c_beam_birth_table_info(restored.get(), &after) ==
            PB11_STATUS_OK,
        "restored info");
  check(same_info(before, after), "roundtrip metadata parity");

  const std::vector<unsigned char> repacked = pack_table(restored.get());
  check(repacked == bytes, "byte-identical repack");

  for (const double fraction : {0.0, .17, .5, .83, 1.0}) {
    const double temperature =
        fraction == 0.0
            ? fixture.lower_kT_J
            : fraction == 1.0
                  ? fixture.upper_kT_J
                  : std::exp((1 - fraction) * std::log(fixture.lower_kT_J) +
                              fraction * std::log(fixture.upper_kT_J));
    exact_evaluation(evaluate(fixture.table.get(), before.cells, temperature),
                     evaluate(restored.get(), after.cells, temperature),
                     "independent unpack evaluation");
  }

  // The two endpoint knots are also checked against a fresh direct source
  // evaluation, which does not inspect the persistence representation.
  for (const double temperature : {fixture.lower_kT_J, fixture.upper_kT_J}) {
    const int cells = before.cells;
    std::vector<double> direct_grid(static_cast<size_t>(7 * cells));
    fusion_beam_birth_v1 direct{};
    check(fusion_c_beam_birth_grid(
              fixture.channel, fixture.slot, fixture.projectile_energy_J,
              temperature, &fixture.options, cells, fixture.edges.data(),
              direct_grid.data(), &direct) == PB11_STATUS_OK,
          "independent direct endpoint");
    const Evaluation table_eval =
        evaluate(restored.get(), cells, temperature);
    check(table_eval.grid == direct_grid, "direct endpoint grid parity");
    const fusion_birth_coefficients_v1 direct_coeff =
        direct_coefficients(direct);
    check(std::memcmp(&table_eval.coefficients, &direct_coeff,
                      sizeof(table_eval.coefficients)) == 0,
          "direct endpoint coefficient parity");
  }
}

void test_null_and_short_buffer(
    const Fixture &fixture, const std::vector<unsigned char> &bytes) {
  size_t required = 123;
  check(fusion_c_beam_birth_table_pack_size(fixture.table.get(), nullptr) ==
            PB11_STATUS_NULL_OUTPUT,
        "pack size null output");
  check(fusion_c_beam_birth_table_pack_size(nullptr, &required) !=
            PB11_STATUS_OK &&
            required == 0,
        "pack size null table clears count");

  size_t written = 123;
  check(fusion_c_beam_birth_table_pack(fixture.table.get(), nullptr, 0,
                                       &written) != PB11_STATUS_OK &&
            written == 0,
        "pack null buffer clears count");
  written = 123;
  check(fusion_c_beam_birth_table_pack(
            nullptr, const_cast<unsigned char *>(bytes.data()), bytes.size(),
            &written) != PB11_STATUS_OK &&
            written == 0,
        "pack null table clears count");
  check(fusion_c_beam_birth_table_pack(
            fixture.table.get(), const_cast<unsigned char *>(bytes.data()),
            bytes.size(), nullptr) == PB11_STATUS_NULL_OUTPUT,
        "pack null count output");

  std::vector<unsigned char> short_buffer(bytes.size() + 16, 0xCD);
  const std::vector<unsigned char> before = short_buffer;
  written = 123;
  check(fusion_c_beam_birth_table_pack(
            fixture.table.get(), short_buffer.data(), bytes.size() - 1,
            &written) != PB11_STATUS_OK &&
            written == 0,
        "short pack rejects and clears count");
  check(short_buffer == before, "short pack leaves buffer unchanged");

  fusion_beam_birth_table_v1 *out =
      reinterpret_cast<fusion_beam_birth_table_v1 *>(uintptr_t(1));
  check(fusion_c_beam_birth_table_unpack(bytes.data(), bytes.size(), nullptr) ==
            PB11_STATUS_NULL_OUTPUT,
        "unpack null output pointer");
  check(fusion_c_beam_birth_table_unpack(nullptr, bytes.size(), &out) !=
            PB11_STATUS_OK &&
            out == nullptr,
        "unpack null buffer clears table");
  out = reinterpret_cast<fusion_beam_birth_table_v1 *>(uintptr_t(1));
  check(fusion_c_beam_birth_table_unpack(bytes.data(), bytes.size() - 1, &out) !=
            PB11_STATUS_OK &&
            out == nullptr,
        "truncated unpack clears table");

  check(fusion_c_beam_birth_table_info(fixture.table.get(), nullptr) ==
            PB11_STATUS_NULL_OUTPUT,
        "info null output");
}

void expect_bad(const std::vector<unsigned char> &bytes,
                const std::string &description) {
  fusion_beam_birth_table_v1 *out =
      reinterpret_cast<fusion_beam_birth_table_v1 *>(uintptr_t(1));
  const int status =
      fusion_c_beam_birth_table_unpack(bytes.data(), bytes.size(), &out);
  check(status != PB11_STATUS_OK && out == nullptr,
        description + " rejects without partial table");
}

template <class Mutator>
void expect_rechecksum_bad(const std::vector<unsigned char> &original,
                           Mutator mutate, const std::string &description) {
  std::vector<unsigned char> candidate = original;
  mutate(candidate);
  rechecksum(candidate);
  expect_bad(candidate, description);
}

void test_wire_rejections(
    const Fixture &fixture, const std::vector<unsigned char> &bytes,
    const WireLayout &layout) {
  std::vector<unsigned char> wrong_magic = bytes;
  store_u64(wrong_magic, 0, 0xDEADBEEF);
  rechecksum(wrong_magic);
  expect_bad(wrong_magic, "wrong magic");

  std::vector<unsigned char> wrong_version = bytes;
  store_u64(wrong_version, 8, 2);
  rechecksum(wrong_version);
  expect_bad(wrong_version, "wrong version");

  std::vector<unsigned char> wrong_identity = bytes;
  wrong_identity[layout.identity_offset] ^= 1;
  rechecksum(wrong_identity);
  expect_bad(wrong_identity, "wrong identity");

  std::vector<unsigned char> bad_checksum = bytes;
  bad_checksum[layout.metadata("projectile_energy_J")] ^= 1;
  expect_bad(bad_checksum, "bad checksum");

  std::vector<unsigned char> truncated = bytes;
  truncated.resize(truncated.size() - 1);
  expect_bad(truncated, "truncated payload");

  // Insert an extra byte before a freshly recomputed footer. The checksum is
  // valid, so this specifically exercises the stored-size/canonical-length
  // check rather than merely the checksum guard.
  std::vector<unsigned char> extra;
  extra.reserve(bytes.size() + 1);
  extra.insert(extra.end(), bytes.begin(), bytes.end() - 8);
  extra.push_back(0xA5);
  extra.insert(extra.end(), bytes.end() - 8, bytes.end());
  rechecksum(extra);
  expect_bad(extra, "extra payload byte");

  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, layout.metadata("projectile_energy_J"),
                     std::numeric_limits<double>::quiet_NaN());
      },
      "NaN fixed metadata");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("channel"), 99);
      },
      "invalid channel metadata");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("projectile_slot"), 2);
      },
      "invalid projectile slot metadata");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) { store_u64(candidate, layout.metadata("cells"), 0); },
      "invalid cell count metadata");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) { store_u64(candidate, layout.metadata("knots"), 1); },
      "invalid knot count metadata");

  // Options and controls use the same metadata stream but have distinct
  // semantic domains, so exercise each family independently.
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, layout.metadata("source.angular_max_exponent"),
                     7);
      },
      "invalid source angular exponent");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("source.continuation"), 0);
      },
      "invalid source continuation");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("source.broad_mode"), 2);
      },
      "invalid source broad mode");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, layout.metadata("source.l1_fraction"), 2);
      },
      "invalid source fraction");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, layout.metadata("control.max_rate_error"), 2);
      },
      "invalid control rate gate");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("control.max_knots"), 1);
      },
      "invalid control knot cap");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("control.max_evaluations"), 4);
      },
      "invalid control evaluation cap");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("control.max_depth"), 25);
      },
      "invalid control depth cap");

  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("spectral_entries_evaluated"),
                  entry_cap + 1);
      },
      "invalid sampled entry count");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, layout.metadata("stored_spectral_entries"), 0);
      },
      "stored entry count mismatch");

  const auto first_knot = layout.knots.front();
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_u64(candidate, first_knot.count_offset,
                  std::numeric_limits<uint64_t>::max());
      },
      "oversized sparse count");

  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, layout.edge_offsets.front(), -1);
      },
      "negative grid edge");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, layout.edge_offsets[2],
                     load_double(candidate, layout.edge_offsets[1]));
      },
      "duplicate grid edge");

  const KnotWire *grid_knot = nullptr;
  for (const auto &knot : layout.knots)
    if (!knot.entries.empty()) {
      grid_knot = &knot;
      break;
    }
  check(grid_knot, "ordinary fixture contains sparse entries");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, grid_knot->entries.front().value_offset, -1);
      },
      "negative sparse value");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, grid_knot->entries.front().value_offset,
                     std::numeric_limits<double>::quiet_NaN());
      },
      "NaN sparse value");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, grid_knot->temperature_offset,
                     fixture.lower_kT_J * .5);
      },
      "out of range knot temperature");

  const KnotWire *duplicate_knot = nullptr;
  for (const auto &knot : layout.knots)
    if (knot.entries.size() >= 2) {
      duplicate_knot = &knot;
      break;
    }
  check(duplicate_knot, "ordinary fixture contains two sparse entries");
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        const uint64_t first_index =
            load_u64(candidate, duplicate_knot->entries[0].index_offset);
        store_u64(candidate, duplicate_knot->entries[1].index_offset,
                  first_index);
      },
      "duplicate sparse index");

  // Changing one reacting-energy debit while retaining a valid checksum is a
  // semantic conservation violation, not a malformed floating-point value.
  expect_rechecksum_bad(
      bytes,
      [&](auto &candidate) {
        store_double(candidate, grid_knot->coefficient_offsets[1],
                     load_double(candidate, grid_knot->coefficient_offsets[1]) +
                         1.0);
      },
      "nonconservative coefficient");
}

void test_subnormal_roundtrip() {
  const Fixture fixture = make_subnormal_fixture();
  size_t required = 0;
  const std::vector<unsigned char> bytes =
      pack_table(fixture.table.get(), &required);
  fusion_beam_birth_table_info_v1 info{};
  check(fusion_c_beam_birth_table_info(fixture.table.get(), &info) ==
            PB11_STATUS_OK,
        "subnormal info");
  const WireLayout layout = locate_wire(bytes, info);
  size_t subnormal_values = 0;
  for (const auto &knot : layout.knots)
    for (const auto &entry : knot.entries)
      if (std::fpclassify(load_double(bytes, entry.value_offset)) == FP_SUBNORMAL)
        ++subnormal_values;
  check(subnormal_values > 0, "subnormal fixture has stored nonzero entries");

  fusion_beam_birth_table_v1 *raw = nullptr;
  check(fusion_c_beam_birth_table_unpack(bytes.data(), bytes.size(), &raw) ==
            PB11_STATUS_OK &&
            raw,
        "subnormal unpack");
  owner<fusion_beam_birth_table_v1> restored(raw,
                                              fusion_c_beam_birth_table_destroy);
  const std::vector<unsigned char> repacked = pack_table(restored.get());
  check(repacked == bytes, "subnormal byte-identical repack");
  const Evaluation original =
      evaluate(fixture.table.get(), info.cells, fixture.lower_kT_J);
  const Evaluation roundtrip =
      evaluate(restored.get(), info.cells, fixture.lower_kT_J);
  exact_evaluation(original, roundtrip, "subnormal endpoint");
  size_t output_subnormal = 0;
  for (double value : roundtrip.grid)
    if (std::fpclassify(value) == FP_SUBNORMAL) ++output_subnormal;
  check(output_subnormal > 0, "subnormal output retained");
  (void)required;
}

} // namespace

int main() {
  try {
    const Fixture fixture = make_ordinary_fixture();
    size_t required = 0;
    const std::vector<unsigned char> bytes =
        pack_table(fixture.table.get(), &required);
    check(required == bytes.size() && required > 96, "packed size");
    fusion_beam_birth_table_info_v1 info{};
    check(fusion_c_beam_birth_table_info(fixture.table.get(), &info) ==
              PB11_STATUS_OK,
          "ordinary info for layout");
    const WireLayout layout = locate_wire(bytes, info);
    test_identity_and_roundtrip(fixture, bytes, layout);
    test_null_and_short_buffer(fixture, bytes);
    test_wire_rejections(fixture, bytes, layout);
    test_subnormal_roundtrip();
    std::cout << "PASS: beam persistence identity, canonical roundtrip, parity, "
                 "atomic failures, malformed wires, semantic guards, subnormals\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
