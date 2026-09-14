#include "fusion_source_state.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kElectronVoltJ = 1.602176634e-19;
constexpr double kMeVJ = 1.0e6 * kElectronVoltJ;
constexpr int kSpecies = FUSION_SPECIES_COUNT;
constexpr int kAlpha = FUSION_HELIUM4;
constexpr std::size_t kLedgerWords = 5 + 12 * kSpecies + 42 + 2;
constexpr std::size_t kSerializedPrefixWords = 8 + 12 + kLedgerWords;

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

bool same_ledger(const fusion_source_ledger_v1& actual,
                const fusion_source_ledger_v1& expected) {
    for (int i = 0; i < FUSION_CHANNEL_COUNT; ++i) {
        if (actual.events_m3[i] != expected.events_m3[i]) return false;
    }
    for (int i = 0; i < kSpecies; ++i) {
        if (actual.nuclear_born_number_m3[i] !=
                expected.nuclear_born_number_m3[i] ||
            actual.nuclear_born_energy_J_m3[i] !=
                expected.nuclear_born_energy_J_m3[i] ||
            actual.external_born_number_m3[i] !=
                expected.external_born_number_m3[i] ||
            actual.external_born_energy_J_m3[i] !=
                expected.external_born_energy_J_m3[i] ||
            actual.thermal_consumed_number_m3[i] !=
                expected.thermal_consumed_number_m3[i] ||
            actual.thermal_consumed_energy_J_m3[i] !=
                expected.thermal_consumed_energy_J_m3[i] ||
            actual.fast_consumed_number_m3[i] !=
                expected.fast_consumed_number_m3[i] ||
            actual.fast_consumed_energy_J_m3[i] !=
                expected.fast_consumed_energy_J_m3[i] ||
            actual.escaped_number_m3[i] != expected.escaped_number_m3[i] ||
            actual.escaped_energy_J_m3[i] != expected.escaped_energy_J_m3[i] ||
            actual.handed_off_number_m3[i] !=
                expected.handed_off_number_m3[i] ||
            actual.handed_off_energy_J_m3[i] !=
                expected.handed_off_energy_J_m3[i]) {
            return false;
        }
    }
    for (int i = 0; i < kSpecies * 7; ++i) {
        if (actual.heat_to_bath_J_m3[i] != expected.heat_to_bath_J_m3[i]) {
            return false;
        }
    }
    return actual.neutron_number_m3 == expected.neutron_number_m3 &&
           actual.neutron_energy_J_m3 == expected.neutron_energy_J_m3;
}

void check_zero_ledger(const fusion_source_ledger_v1& ledger,
                       const std::string& message) {
    const fusion_source_ledger_v1 zero{};
    check(same_ledger(ledger, zero), message);
}

bool same_heat(const std::array<double, kSpecies>& actual,
               const std::array<double, kSpecies>& expected) {
    for (int i = 0; i < kSpecies; ++i) {
        if (actual[static_cast<std::size_t>(i)] !=
            expected[static_cast<std::size_t>(i)]) {
            return false;
        }
    }
    return true;
}

struct StateCase {
    int cells;
    std::vector<double> edges;
    std::vector<double> initial_s;
    std::vector<double> initial_t;
    std::uint64_t tag;
    fusion_source_state_v1* state = nullptr;

    explicit StateCase(int cell_count = 2,
                       std::uint64_t model_tag = UINT64_C(0x1122334455667788))
        : cells(cell_count),
          edges(static_cast<std::size_t>(cell_count) + 1),
          initial_s(static_cast<std::size_t>(kSpecies) * cell_count, 0.0),
          initial_t(static_cast<std::size_t>(kSpecies) * cell_count, 0.0),
          tag(model_tag) {
        check(cells >= 2, "test case needs at least two energy cells");
        for (int i = 0; i <= cells; ++i) {
            edges[static_cast<std::size_t>(i)] =
                static_cast<double>(i) * 2.0 * kMeVJ;
        }
        initial_s[index(kAlpha, cells - 1)] = 1.0e12;
    }

    ~StateCase() { fusion_c_source_state_destroy(state); }

    StateCase(const StateCase&) = delete;
    StateCase& operator=(const StateCase&) = delete;

    std::size_t index(int species, int cell) const {
        return static_cast<std::size_t>(species) *
                   static_cast<std::size_t>(cells) +
               static_cast<std::size_t>(cell);
    }

    void create(double initial_time = 3.25) {
        expect_ok(fusion_c_source_state_create(
                      cells, edges.data(), initial_s.data(), initial_t.data(),
                      initial_time, tag, &state),
                  "create inert source-state fixture");
        check(state != nullptr, "source-state fixture has a context");
    }

    std::vector<double> alpha_in_cell(int cell) const {
        check(cell >= 0 && cell < cells - 1,
              "fixture target cell is below initial fast cell");
        std::vector<double> candidate = initial_s;
        candidate[index(kAlpha, cells - 1)] = 0.0;
        candidate[index(kAlpha, cell)] = 1.0e12;
        return candidate;
    }

    std::vector<double> all_species_in_cell(int cell) const {
        check(cell >= 0 && cell < cells - 1,
              "fixture target cell is below initial fast cell");
        std::vector<double> candidate = initial_s;
        for (int species = 0; species < kSpecies; ++species) {
            candidate[index(species, cells - 1)] = 0.0;
            candidate[index(species, cell)] =
                initial_s[index(species, cells - 1)];
        }
        return candidate;
    }

    double heat_for_move(int species, int from_cell, int to_cell) const {
        const long double population =
            static_cast<long double>(initial_s[index(species, from_cell)]) +
            static_cast<long double>(initial_t[index(species, from_cell)]);
        const long double from_energy =
            (static_cast<long double>(edges[static_cast<std::size_t>(from_cell)]) +
             static_cast<long double>(edges[static_cast<std::size_t>(from_cell + 1)])) /
            2.0L;
        const long double to_energy =
            (static_cast<long double>(edges[static_cast<std::size_t>(to_cell)]) +
             static_cast<long double>(edges[static_cast<std::size_t>(to_cell + 1)])) /
            2.0L;
        return static_cast<double>(population * (from_energy - to_energy));
    }

    double heat_for_move(int from_cell, int to_cell) const {
        return heat_for_move(kAlpha, from_cell, to_cell);
    }
};

fusion_source_ledger_v1 zero_step() { return fusion_source_ledger_v1{}; }

void snapshot_inert(const fusion_source_state_v1* state, int cells,
                    std::vector<double>& accepted_s,
                    std::vector<double>& accepted_t,
                    fusion_source_ledger_v1& cumulative,
                    std::array<double, kSpecies>& inert_heat,
                    double& accepted_time, std::uint64_t& epoch) {
    accepted_s.assign(static_cast<std::size_t>(kSpecies) * cells, 0.0);
    accepted_t.assign(static_cast<std::size_t>(kSpecies) * cells, 0.0);
    expect_ok(fusion_c_source_state_snapshot_inert(
                  state, accepted_s.data(), accepted_t.data(), &cumulative,
                  inert_heat.data(), &accepted_time, &epoch),
              "snapshot accepted inert source-state");
}

std::vector<unsigned char> pack_state(const fusion_source_state_v1* state) {
    std::size_t required = 0;
    expect_ok(fusion_c_source_state_pack_size(state, &required),
              "query source-state restart size");
    check(required != 0, "restart packet has a nonzero size");
    std::vector<unsigned char> bytes(required, 0);
    std::size_t written = 0;
    expect_ok(fusion_c_source_state_pack(state, bytes.data(), bytes.size(),
                                          &written),
              "pack accepted source-state");
    check(written == required, "restart packet writes its required size");
    return bytes;
}

std::uint64_t read_u64_little_endian(const unsigned char* data) {
    std::uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        value |= static_cast<std::uint64_t>(data[i]) << (8 * i);
    }
    return value;
}

void write_u64_little_endian(unsigned char* data, std::uint64_t value) {
    for (int i = 0; i < 8; ++i) {
        data[i] = static_cast<unsigned char>(value & UINT64_C(0xff));
        value >>= 8;
    }
}

double read_double_little_endian(const unsigned char* data) {
    const std::uint64_t bits = read_u64_little_endian(data);
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

void write_double_little_endian(unsigned char* data, double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    write_u64_little_endian(data, bits);
}

std::uint64_t restart_checksum(const unsigned char* data, std::size_t count) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (std::size_t i = 0; i < count; ++i) {
        hash ^= static_cast<std::uint64_t>(data[i]);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

void refresh_restart_checksum(std::vector<unsigned char>& bytes) {
    check(bytes.size() >= sizeof(std::uint64_t),
          "restart packet has checksum storage");
    write_u64_little_endian(
        bytes.data() + bytes.size() - sizeof(std::uint64_t),
        restart_checksum(bytes.data(), bytes.size() - sizeof(std::uint64_t)));
}

void test_signed_inert_heat_and_snapshot_contract() {
    StateCase fixture;
    for (int species = 0; species < kSpecies; ++species) {
        fixture.initial_s[fixture.index(species, fixture.cells - 1)] =
            1.0e12 * static_cast<double>(species + 1);
    }
    fixture.create();
    const std::vector<double> lower = fixture.all_species_in_cell(0);
    std::array<double, kSpecies> expected_heat{};
    for (int species = 0; species < kSpecies; ++species) {
        expected_heat[static_cast<std::size_t>(species)] =
            fixture.heat_for_move(species, fixture.cells - 1, 0);
        check(expected_heat[static_cast<std::size_t>(species)] > 0.0 &&
                  std::isfinite(expected_heat[static_cast<std::size_t>(species)]),
              "fixture has a finite positive fast-energy drop for every species");
    }

    const auto step = zero_step();
    std::uint64_t ticket = 0;
    expect_ok(fusion_c_source_state_begin(fixture.state, 0.25, &ticket),
              "begin energy-drop trial");

    /* The original stage API supplies no non-network-bath account.  The
     * candidate therefore cannot pass its energy identity. */
    expect_rejected(fusion_c_source_state_stage(
                        fixture.state, ticket, lower.data(),
                        fixture.initial_t.data(), &step),
                    "legacy stage rejects an unaccounted fast-energy drop");
    expect_rejected(fusion_c_source_state_commit(fixture.state, ticket),
                    "failed legacy stage cannot commit");

    expect_ok(fusion_c_source_state_stage_inert(
                  fixture.state, ticket, lower.data(), fixture.initial_t.data(),
                  &step, expected_heat.data()),
              "extended stage accepts the complete inert heat account");
    expect_ok(fusion_c_source_state_commit(fixture.state, ticket),
              "commit complete inert heat account");

    std::vector<double> accepted_s, accepted_t;
    fusion_source_ledger_v1 cumulative{};
    std::array<double, kSpecies> actual_heat{};
    double accepted_time = 0.0;
    std::uint64_t epoch = 0;
    snapshot_inert(fixture.state, fixture.cells, accepted_s, accepted_t,
                   cumulative, actual_heat, accepted_time, epoch);
    check(accepted_s == lower && accepted_t == fixture.initial_t,
          "snapshot_inert stores the accepted lower-energy fast inventory");
    check(same_heat(actual_heat, expected_heat),
          "snapshot_inert preserves all six cumulative inert heat species");
    check_zero_ledger(cumulative,
                      "inert-only step leaves the base cumulative ledger zero");
    check(accepted_time == 3.25 + 0.25 && epoch == 1,
          "inert-only commit advances accepted time and epoch once");

    /* The old snapshot cannot silently discard the promoted account.  Its
     * kinetic outputs stay untouched, while scalar and ledger outputs clear. */
    std::vector<double> old_s(accepted_s.size(), -17.0);
    std::vector<double> old_t(accepted_t.size(), -19.0);
    const std::vector<double> old_s_before = old_s;
    const std::vector<double> old_t_before = old_t;
    fusion_source_ledger_v1 old_ledger{};
    old_ledger.events_m3[0] = 4.0;
    old_ledger.heat_to_bath_J_m3[0] = 5.0;
    old_ledger.neutron_energy_J_m3 = 6.0;
    double old_time = 7.0;
    std::uint64_t old_epoch = 8;
    expect_rejected(fusion_c_source_state_snapshot(
                        fixture.state, old_s.data(), old_t.data(), &old_ledger,
                        &old_time, &old_epoch),
                    "old snapshot rejects an extended context");
    check(old_s == old_s_before && old_t == old_t_before,
          "failed old snapshot leaves kinetic output arrays untouched");
    check_zero_ledger(old_ledger,
                      "failed old snapshot clears its ledger output");
    check(old_time == 0.0 && old_epoch == 0,
          "failed old snapshot clears scalar outputs");

    /* Once promoted, the original API remains usable for a zero-extra-heat
     * step and carries the already accepted inert total forward. */
    std::uint64_t next_ticket = 0;
    expect_ok(fusion_c_source_state_begin(fixture.state, 0.1, &next_ticket),
              "begin zero-extra step on extended context");
    expect_ok(fusion_c_source_state_stage(
                  fixture.state, next_ticket, lower.data(),
                  fixture.initial_t.data(), &step),
              "legacy stage supplies zero extra heat on extended context");
    expect_ok(fusion_c_source_state_commit(fixture.state, next_ticket),
              "commit zero-extra step on extended context");
    snapshot_inert(fixture.state, fixture.cells, accepted_s, accepted_t,
                   cumulative, actual_heat, accepted_time, epoch);
    check(same_heat(actual_heat, expected_heat) && epoch == 2,
          "legacy stage preserves the cumulative inert heat after promotion");
}

void test_discard_and_rejected_extended_stage_do_not_promote() {
    {
        StateCase fixture;
        fixture.create();
        const std::vector<double> lower = fixture.alpha_in_cell(0);
        const double heat = fixture.heat_for_move(1, 0);
        const auto step = zero_step();
        std::array<double, kSpecies> discard_heat{};
        discard_heat[static_cast<std::size_t>(kAlpha)] = heat;
        std::uint64_t ticket = 0;
        expect_ok(fusion_c_source_state_begin(fixture.state, 0.2, &ticket),
                  "begin trial that will be discarded");
        expect_ok(fusion_c_source_state_stage_inert(
                      fixture.state, ticket, lower.data(), fixture.initial_t.data(),
                      &step, discard_heat.data()),
                  "stage valid extended trial before discard");
        expect_ok(fusion_c_source_state_discard(fixture.state, ticket),
                  "discard valid extended trial");
        expect_rejected(fusion_c_source_state_commit(fixture.state, ticket),
                        "discarded extended trial cannot commit");

        std::vector<double> accepted_s(
            static_cast<std::size_t>(kSpecies) * fixture.cells, 0.0);
        std::vector<double> accepted_t(
            static_cast<std::size_t>(kSpecies) * fixture.cells, 0.0);
        fusion_source_ledger_v1 cumulative{};
        std::array<double, kSpecies> inert{};
        double accepted_time = 0.0;
        std::uint64_t epoch = 0;
        expect_ok(fusion_c_source_state_snapshot(
                      fixture.state, accepted_s.data(), accepted_t.data(),
                      &cumulative, &accepted_time, &epoch),
                  "discarded trial leaves a legacy snapshot");
        snapshot_inert(fixture.state, fixture.cells, accepted_s, accepted_t,
                       cumulative, inert, accepted_time, epoch);
        check(accepted_s == fixture.initial_s && accepted_t == fixture.initial_t,
              "discard leaves the accepted kinetic inventory unchanged");
        check(same_heat(inert, std::array<double, kSpecies>{}),
              "discard prevents inert promotion");
        check_zero_ledger(cumulative,
                          "discard leaves the cumulative ledger unchanged");
    }

    {
        StateCase fixture;
        fixture.create();
        const std::vector<double> lower = fixture.alpha_in_cell(0);
        const auto step = zero_step();
        std::uint64_t ticket = 0;
        expect_ok(fusion_c_source_state_begin(fixture.state, 0.2, &ticket),
                  "begin trial with invalid inert heat");
        std::array<double, kSpecies> invalid_heat{};
        invalid_heat[static_cast<std::size_t>(kAlpha)] =
            std::numeric_limits<double>::quiet_NaN();
        expect_rejected(fusion_c_source_state_stage_inert(
                            fixture.state, ticket, lower.data(),
                            fixture.initial_t.data(), &step,
                            invalid_heat.data()),
                        "nonfinite inert heat is rejected");
        expect_rejected(fusion_c_source_state_commit(fixture.state, ticket),
                        "failed nonfinite inert stage cannot commit");
        expect_ok(fusion_c_source_state_discard(fixture.state, ticket),
                  "discard invalid inert trial");

        std::vector<double> accepted_s, accepted_t;
        fusion_source_ledger_v1 cumulative{};
        std::array<double, kSpecies> inert{};
        double accepted_time = 0.0;
        std::uint64_t epoch = 0;
        snapshot_inert(fixture.state, fixture.cells, accepted_s, accepted_t,
                       cumulative, inert, accepted_time, epoch);
        check(accepted_s == fixture.initial_s && accepted_t == fixture.initial_t,
              "rejected inert stage leaves accepted inventory unchanged");
        check(same_heat(inert, std::array<double, kSpecies>{}),
              "rejected inert stage cannot promote the context");
        check_zero_ledger(cumulative,
                          "rejected inert stage leaves the ledger unchanged");
    }
}

void test_repeated_extended_stage_replaces_without_accumulation() {
    StateCase fixture(3);
    fixture.create();
    const std::vector<double> middle = fixture.alpha_in_cell(1);
    const std::vector<double> lower = fixture.alpha_in_cell(0);
    const double first_heat = fixture.heat_for_move(2, 1);
    const double second_heat = fixture.heat_for_move(2, 0);
    check(second_heat > first_heat && first_heat > 0.0,
          "three-cell fixture gives two distinct energy drops");
    const auto step = zero_step();
    std::uint64_t ticket = 0;
    expect_ok(fusion_c_source_state_begin(fixture.state, 0.125, &ticket),
              "begin repeated extended trial");

    std::array<double, kSpecies> first_inert{};
    first_inert[static_cast<std::size_t>(kAlpha)] = first_heat;
    expect_ok(fusion_c_source_state_stage_inert(
                  fixture.state, ticket, middle.data(), fixture.initial_t.data(),
                  &step, first_inert.data()),
              "stage first extended nonlinear evaluation");
    std::array<double, kSpecies> second_inert{};
    second_inert[static_cast<std::size_t>(kAlpha)] = second_heat;
    expect_ok(fusion_c_source_state_stage_inert(
                  fixture.state, ticket, lower.data(), fixture.initial_t.data(),
                  &step, second_inert.data()),
              "stage replacement extended nonlinear evaluation");
    expect_ok(fusion_c_source_state_commit(fixture.state, ticket),
              "commit replacement extended evaluation");

    std::vector<double> accepted_s, accepted_t;
    fusion_source_ledger_v1 cumulative{};
    std::array<double, kSpecies> actual_inert{};
    double accepted_time = 0.0;
    std::uint64_t epoch = 0;
    snapshot_inert(fixture.state, fixture.cells, accepted_s, accepted_t,
                   cumulative, actual_inert, accepted_time, epoch);
    check(accepted_s == lower && accepted_t == fixture.initial_t,
          "replacement stage commits only the second kinetic candidate");
    check(actual_inert[static_cast<std::size_t>(kAlpha)] == second_heat,
          "replacement stage stores only the second inert heat amount");
    check(actual_inert[static_cast<std::size_t>(kAlpha)] !=
              first_heat + second_heat,
          "replacement stage does not accumulate nonlinear evaluations");
    check_zero_ledger(cumulative,
                      "replacement inert trial keeps the expected base ledger");
    check(accepted_time == 3.25 + 0.125 && epoch == 1,
          "replacement extended trial commits one step");
}

void test_v1_and_v2_restart_compatibility() {
    StateCase legacy;
    legacy.create();
    const std::vector<unsigned char> v1 = pack_state(legacy.state);
    check(read_u64_little_endian(v1.data() + sizeof(std::uint64_t)) == 1,
          "unpromoted context writes restart version 1");

    std::vector<double> accepted_s, accepted_t;
    fusion_source_ledger_v1 cumulative{};
    std::array<double, kSpecies> inert{};
    double accepted_time = 0.0;
    std::uint64_t epoch = 0;
    snapshot_inert(legacy.state, legacy.cells, accepted_s, accepted_t,
                   cumulative, inert, accepted_time, epoch);
    check(same_heat(inert, std::array<double, kSpecies>{}),
          "snapshot_inert reports zero for a legacy context");
    check(pack_state(legacy.state) == v1,
          "snapshot_inert does not promote or rewrite legacy v1 bytes");

    fusion_source_state_v1* legacy_restored = nullptr;
    expect_ok(fusion_c_source_state_unpack(v1.data(), v1.size(), legacy.tag,
                                            &legacy_restored),
              "unpack legacy v1 bytes");
    check(legacy_restored != nullptr, "unpack v1 returns a context");
    snapshot_inert(legacy_restored, legacy.cells, accepted_s, accepted_t,
                   cumulative, inert, accepted_time, epoch);
    check(same_heat(inert, std::array<double, kSpecies>{}),
          "unpacked v1 context has zero inert heat");
    const std::vector<unsigned char> v1_round_trip =
        pack_state(legacy_restored);
    check(v1_round_trip == v1,
          "legacy v1 bytes survive unpack and repack byte-for-byte");
    fusion_c_source_state_destroy(legacy_restored);

    StateCase extended;
    extended.create();
    const std::vector<double> lower = extended.alpha_in_cell(0);
    const double heat = extended.heat_for_move(1, 0);
    const auto step = zero_step();
    std::array<double, kSpecies> first_inert{};
    first_inert[static_cast<std::size_t>(kAlpha)] = heat;
    std::uint64_t ticket = 0;
    expect_ok(fusion_c_source_state_begin(extended.state, 0.25, &ticket),
              "begin context that will be promoted");
    expect_ok(fusion_c_source_state_stage_inert(
                  extended.state, ticket, lower.data(),
                  extended.initial_t.data(), &step, first_inert.data()),
              "stage context that will be promoted");
    expect_ok(fusion_c_source_state_commit(extended.state, ticket),
              "commit context that will be promoted");

    const std::vector<unsigned char> v2 = pack_state(extended.state);
    check(read_u64_little_endian(v2.data() + sizeof(std::uint64_t)) == 2,
          "promoted context writes restart version 2");
    check(v2.size() == v1.size() + 6 * sizeof(double),
          "version 2 adds exactly six inert heat words");

    /* A checksum-valid change to inert heat must still fail the energy
     * identity during unpack. */
    std::vector<unsigned char> semantic = v2;
    const std::size_t inert_offset = kSerializedPrefixWords * sizeof(double);
    check(semantic.size() >= inert_offset +
                            static_cast<std::size_t>(kSpecies) * sizeof(double) +
                            sizeof(std::uint64_t),
          "version 2 packet exposes all inert heat words");
    const std::size_t alpha_offset =
        inert_offset + static_cast<std::size_t>(kAlpha) * sizeof(double);
    const double serialized_heat =
        read_double_little_endian(semantic.data() + alpha_offset);
    check(serialized_heat == heat,
          "version 2 packet places alpha inert heat after the base ledger");
    write_double_little_endian(semantic.data() + alpha_offset,
                               serialized_heat + 0.25);
    refresh_restart_checksum(semantic);
    fusion_source_state_v1* semantically_invalid = nullptr;
    expect_rejected(fusion_c_source_state_unpack(
                        semantic.data(), semantic.size(), extended.tag,
                        &semantically_invalid),
                    "checksum-valid but nonconserving inert heat is rejected");
    check(semantically_invalid == nullptr,
          "semantic restart rejection returns no context");

    fusion_source_state_v1* restored = nullptr;
    expect_ok(fusion_c_source_state_unpack(v2.data(), v2.size(), extended.tag,
                                            &restored),
              "unpack promoted version 2 context");
    check(restored != nullptr, "unpack v2 returns a context");
    snapshot_inert(extended.state, extended.cells, accepted_s, accepted_t,
                   cumulative, inert, accepted_time, epoch);
    std::vector<double> restored_s, restored_t;
    fusion_source_ledger_v1 restored_cumulative{};
    std::array<double, kSpecies> restored_inert{};
    double restored_time = 0.0;
    std::uint64_t restored_epoch = 0;
    snapshot_inert(restored, extended.cells, restored_s, restored_t,
                   restored_cumulative, restored_inert, restored_time,
                   restored_epoch);
    check(restored_s == accepted_s && restored_t == accepted_t &&
              same_ledger(restored_cumulative, cumulative) &&
              same_heat(restored_inert, inert) &&
              restored_time == accepted_time && restored_epoch == epoch,
          "unpacked v2 state matches accepted populations and all accounts");

    /* Exercise a negative extra-heat step after restart and compare the
     * complete future packet.  The higher-energy candidate needs signed
     * negative inert heat to close the same inventory identity. */
    const std::vector<double> higher = extended.initial_s;
    std::array<double, kSpecies> negative_inert{};
    negative_inert[static_cast<std::size_t>(kAlpha)] = -heat;
    std::uint64_t original_ticket = 0;
    std::uint64_t restored_ticket = 0;
    expect_ok(fusion_c_source_state_begin(extended.state, 0.125,
                                           &original_ticket),
              "begin future step on original v2 context");
    expect_ok(fusion_c_source_state_begin(restored, 0.125, &restored_ticket),
              "begin future step on restored v2 context");
    expect_ok(fusion_c_source_state_stage_inert(
                  extended.state, original_ticket, higher.data(),
                  extended.initial_t.data(), &step, negative_inert.data()),
              "stage future negative inert heat on original context");
    expect_ok(fusion_c_source_state_stage_inert(
                  restored, restored_ticket, higher.data(),
                  extended.initial_t.data(), &step, negative_inert.data()),
              "stage future negative inert heat on restored context");
    expect_ok(fusion_c_source_state_commit(extended.state, original_ticket),
              "commit future step on original v2 context");
    expect_ok(fusion_c_source_state_commit(restored, restored_ticket),
              "commit future step on restored v2 context");

    snapshot_inert(extended.state, extended.cells, accepted_s, accepted_t,
                   cumulative, inert, accepted_time, epoch);
    snapshot_inert(restored, extended.cells, restored_s, restored_t,
                   restored_cumulative, restored_inert, restored_time,
                   restored_epoch);
    check(restored_s == accepted_s && restored_t == accepted_t &&
              same_ledger(restored_cumulative, cumulative) &&
              same_heat(restored_inert, inert) && restored_time == accepted_time &&
              restored_epoch == epoch,
          "restart contexts follow the same future accepted trajectory");
    check(inert[static_cast<std::size_t>(kAlpha)] == 0.0,
          "positive and negative inert heat cancel cumulatively");
    check(pack_state(extended.state) == pack_state(restored),
          "restart contexts produce identical future v2 bytes");
    expect_rejected(fusion_c_source_state_snapshot(
                        restored, accepted_s.data(), accepted_t.data(),
                        &cumulative, &accepted_time, &epoch),
                    "old snapshot rejects restored promoted context");
    fusion_c_source_state_destroy(restored);
}

void test_snapshot_inert_failure_clears_nonkinetic_outputs() {
    constexpr int cells = 2;
    std::vector<double> accepted_s(static_cast<std::size_t>(kSpecies) * cells,
                                   -3.0);
    std::vector<double> accepted_t(static_cast<std::size_t>(kSpecies) * cells,
                                   -4.0);
    const std::vector<double> s_before = accepted_s;
    const std::vector<double> t_before = accepted_t;
    fusion_source_ledger_v1 cumulative{};
    cumulative.events_m3[0] = 1.0;
    std::array<double, kSpecies> inert{};
    inert.fill(2.0);
    double accepted_time = 5.0;
    std::uint64_t epoch = 6;
    expect_rejected(fusion_c_source_state_snapshot_inert(
                        nullptr, accepted_s.data(), accepted_t.data(),
                        &cumulative, inert.data(), &accepted_time, &epoch),
                    "snapshot_inert rejects a null context");
    check(accepted_s == s_before && accepted_t == t_before,
          "failed snapshot_inert leaves kinetic arrays untouched");
    check_zero_ledger(cumulative,
                      "failed snapshot_inert clears cumulative ledger");
    check(same_heat(inert, std::array<double, kSpecies>{}),
          "failed snapshot_inert clears inert output");
    check(accepted_time == 0.0 && epoch == 0,
          "failed snapshot_inert clears scalar outputs");
}

}  // namespace

int main() {
    try {
        test_signed_inert_heat_and_snapshot_contract();
        test_discard_and_rejected_extended_stage_do_not_promote();
        test_repeated_extended_stage_replaces_without_accumulation();
        test_v1_and_v2_restart_compatibility();
        test_snapshot_inert_failure_clears_nonkinetic_outputs();
        std::cout << "All fusion source-state inert tests passed\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "FAIL: " << ex.what() << '\n';
        return 1;
    }
}
