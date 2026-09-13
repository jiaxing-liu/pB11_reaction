#include "fusion_source_state.h"
#include "fusion_kinetics.h"
#include "fusion_nuclear_data.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kElectronVoltJ = 1.602176634e-19;
constexpr double kMeVJ = 1.0e6 * kElectronVoltJ;
constexpr int kSpecies = FUSION_SPECIES_COUNT;
constexpr int kBaths = 7;
constexpr int kLedgerHeat = kSpecies * kBaths;
constexpr int kAlpha = FUSION_HELIUM4;
constexpr int kDeuteron = FUSION_DEUTERON;
constexpr int kTriton = FUSION_TRITON;

void check(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void check_status(int status, const std::string& message) {
    check(status == PB11_STATUS_OK,
          message + " (status=" + std::to_string(status) + ")");
}

void check_rejected(int status, const std::string& message) {
    check(status != PB11_STATUS_OK,
          message + " unexpectedly returned OK");
}

using StateArray = std::array<double, kSpecies * 2>;

bool same_array(const double* actual, const double* expected, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        if (actual[i] != expected[i]) {
            return false;
        }
    }
    return true;
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
    for (int i = 0; i < kLedgerHeat; ++i) {
        if (actual.heat_to_bath_J_m3[i] != expected.heat_to_bath_J_m3[i]) {
            return false;
        }
    }
    return actual.neutron_number_m3 == expected.neutron_number_m3 &&
           actual.neutron_energy_J_m3 == expected.neutron_energy_J_m3;
}

void check_zero_ledger(const fusion_source_ledger_v1& ledger,
                       const std::string& label) {
    fusion_source_ledger_v1 zero{};
    check(same_ledger(ledger, zero), label + " is zero");
}

struct BasicState {
    int cells = 2;
    std::array<double, 3> edges{{0.0, 2.0 * kMeVJ, 4.0 * kMeVJ}};
    StateArray initial_s{};
    StateArray initial_t{};
    double initial_time = 3.25;
    std::uint64_t tag = UINT64_C(0x1122334455667788);
    fusion_source_state_v1* state = nullptr;

    BasicState() {
        initial_s[static_cast<std::size_t>(kAlpha * cells)] = 1.0e12;
    }

    ~BasicState() { fusion_c_source_state_destroy(state); }

    BasicState(const BasicState&) = delete;
    BasicState& operator=(const BasicState&) = delete;
};

void create_basic(BasicState& basic) {
    check_status(fusion_c_source_state_create(
                     basic.cells, basic.edges.data(), basic.initial_s.data(),
                     basic.initial_t.data(), basic.initial_time, basic.tag,
                     &basic.state),
                 "create two-cell source state");
    check(basic.state != nullptr, "create returns a context");
}

void snapshot(const BasicState& basic, StateArray& s, StateArray& t,
              fusion_source_ledger_v1& cumulative, double& time,
              std::uint64_t& epoch) {
    check_status(fusion_c_source_state_snapshot(
                     basic.state, s.data(), t.data(), &cumulative, &time,
                     &epoch),
                 "snapshot accepted source state");
}

fusion_source_ledger_v1 zero_step() { return fusion_source_ledger_v1{}; }

void test_unrepresentable_restart_anchor() {
    double edges[2]{0.0, 1e-300};
    double s[6]{1e-300, 0, 0, 0, 0, 0}, t[6]{};
    fusion_source_state_v1* state = nullptr;
    check_rejected(fusion_c_source_state_create(1, edges, s, t, 0, 123, &state),
                   "unrepresentable initial energy anchor is rejected at create");
    check(state == nullptr, "failed anchor create clears context");
}

void test_create_and_snapshot() {
    BasicState basic;
    create_basic(basic);
    int cells = 0;
    check_status(fusion_c_source_state_cells(basic.state, &cells),
                 "query initial source-state cell count");
    check(cells == basic.cells, "initial cell-count query matches create");
    cells = 123;
    check_rejected(fusion_c_source_state_cells(nullptr, &cells),
                   "cell-count query rejects null state");
    check(cells == 0, "null-state cell-count query clears output");
    check_rejected(fusion_c_source_state_cells(basic.state, nullptr),
                   "cell-count query rejects null output");

    StateArray s{}, t{};
    fusion_source_ledger_v1 cumulative{};
    double time = -1.0;
    std::uint64_t epoch = UINT64_MAX;
    snapshot(basic, s, t, cumulative, time, epoch);
    check(same_array(s.data(), basic.initial_s.data(), s.size()),
          "snapshot returns initial S populations");
    check(same_array(t.data(), basic.initial_t.data(), t.size()),
          "snapshot returns initial T populations");
    check_zero_ledger(cumulative, "initial cumulative ledger");
    check(time == basic.initial_time, "snapshot returns initial time");
    check(epoch == 0, "initial accepted epoch is zero");
}

void test_begin_stage_commit_once() {
    BasicState basic;
    create_basic(basic);
    const auto step = zero_step();
    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.25, &ticket),
                 "begin no-op step");
    check(ticket != 0, "begin returns a nonzero ticket");
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, basic.initial_s.data(),
                     basic.initial_t.data(), &step),
                 "stage no-op step");
    check_status(fusion_c_source_state_commit(basic.state, ticket),
                 "commit no-op step once");

    StateArray s{}, t{};
    fusion_source_ledger_v1 cumulative{};
    double time = 0.0;
    std::uint64_t epoch = 0;
    snapshot(basic, s, t, cumulative, time, epoch);
    check(same_array(s.data(), basic.initial_s.data(), s.size()),
          "no-op commit leaves S unchanged");
    check(same_array(t.data(), basic.initial_t.data(), t.size()),
          "no-op commit leaves T unchanged");
    check(time == basic.initial_time + 0.25, "commit advances time once");
    check(epoch == 1, "commit increments epoch once");
    check_zero_ledger(cumulative, "no-op commit cumulative ledger");
    check_rejected(fusion_c_source_state_commit(basic.state, ticket),
                   "repeated commit is rejected");

    snapshot(basic, s, t, cumulative, time, epoch);
    check(time == basic.initial_time + 0.25 && epoch == 1,
          "repeated commit does not advance accepted state");
}

void test_stale_superseded_ticket() {
    BasicState basic;
    create_basic(basic);
    const auto step = zero_step();
    std::uint64_t old_ticket = 0;
    std::uint64_t current_ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.1, &old_ticket),
                 "begin first superseded trial");
    check_status(fusion_c_source_state_begin(basic.state, 0.2, &current_ticket),
                 "begin second superseding trial");
    check(old_ticket != current_ticket, "superseding begin returns new ticket");
    check_rejected(fusion_c_source_state_stage(
                       basic.state, old_ticket, basic.initial_s.data(),
                       basic.initial_t.data(), &step),
                   "stale ticket cannot stage");
    check_status(fusion_c_source_state_stage(
                     basic.state, current_ticket, basic.initial_s.data(),
                     basic.initial_t.data(), &step),
                 "current ticket stages after stale rejection");
    check_rejected(fusion_c_source_state_commit(basic.state, old_ticket),
                   "stale ticket cannot commit");
    check_status(fusion_c_source_state_commit(basic.state, current_ticket),
                 "current superseding ticket commits");

    StateArray s{}, t{};
    fusion_source_ledger_v1 cumulative{};
    double time = 0.0;
    std::uint64_t epoch = 0;
    snapshot(basic, s, t, cumulative, time, epoch);
    check(time == basic.initial_time + 0.2 && epoch == 1,
          "superseded trial uses only current dt");
}

void test_failed_stage_invalidates_pending() {
    BasicState basic;
    create_basic(basic);
    const auto step = zero_step();
    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.1, &ticket),
                 "begin trial before failed stage");
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, basic.initial_s.data(),
                     basic.initial_t.data(), &step),
                 "stage valid candidate before failed stage");
    StateArray bad_s = basic.initial_s;
    bad_s[0] = -1.0;
    check_rejected(fusion_c_source_state_stage(
                       basic.state, ticket, bad_s.data(),
                       basic.initial_t.data(), &step),
                   "invalid stage is rejected");
    check_rejected(fusion_c_source_state_commit(basic.state, ticket),
                   "failed stage invalidates pending commit");

    StateArray s{}, t{};
    fusion_source_ledger_v1 cumulative{};
    double time = 0.0;
    std::uint64_t epoch = 0;
    snapshot(basic, s, t, cumulative, time, epoch);
    check(same_array(s.data(), basic.initial_s.data(), s.size()) &&
              same_array(t.data(), basic.initial_t.data(), t.size()) &&
              time == basic.initial_time && epoch == 0,
          "failed stage leaves accepted state unchanged");
}

void test_discard_unchanged() {
    BasicState basic;
    create_basic(basic);
    auto step = zero_step();
    step.escaped_number_m3[kAlpha] = 2.0e11;
    step.escaped_energy_J_m3[kAlpha] = 2.0e11 * kMeVJ;
    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.5, &ticket),
                 "begin trial to discard");
    StateArray candidate = basic.initial_s;
    candidate[static_cast<std::size_t>(kAlpha * basic.cells)] -= 2.0e11;
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, candidate.data(),
                     basic.initial_t.data(), &step),
                 "stage candidate to discard");
    check_status(fusion_c_source_state_discard(basic.state, ticket),
                 "discard pending candidate");
    check_rejected(fusion_c_source_state_commit(basic.state, ticket),
                   "discarded ticket cannot commit");

    StateArray s{}, t{};
    fusion_source_ledger_v1 cumulative{};
    double time = 0.0;
    std::uint64_t epoch = 0;
    snapshot(basic, s, t, cumulative, time, epoch);
    check(same_array(s.data(), basic.initial_s.data(), s.size()) &&
              same_array(t.data(), basic.initial_t.data(), t.size()) &&
              time == basic.initial_time && epoch == 0,
          "discard leaves accepted state unchanged");
    check_zero_ledger(cumulative, "discard leaves cumulative ledger unchanged");
}

void test_repeated_nonlinear_stage_replaces() {
    BasicState basic;
    create_basic(basic);
    auto first_step = zero_step();
    first_step.escaped_number_m3[kAlpha] = 1.0e11;
    first_step.escaped_energy_J_m3[kAlpha] = 1.0e11 * kMeVJ;
    auto second_step = zero_step();
    second_step.escaped_number_m3[kAlpha] = 2.0e11;
    second_step.escaped_energy_J_m3[kAlpha] = 2.0e11 * kMeVJ;
    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.1, &ticket),
                 "begin repeated nonlinear trial");
    StateArray first = basic.initial_s;
    first[static_cast<std::size_t>(kAlpha * basic.cells)] -= 1.0e11;
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, first.data(), basic.initial_t.data(),
                     &first_step),
                 "stage first nonlinear evaluation");
    StateArray second = basic.initial_s;
    second[static_cast<std::size_t>(kAlpha * basic.cells)] -= 2.0e11;
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, second.data(),
                     basic.initial_t.data(), &second_step),
                 "stage replacement nonlinear evaluation");
    check_status(fusion_c_source_state_commit(basic.state, ticket),
                 "commit replacement nonlinear evaluation");

    StateArray s{}, t{};
    fusion_source_ledger_v1 cumulative{};
    double time = 0.0;
    std::uint64_t epoch = 0;
    snapshot(basic, s, t, cumulative, time, epoch);
    check(same_array(s.data(), second.data(), s.size()),
          "repeated nonlinear stage commits replacement state");
    check(s[static_cast<std::size_t>(kAlpha * basic.cells)] ==
              basic.initial_s[static_cast<std::size_t>(kAlpha * basic.cells)] -
                  2.0e11,
          "repeated nonlinear stage does not add first candidate");
    check(epoch == 1 && time == basic.initial_time + 0.1,
          "replacement stage commits one time step");
}

struct FpTrial {
    std::array<double, 2> trial{};
    std::array<double, 1> heat{};
    fusion_kinetic_ledger_v1 ledger{};
};

FpTrial actual_alpha_fp_trial() {
    const double edges[3]{0.0, 2.0 * kMeVJ, 4.0 * kMeVJ};
    const double old[2]{1.0e12, 0.0};
    const double kT[1]{0.1 * kMeVJ};
    const double diffusion[1]{1.0e-26};
    const double birth[2]{0.0, 0.0};
    const double escape[2]{0.1, 0.2};
    FpTrial result;
    check_status(fusion_c_energy_fp_trial(
                     2, 1, 0.01, edges, old, kT, diffusion, birth, escape,
                     0.0, result.trial.data(), result.heat.data(),
                     &result.ledger),
                 "actual alpha finite-volume trial");
    check(std::isfinite(result.heat[0]) && result.heat[0] != 0.0,
          "actual alpha trial has signed heat row " +
              std::to_string(result.heat[0]));
    check(result.ledger.escaped_number_m3 > 0.0 &&
              result.ledger.escaped_energy_J_m3 > 0.0,
          "actual alpha trial produces escape ledger");
    return result;
}

void test_actual_fp_integration() {
    BasicState basic;
    create_basic(basic);
    const FpTrial fp = actual_alpha_fp_trial();
    StateArray trial_s = basic.initial_s;
    trial_s[static_cast<std::size_t>(kAlpha * basic.cells)] = fp.trial[0];
    trial_s[static_cast<std::size_t>(kAlpha * basic.cells + 1)] = fp.trial[1];
    StateArray trial_t{};
    fusion_source_ledger_v1 step{};
    step.escaped_number_m3[kAlpha] = fp.ledger.escaped_number_m3;
    step.escaped_energy_J_m3[kAlpha] = fp.ledger.escaped_energy_J_m3;
    step.heat_to_bath_J_m3[kAlpha * kBaths] = fp.heat[0];

    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.01, &ticket),
                 "begin actual FP state step");
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, trial_s.data(), trial_t.data(), &step),
                 "stage actual FP state step");
    check_status(fusion_c_source_state_commit(basic.state, ticket),
                 "commit actual FP state step");

    StateArray accepted_s{}, accepted_t{};
    fusion_source_ledger_v1 cumulative{};
    double time = 0.0;
    std::uint64_t epoch = 0;
    snapshot(basic, accepted_s, accepted_t, cumulative, time, epoch);
    check(same_array(accepted_s.data(), trial_s.data(), accepted_s.size()),
          "accepted state equals actual FP trial");
    check(same_array(accepted_t.data(), trial_t.data(), accepted_t.size()),
          "actual FP integration keeps T component empty");
    check(cumulative.escaped_number_m3[kAlpha] ==
              fp.ledger.escaped_number_m3 &&
              cumulative.escaped_energy_J_m3[kAlpha] ==
                  fp.ledger.escaped_energy_J_m3,
          "state escape ledger uses FP particle and energy amounts");
    check(cumulative.heat_to_bath_J_m3[kAlpha * kBaths] == fp.heat[0],
          "state heat ledger uses He4 to electron bath row");
    check(time == basic.initial_time + 0.01 && epoch == 1,
          "actual FP step advances accepted time and epoch once");
}

struct DtLedgerCase {
    std::array<double, 2> edges{{0.0, 2.0 * kMeVJ}};
    std::array<double, kSpecies> initial_s{};
    std::array<double, kSpecies> initial_t{};
    std::array<double, kSpecies> trial_s{};
    std::array<double, kSpecies> trial_t{};
    fusion_source_ledger_v1 step{};
    double alpha_energy_J = 0.0;
    double neutron_energy_J = 0.0;
    double q_J = 0.0;
};

DtLedgerCase make_dt_case() {
    DtLedgerCase result;
    fusion_nuclear_channel_v1 channel{};
    check_status(fusion_c_nuclear_channel(FUSION_DT_ALPHAN, &channel),
                 "read new DT nuclear channel");
    result.q_J = channel.q_J;
    check(result.q_J > 0.0, "new DT Q is positive");

    /* One DT event consumes one thermal D and one thermal T.  Their
     * supplied kinetic energies are part of the available product energy. */
    const double consumed_d = 0.5 * kMeVJ;
    const double consumed_t = 0.25 * kMeVJ;
    result.alpha_energy_J = 1.0 * kMeVJ;
    result.neutron_energy_J = result.q_J + consumed_d + consumed_t -
                              result.alpha_energy_J;
    check(result.neutron_energy_J > 0.0,
          "DT one-event neutron remainder is positive");

    result.trial_s[kAlpha] = 1.0;
    result.step.events_m3[FUSION_DT_ALPHAN] = 1.0;
    result.step.nuclear_born_number_m3[kAlpha] = 1.0;
    result.step.nuclear_born_energy_J_m3[kAlpha] = result.alpha_energy_J;
    result.step.neutron_number_m3 = 1.0;
    result.step.neutron_energy_J_m3 = result.neutron_energy_J;
    result.step.thermal_consumed_number_m3[kDeuteron] = 1.0;
    result.step.thermal_consumed_number_m3[kTriton] = 1.0;
    result.step.thermal_consumed_energy_J_m3[kDeuteron] = consumed_d;
    result.step.thermal_consumed_energy_J_m3[kTriton] = consumed_t;
    return result;
}

void test_nuclear_dt_ledger_validation() {
    DtLedgerCase data = make_dt_case();
    fusion_source_state_v1* state = nullptr;
    check_status(fusion_c_source_state_create(
                     1, data.edges.data(), data.initial_s.data(),
                     data.initial_t.data(), 0.0, UINT64_C(0xD7), &state),
                 "create empty state for DT ledger");
    check(state != nullptr, "DT ledger context exists");

    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(state, 0.37, &ticket),
                 "begin DT ledger step");
    check_status(fusion_c_source_state_stage(
                     state, ticket, data.trial_s.data(), data.trial_t.data(),
                     &data.step),
                 "stage consistent DT one-event ledger");
    check_status(fusion_c_source_state_commit(state, ticket),
                 "commit consistent DT one-event ledger");
    std::array<double, kSpecies> accepted_s{}, accepted_t{};
    fusion_source_ledger_v1 cumulative{};
    double accepted_time = 0.0;
    std::uint64_t accepted_epoch = 0;
    check_status(fusion_c_source_state_snapshot(
                     state, accepted_s.data(), accepted_t.data(), &cumulative,
                     &accepted_time, &accepted_epoch),
                 "snapshot committed DT one-event ledger");
    check(same_array(accepted_s.data(), data.trial_s.data(), kSpecies) &&
              same_array(accepted_t.data(), data.trial_t.data(), kSpecies),
          "DT one-event snapshot stores newborn alpha in S");
    check(same_ledger(cumulative, data.step),
          "DT one-event cumulative ledger stores every explicit field");
    check(accepted_time == 0.37 && accepted_epoch == 1,
          "DT one-event snapshot stores time and epoch");
    fusion_c_source_state_destroy(state);

    /* Wrong reaction stoichiometry is rejected even though the fast alpha
     * population and all other ledger fields are unchanged. */
    data = make_dt_case();
    data.step.thermal_consumed_number_m3[kTriton] = 2.0;
    check_status(fusion_c_source_state_create(
                     1, data.edges.data(), data.initial_s.data(),
                     data.initial_t.data(), 0.0, UINT64_C(0xD8), &state),
                 "create state for wrong DT stoichiometry");
    check_status(fusion_c_source_state_begin(state, 0.2, &ticket),
                 "begin wrong stoichiometry step");
    check_rejected(fusion_c_source_state_stage(
                       state, ticket, data.trial_s.data(), data.trial_t.data(),
                       &data.step),
                   "wrong DT stoichiometry is rejected");
    check_rejected(fusion_c_source_state_commit(state, ticket),
                   "wrong DT stoichiometry invalidates trial");
    fusion_c_source_state_destroy(state);

    /* The new Q is tied to the canonical masses.  Perturbing only the neutron
     * kinetic energy must fail the reaction energy identity. */
    data = make_dt_case();
    data.step.neutron_energy_J_m3 += 0.5 * kMeVJ;
    check_status(fusion_c_source_state_create(
                     1, data.edges.data(), data.initial_s.data(),
                     data.initial_t.data(), 0.0, UINT64_C(0xD9), &state),
                 "create state for wrong DT Q");
    check_status(fusion_c_source_state_begin(state, 0.2, &ticket),
                 "begin wrong Q step");
    check_rejected(fusion_c_source_state_stage(
                       state, ticket, data.trial_s.data(), data.trial_t.data(),
                       &data.step),
                   "wrong DT Q energy is rejected");
    check_rejected(fusion_c_source_state_commit(state, ticket),
                   "wrong DT Q invalidates trial");
    fusion_c_source_state_destroy(state);
}

struct PackedState {
    std::vector<unsigned char> bytes;
    std::size_t required = 0;
};

std::uint64_t restart_checksum(const unsigned char* data, std::size_t count) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (std::size_t i = 0; i < count; ++i) {
        hash ^= static_cast<std::uint64_t>(data[i]);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

void put_u64_little_endian(unsigned char* data, std::uint64_t value) {
    for (int i = 0; i < 8; ++i) {
        data[i] = static_cast<unsigned char>(value & UINT64_C(0xff));
        value >>= 8;
    }
}

void refresh_restart_checksum(std::vector<unsigned char>& bytes) {
    check(bytes.size() >= sizeof(std::uint64_t),
          "restart packet has room for checksum");
    put_u64_little_endian(bytes.data() + bytes.size() - sizeof(std::uint64_t),
                          restart_checksum(bytes.data(),
                                           bytes.size() - sizeof(std::uint64_t)));
}

PackedState pack_accepted(const BasicState& basic) {
    PackedState packed;
    check_status(fusion_c_source_state_pack_size(basic.state,
                                                 &packed.required),
                 "get accepted state pack size");
    check(packed.required > 0, "accepted pack size is positive");
    packed.bytes.resize(packed.required);
    std::size_t written = 0;
    check_status(fusion_c_source_state_pack(
                     basic.state, packed.bytes.data(), packed.bytes.size(),
                     &written),
                 "pack accepted state");
    check(written == packed.required, "pack writes required bytes");
    return packed;
}

void test_pack_restart_and_rejections() {
    BasicState basic;
    create_basic(basic);
    const auto step = zero_step();
    std::uint64_t ticket = 0;
    check_status(fusion_c_source_state_begin(basic.state, 0.2, &ticket),
                 "begin pending pack trial");
    StateArray candidate = basic.initial_s;
    check_status(fusion_c_source_state_stage(
                     basic.state, ticket, candidate.data(),
                     basic.initial_t.data(), &step),
                 "stage pending pack trial");

    std::size_t required = 123;
    check_rejected(fusion_c_source_state_pack_size(basic.state, &required),
                   "pack size rejects pending state");
    check(required == 0, "pending pack size clears required bytes");
    std::array<unsigned char, 8> pending_buffer{};
    std::size_t written = 123;
    check_rejected(fusion_c_source_state_pack(
                       basic.state, pending_buffer.data(),
                       pending_buffer.size(), &written),
                   "pack rejects pending state");
    check(written == 0, "pending pack clears bytes written");
    check_status(fusion_c_source_state_discard(basic.state, ticket),
                 "discard pending pack trial");

    PackedState packed = pack_accepted(basic);
    std::vector<unsigned char> short_buffer(packed.required - 1);
    written = 123;
    check_rejected(fusion_c_source_state_pack(
                       basic.state, short_buffer.data(), short_buffer.size(),
                       &written),
                   "pack rejects truncated output capacity");
    check(written == 0, "short output capacity clears bytes written");

    fusion_source_state_v1* restored = nullptr;
    check_status(fusion_c_source_state_unpack(
                     packed.bytes.data(), packed.bytes.size(), basic.tag,
                     &restored),
                 "unpack accepted state with matching tag");
    check(restored != nullptr, "unpack returns a new context");
    StateArray original_s{}, original_t{}, restored_s{}, restored_t{};
    fusion_source_ledger_v1 original_ledger{}, restored_ledger{};
    double original_time = 0.0, restored_time = 0.0;
    std::uint64_t original_epoch = 0, restored_epoch = 0;
    snapshot(basic, original_s, original_t, original_ledger, original_time,
             original_epoch);
    BasicState restored_view;
    restored_view.state = restored;
    snapshot(restored_view, restored_s, restored_t, restored_ledger,
             restored_time, restored_epoch);
    check(same_array(restored_s.data(), original_s.data(), restored_s.size()) &&
              same_array(restored_t.data(), original_t.data(), restored_t.size()),
          "unpacked accepted populations equal original");
    check(same_ledger(restored_ledger, original_ledger) &&
              restored_time == original_time && restored_epoch == original_epoch,
          "unpacked accepted ledger, time and epoch equal original");
    int restored_cells = 0;
    check_status(fusion_c_source_state_cells(restored, &restored_cells),
                 "query unpacked source-state cell count");
    check(restored_cells == basic.cells,
          "unpacked cell-count query matches original");

    const auto restored_step = zero_step();
    std::uint64_t restored_ticket = 0;
    check_status(fusion_c_source_state_begin(restored, 0.125,
                                             &restored_ticket),
                 "begin next step after unpack");
    check_status(fusion_c_source_state_stage(
                     restored, restored_ticket, restored_s.data(),
                     restored_t.data(), &restored_step),
                 "stage next step after unpack");
    check_status(fusion_c_source_state_commit(restored, restored_ticket),
                 "commit next step after unpack");
    snapshot(restored_view, restored_s, restored_t, restored_ledger,
             restored_time, restored_epoch);
    check(restored_time == original_time + 0.125 && restored_epoch ==
              original_epoch + 1,
          "unpacked context continues from same accepted step");
    restored_view.state = nullptr;
    fusion_c_source_state_destroy(restored);

    fusion_source_state_v1* mismatch = nullptr;
    check_rejected(fusion_c_source_state_unpack(
                       packed.bytes.data(), packed.bytes.size(),
                       basic.tag + 1, &mismatch),
                   "unpack rejects model tag mismatch");
    check(mismatch == nullptr, "tag mismatch returns no context");

    std::vector<unsigned char> corrupt = packed.bytes;
    corrupt[corrupt.size() / 2] ^= 0x01u;
    fusion_source_state_v1* corrupted = nullptr;
    check_rejected(fusion_c_source_state_unpack(
                       corrupt.data(), corrupt.size(), basic.tag, &corrupted),
                   "unpack rejects checksum corruption");
    check(corrupted == nullptr, "corruption returns no context");

    fusion_source_state_v1* truncated = nullptr;
    check_rejected(fusion_c_source_state_unpack(
                       packed.bytes.data(), packed.bytes.size() - 1, basic.tag,
                       &truncated),
                   "unpack rejects truncated input");
    check(truncated == nullptr, "truncation returns no context");

    /* Recompute the specified FNV-1a checksum after changing the serialized
     * time field.  This reaches semantic validation after the corruption
     * check, rather than merely exercising checksum failure again. */
    std::vector<unsigned char> semantic = packed.bytes;
    constexpr std::size_t time_offset = 6 * sizeof(std::uint64_t);
    check(semantic.size() > time_offset + sizeof(std::uint64_t),
          "restart packet exposes a serialized time field");
    put_u64_little_endian(semantic.data() + time_offset,
                          UINT64_C(0x7ff8000000000000));
    refresh_restart_checksum(semantic);
    fusion_source_state_v1* invalid = nullptr;
    check_rejected(fusion_c_source_state_unpack(
                       semantic.data(), semantic.size(), basic.tag, &invalid),
                   "unpack rejects semantically invalid restart");
    check(invalid == nullptr, "semantic invalidity returns no context");
}

void test_two_independent_contexts() {
    BasicState first;
    BasicState second;
    first.tag = UINT64_C(0xA1);
    second.tag = UINT64_C(0xB2);
    create_basic(first);
    create_basic(second);
    const auto step = zero_step();
    std::uint64_t first_ticket = 0, second_ticket = 0;
    check_status(fusion_c_source_state_begin(first.state, 0.1, &first_ticket),
                 "begin first independent context");
    check_status(fusion_c_source_state_begin(second.state, 0.4, &second_ticket),
                 "begin second independent context");
    check_status(fusion_c_source_state_stage(
                     first.state, first_ticket, first.initial_s.data(),
                     first.initial_t.data(), &step),
                 "stage first independent context");
    check_status(fusion_c_source_state_stage(
                     second.state, second_ticket, second.initial_s.data(),
                     second.initial_t.data(), &step),
                 "stage second independent context");
    check_status(fusion_c_source_state_commit(first.state, first_ticket),
                 "commit first independent context");
    check_status(fusion_c_source_state_commit(second.state, second_ticket),
                 "commit second independent context");

    StateArray first_s{}, first_t{}, second_s{}, second_t{};
    fusion_source_ledger_v1 first_ledger{}, second_ledger{};
    double first_time = 0.0, second_time = 0.0;
    std::uint64_t first_epoch = 0, second_epoch = 0;
    snapshot(first, first_s, first_t, first_ledger, first_time, first_epoch);
    snapshot(second, second_s, second_t, second_ledger, second_time,
             second_epoch);
    check(first_time == first.initial_time + 0.1 && first_epoch == 1,
          "first context keeps its own time and epoch");
    check(second_time == second.initial_time + 0.4 && second_epoch == 1,
          "second context keeps its own time and epoch");
    check(same_array(first_s.data(), first.initial_s.data(), first_s.size()) &&
              same_array(second_s.data(), second.initial_s.data(),
                         second_s.size()),
          "independent contexts keep independent populations");
}

}  // namespace

int main() {
    try {
        test_unrepresentable_restart_anchor();
        test_create_and_snapshot();
        test_begin_stage_commit_once();
        test_stale_superseded_ticket();
        test_failed_stage_invalidates_pending();
        test_discard_unchanged();
        test_repeated_nonlinear_stage_replaces();
        test_actual_fp_integration();
        test_nuclear_dt_ledger_validation();
        test_pack_restart_and_rejections();
        test_two_independent_contexts();
        std::cout << "All fusion source-state tests passed\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "FAIL: " << ex.what() << '\n';
        return 1;
    }
}
