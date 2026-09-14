#include "../src/fusion_source_rounding_internal.h"
using fusion_detail::source_rounding::floor_source_rate;
using fusion_detail::source_rounding::source_roundoff_accumulate;

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace {

using Rounding = long double;
constexpr Rounding kHalfQuantum =
    Rounding(std::numeric_limits<double>::denorm_min()) / 2;

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& message) {
    ++checks;
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

bool same(Rounding left, Rounding right) {
    return left == right;
}

bool bounded_counter(Rounding value, Rounding exact) {
    return std::isfinite(value) && value >= exact && value < kHalfQuantum;
}

void test_actual_tail_and_optional_accounting() {
    // Captured He4-cell-214 packet from the D109 source replay.  Keep the
    // hexadecimal literal exact so this test cannot drift through decimal I/O.
    const Rounding amount = 0xb.843a57ab5777e7fp-1234L;
    const double dt = 1.25e-5;
    const Rounding exact_rate = amount / Rounding(dt);
    const Rounding denorm = Rounding(std::numeric_limits<double>::denorm_min());

    check(amount > 0 && amount < kHalfQuantum,
          "actual packet is positive and below half a binary64 quantum");
    check(exact_rate > 0 && exact_rate < kHalfQuantum,
          "actual source rate is positive and below half a binary64 quantum");

    double rate = -1.0;
    check(!floor_source_rate(amount, dt, rate),
          "strict default rejects the positive actual tail");
    check(rate == 0.0, "actual tail naturally converts to a zero double rate");

    Rounding lost = -1.0L;
    rate = -1.0;
    check(floor_source_rate(amount, dt, rate, &lost),
          "optional accounting accepts only the naturally unrepresentable tail");
    check(rate == 0.0 && same(lost, amount),
          "optional accounting returns zero rate and exact lost packet");
    check(std::abs(Rounding(rate) * Rounding(dt) - amount) <=
              2 * std::numeric_limits<double>::epsilon() * std::abs(amount) +
                  denorm * Rounding(dt),
          "optional actual tail still satisfies source-scale reconstruction bound");

    // The strict inequality is part of the optional contract.  Just below the
    // half quantum is eligible; equality is not.
    const Rounding below = std::nextafter(kHalfQuantum, Rounding(0));
    lost = -1.0L;
    rate = -1.0;
    check(floor_source_rate(below, 1.0, rate, &lost),
          "packet strictly below half quantum may be accounted as lost");
    check(rate == 0.0 && same(lost, below),
          "below-half packet is returned atomically as lost");

    lost = -1.0L;
    rate = -1.0;
    check(!floor_source_rate(kHalfQuantum, 1.0, rate, &lost),
          "packet exactly at half quantum is rejected by optional mode");
    check(lost == 0.0,
          "rejected half-quantum packet does not publish optional loss");
}

void test_representable_tiny_rate_and_large_dt() {
    const Rounding denorm = Rounding(std::numeric_limits<double>::denorm_min());

    // The packet itself is a long-double amount, but its exact rate is a
    // representable binary64 subnormal.  It must remain visible to callers.
    const double dt = 1.25e-5;
    const Rounding amount = 3 * denorm * Rounding(dt);
    const Rounding expected_rate = 3 * denorm;
    double rate = 0.0;
    Rounding lost = -1.0L;
    check(floor_source_rate(amount, dt, rate, &lost),
          "representable tiny source rate is accepted");
    check(rate == double(expected_rate) && rate != 0.0,
          "representable tiny source rate is preserved as nonzero");
    check(lost == 0.0,
          "representable tiny rate does not enter the unrepresented account");

    // A packet that is itself at least half a quantum cannot be hidden merely
    // because a very large dt makes its per-second rate underflow.
    const Rounding packet = denorm;
    const double huge_dt = std::numeric_limits<double>::max();
    rate = -1.0;
    lost = -1.0L;
    check(!floor_source_rate(packet, huge_dt, rate, &lost),
          "representable packet with underflowed huge-dt rate is rejected");
    check(rate == 0.0 && lost == 0.0,
          "huge-dt rejection exposes no positive packet as a zero source");
}

void test_roundoff_accumulator_boundaries() {
    const Rounding q = kHalfQuantum / 8;

    Rounding number = 0.0L;
    Rounding energy = 0.0L;
    check(source_roundoff_accumulate(q, 1.0L, number, energy),
          "accumulator accepts totals strictly below both half-quantum bounds");
    check(bounded_counter(number, q) && bounded_counter(energy, q),
          "first accepted loss conservatively bounds number and energy");

    const Rounding old_number = number;
    const Rounding old_energy = energy;
    // number would become 8q = half, so the strict N boundary rejects and
    // leaves both totals untouched.
    check(!source_roundoff_accumulate(7 * q, 1.0L, number, energy),
          "accumulator rejects aggregate number exactly at half quantum");
    check(same(number, old_number) && same(energy, old_energy),
          "number-boundary rejection is atomic");

    // This update keeps N below half (2q) but takes E to exactly half (8q).
    check(!source_roundoff_accumulate(q, 7.0L, number, energy),
          "accumulator rejects aggregate energy exactly at half quantum");
    check(same(number, old_number) && same(energy, old_energy),
          "energy-boundary rejection is atomic");

    number = 0.0L;
    energy = 0.0L;
    check(source_roundoff_accumulate(q, 1.0L, number, energy),
          "first aggregate contribution is accepted");
    check(source_roundoff_accumulate(q, 1.0L, number, energy),
          "aggregation remains valid while both totals stay below half");
    check(bounded_counter(number, 2 * q) && bounded_counter(energy, 2 * q),
          "successful aggregation conservatively bounds both totals");
    const Rounding aggregate_number = number;
    const Rounding aggregate_energy = energy;
    check(!source_roundoff_accumulate(6 * q, 1.0L, number, energy),
          "later aggregation crossing N boundary is rejected");
    check(same(number, aggregate_number) && same(energy, aggregate_energy),
          "crossing N boundary preserves prior aggregate");

    const Rounding before_zero_number = number;
    const Rounding before_zero_energy = energy;
    check(source_roundoff_accumulate(0.0L, 1.0L, number, energy),
          "zero loss is harmless while prior totals remain below half");
    check(same(number, before_zero_number) && same(energy, before_zero_energy),
          "zero loss does not inflate conservative counters");
}

void test_invalid_and_overflow_inputs() {
    const Rounding nan = std::numeric_limits<Rounding>::quiet_NaN();
    const Rounding inf = std::numeric_limits<Rounding>::infinity();
    double rate = 7.0;
    Rounding lost = 7.0L;

    check(!floor_source_rate(-1.0L, 1.0, rate, &lost),
          "negative packet is rejected");
    check(!floor_source_rate(nan, 1.0, rate, &lost),
          "NaN packet is rejected");
    check(!floor_source_rate(1.0L, nan, rate, &lost),
          "NaN dt is rejected");
    check(!floor_source_rate(1.0L, 0.0, rate, &lost),
          "zero dt is rejected");
    check(!floor_source_rate(1.0L, -1.0, rate, &lost),
          "negative dt is rejected");
    check(!floor_source_rate(inf, 1.0, rate, &lost),
          "infinite packet is rejected");
    check(!floor_source_rate(1.0L, std::numeric_limits<double>::infinity(), rate,
                             &lost),
          "infinite dt is rejected");

    const Rounding over_rate = Rounding(std::numeric_limits<double>::max()) * 2;
    check(!floor_source_rate(over_rate, 1.0, rate, &lost),
          "finite long-double packet whose rate overflows binary64 is rejected");
    check(!floor_source_rate(1.0L,
                             std::numeric_limits<double>::denorm_min(), rate,
                             &lost),
          "finite packet with binary64-overflowing rate is rejected");

    Rounding number = 0.0L;
    Rounding energy = 0.0L;
    Rounding bad_number = nan;
    check(!source_roundoff_accumulate(-1.0L, 1.0L, number, energy),
          "negative lost amount is rejected");
    check(!source_roundoff_accumulate(nan, 1.0L, number, energy),
          "NaN lost amount is rejected");
    check(!source_roundoff_accumulate(1.0L, nan, number, energy),
          "NaN center is rejected");
    check(!source_roundoff_accumulate(1.0L, 0.0L, number, energy),
          "zero center is rejected");
    check(!source_roundoff_accumulate(1.0L, 1.0L, bad_number, energy),
          "NaN number total is rejected");

    number = 0.0L;
    energy = 0.0L;
    check(!source_roundoff_accumulate(std::numeric_limits<Rounding>::max(),
                                      std::numeric_limits<Rounding>::max(),
                                      number, energy),
          "overflowing roundoff contribution is rejected");
    check(number == 0.0L && energy == 0.0L,
          "overflow rejection is atomic");
}

void test_regular_finite_roundtrip() {
    const Rounding amount = 0.1234567890123456789L;
    const double dt = 1.25e-5;
    double rate = 0.0;
    Rounding lost = -1.0L;
    check(floor_source_rate(amount, dt, rate, &lost),
          "finite regular packet completes source-rate roundtrip");
    check(std::isfinite(rate) && rate > 0.0 && lost == 0.0,
          "finite regular roundtrip returns finite positive rate without loss");
    const Rounding recovered = Rounding(rate) * Rounding(dt);
    const Rounding allowance =
        2 * std::numeric_limits<double>::epsilon() * std::abs(amount) +
        Rounding(std::numeric_limits<double>::denorm_min()) * Rounding(dt);
    check(std::abs(recovered - amount) <= allowance,
          "finite regular roundtrip obeys source-scale reconstruction bound");
}

}  // namespace

int main() {
    test_actual_tail_and_optional_accounting();
    test_representable_tiny_rate_and_large_dt();
    test_roundoff_accumulator_boundaries();
    test_invalid_and_overflow_inputs();
    test_regular_finite_roundtrip();

    if (failures != 0) {
        std::cerr << "FAIL source rounding: " << failures << " of " << checks
                  << " checks\n";
        return 1;
    }
    std::cout << "PASS source rounding: " << checks << " checks\n";
    return 0;
}
