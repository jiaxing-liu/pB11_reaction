#include "fusion_kinetic_rounding_internal.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>

int main() {
    long double accumulated = 0;
    assert(fusion_detail::accumulate_rounding(0, accumulated) && accumulated == 0);
    assert(!fusion_detail::accumulate_rounding(-1, accumulated) && accumulated == 0);
    assert(!fusion_detail::accumulate_rounding(
        std::numeric_limits<long double>::infinity(), accumulated));
    accumulated = std::numeric_limits<long double>::max();
    assert(!fusion_detail::accumulate_rounding(accumulated, accumulated));
    const double edges[] = {0.0, 1.0};
    const double zero[] = {0.0};
    const double birth[] = {2.02e-322};
    double s[1] = {7.0}, t[1] = {7.0};
    fusion_two_component_ledger_v1 ledger{};
    fusion_detail::kinetic_rounding_budget budget{7, 7};
    int status = fusion_detail::two_component_trial_precise(
        1, 0, 6.25e-6, edges, zero, zero, nullptr, nullptr,
        birth, zero, zero, zero, s, t, nullptr, &ledger, &budget);
    assert(status == PB11_STATUS_OK && s[0] == 0.0 && t[0] == 0.0);
    const long double exact_birth = static_cast<long double>(6.25e-6) * birth[0];
    assert(exact_birth > 0 && budget.number >= exact_birth);
    assert(budget.energy >= exact_birth / 2);
    auto ulps_above = [](long double value, int count) {
        while (count--) value = std::nextafter(value, std::numeric_limits<long double>::infinity());
        return value;
    };
    assert(budget.number <= ulps_above(exact_birth, 4));
    assert(budget.energy <= ulps_above(exact_birth / 2, 4));

    // Multiple tiny births exercise accumulation across cells and components.
    const double multi_edges[] = {0, 1, 2, 3};
    const double multi_zero[] = {0, 0, 0};
    const double multi_birth[] = {2.02e-322, 1.04e-322, 3.01e-322};
    double multi_s[3], multi_t[3];
    status = fusion_detail::two_component_trial_precise(3, 0, 6.25e-6,
        multi_edges, multi_zero, multi_zero, nullptr, nullptr,
        multi_birth, multi_birth, multi_zero, multi_zero,
        multi_s, multi_t, nullptr, &ledger, &budget);
    assert(status == PB11_STATUS_OK);
    long double number_loss = 0, energy_loss = 0;
    for (int i = 0; i < 3; ++i) {
        assert(multi_s[i] == 0 && multi_t[i] == 0);
        const long double loss = static_cast<long double>(6.25e-6) * multi_birth[i];
        number_loss += 2 * loss;
        energy_loss += 2 * loss * (i + 0.5L);
    }
    assert(budget.number >= number_loss && budget.energy >= energy_loss);
    assert(budget.number <= ulps_above(number_loss, 24));
    assert(budget.energy <= ulps_above(energy_loss, 24));
    // The real residual is the lost positive birth; a relative-only gate fails.
    assert(exact_birth > 1e-10L * exact_birth);

    const double old[] = {3.0}, normal_birth[] = {2.0};
    const double escape[] = {0.5}, transfer[] = {0.25};
    double public_s[1], public_t[1];
    fusion_two_component_ledger_v1 public_ledger{};
    status = fusion_detail::two_component_trial_precise(
        1, 0, 0.125, edges, old, old, nullptr, nullptr,
        normal_birth, normal_birth, escape, transfer,
        s, t, nullptr, &ledger, &budget);
    assert(status == PB11_STATUS_OK && budget.number == 0 && budget.energy == 0);
    assert(fusion_c_two_component_trial(1, 0, 0.125, edges, old, old,
        nullptr, nullptr, normal_birth, normal_birth, escape, transfer,
        public_s, public_t, nullptr, &public_ledger) == PB11_STATUS_OK);
    assert(std::memcmp(s, public_s, sizeof(s)) == 0);
    assert(std::memcmp(t, public_t, sizeof(t)) == 0);
    assert(std::memcmp(&ledger, &public_ledger, sizeof(ledger)) == 0);

    budget = {7, 7}; s[0] = t[0] = 7;
    status = fusion_detail::two_component_trial_precise(
        1, 0, -1, edges, old, old, nullptr, nullptr,
        normal_birth, normal_birth, escape, transfer,
        s, t, nullptr, &ledger, &budget);
    assert(status == PB11_STATUS_OUT_OF_RANGE && s[0] == 0 && t[0] == 0);
    assert(budget.number == 0 && budget.energy == 0);
    fusion_two_component_ledger_v1 empty{};
    assert(std::memcmp(&ledger, &empty, sizeof(ledger)) == 0);

    fusion_kinetic_ledger_v1 fp_ledger{};
    budget = {7, 7}; s[0] = 7;
    status = fusion_detail::energy_fp_trial_precise(1, 0, -1,
        edges, old, nullptr, nullptr, normal_birth, escape, 0,
        s, nullptr, &fp_ledger, &budget);
    assert(status == PB11_STATUS_OUT_OF_RANGE && s[0] == 0);
    assert(budget.number == 0 && budget.energy == 0);
}
