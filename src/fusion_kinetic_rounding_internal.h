#ifndef FUSION_KINETIC_ROUNDING_INTERNAL_H
#define FUSION_KINETIC_ROUNDING_INTERNAL_H
#include "fusion_two_component.h"
#include <cmath>
#include <limits>
namespace fusion_detail {
struct kinetic_rounding_budget { long double number = 0, energy = 0; };
// Round accumulated measured allowances outward; exact zero adds no allowance.
inline bool accumulate_rounding(long double value, long double& budget) {
    if (!std::isfinite(value) || value < 0 ||
        !std::isfinite(budget) || budget < 0) return false;
    if (value == 0) return true;
    const long double sum = budget + value;
    if (!std::isfinite(sum)) return false;
    const long double outward = std::nextafter(sum,
        std::numeric_limits<long double>::infinity());
    if (!std::isfinite(outward)) return false;
    budget = outward;
    return true;
}
int energy_fp_trial_precise(int cells, int baths, double dt_s,
    const double* edges_J, const double* old_number_m3,
    const double* bath_kT_J, const double* diffusion_J2_s,
    const double* birth_m3_s, const double* escape_s_inv,
    double thermalization_s_inv, double* trial_number_m3,
    double* heat_to_bath_J_m3, fusion_kinetic_ledger_v1* ledger,
    kinetic_rounding_budget* budget);
int two_component_trial_precise(int cells, int baths, double dt_s,
    const double* edges_J, const double* old_s_m3, const double* old_t_m3,
    const double* bath_kT_J, const double* diffusion_J2_s,
    const double* birth_s_m3_s, const double* birth_t_m3_s,
    const double* escape_s_inv, const double* transfer_s_inv,
    double* trial_s_m3, double* trial_t_m3, double* heat_to_bath_J_m3,
    fusion_two_component_ledger_v1* ledger, kinetic_rounding_budget* budget);
}
#endif
