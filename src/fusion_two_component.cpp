#include "fusion_two_component.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

using Real = long double;

bool finite_nonnegative(double value) {
    return std::isfinite(value) && value >= 0.0;
}

bool finite_positive(double value) {
    return std::isfinite(value) && value > 0.0;
}

/* Keep the same representability rule as fusion_c_energy_fp_trial's ledger
 * conversion: a nonzero long-double quantity must remain nonzero in double.
 * This also rejects overflow and all nonfinite values before an output is
 * committed. */
bool put_double(Real value, double& destination) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<Real>(std::numeric_limits<double>::max())) {
        return false;
    }
    destination = static_cast<double>(value);
    return value == 0.0L || destination != 0.0;
}

void clear_outputs(int cells, int baths, double* trial_s, double* trial_t,
                   double* heat, fusion_two_component_ledger_v1* ledger) noexcept {
    if (trial_s != nullptr && cells > 0) {
        std::fill(trial_s, trial_s + cells, 0.0);
    }
    if (trial_t != nullptr && cells > 0) {
        std::fill(trial_t, trial_t + cells, 0.0);
    }
    if (heat != nullptr && baths > 0) {
        std::fill(heat, heat + baths, 0.0);
    }
    if (ledger != nullptr) {
        *ledger = {};
    }
}

}  // namespace

extern "C" int fusion_c_two_component_trial(
    int cells, int baths, double dt_s, const double* edges_J,
    const double* old_s_m3, const double* old_t_m3, const double* bath_kT_J,
    const double* diffusion_J2_s, const double* birth_s_m3_s,
    const double* birth_t_m3_s, const double* escape_s_inv,
    const double* transfer_s_inv, double* trial_s_m3, double* trial_t_m3,
    double* heat_to_bath_J_m3, fusion_two_component_ledger_v1* ledger) {
    if (ledger != nullptr) {
        *ledger = {};
    }

    if (cells < 1 || baths < 0) {
        return PB11_STATUS_INVALID_ARGUMENT;
    }

    const std::size_t n = static_cast<std::size_t>(cells);
    const std::size_t b = static_cast<std::size_t>(baths);
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(double) ||
        n + 1 < n ||
        b > std::numeric_limits<std::size_t>::max() / sizeof(double) ||
        (cells > 1 &&
         b > std::numeric_limits<std::size_t>::max() /
                     static_cast<std::size_t>(cells - 1))) {
        return PB11_STATUS_INVALID_ARGUMENT;
    }

    /* Do this before every validation return, matching the existing trial
     * kernel's valid-dimension output contract. */
    clear_outputs(cells, baths, trial_s_m3, trial_t_m3, heat_to_bath_J_m3,
                  ledger);

    if (trial_s_m3 == nullptr || trial_t_m3 == nullptr || ledger == nullptr ||
        (baths > 0 && heat_to_bath_J_m3 == nullptr)) {
        return PB11_STATUS_NULL_OUTPUT;
    }

    if (edges_J == nullptr || old_s_m3 == nullptr || old_t_m3 == nullptr ||
        birth_s_m3_s == nullptr || birth_t_m3_s == nullptr ||
        escape_s_inv == nullptr || transfer_s_inv == nullptr ||
        (baths > 0 && bath_kT_J == nullptr) ||
        (baths > 0 && cells > 1 && diffusion_J2_s == nullptr)) {
        return PB11_STATUS_INVALID_ARGUMENT;
    }

    if (!std::isfinite(dt_s)) {
        return PB11_STATUS_INVALID_ARGUMENT;
    }
    if (dt_s <= 0.0) {
        return PB11_STATUS_OUT_OF_RANGE;
    }

    for (std::size_t i = 0; i <= n; ++i) {
        if (!std::isfinite(edges_J[i]) || edges_J[i] < 0.0) {
            return PB11_STATUS_INVALID_ARGUMENT;
        }
        if (i > 0 && edges_J[i] <= edges_J[i - 1]) {
            return PB11_STATUS_OUT_OF_RANGE;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!finite_nonnegative(old_s_m3[i]) ||
            !finite_nonnegative(old_t_m3[i]) ||
            !finite_nonnegative(birth_s_m3_s[i]) ||
            !finite_nonnegative(birth_t_m3_s[i]) ||
            !finite_nonnegative(escape_s_inv[i]) ||
            !finite_nonnegative(transfer_s_inv[i])) {
            return PB11_STATUS_INVALID_ARGUMENT;
        }
    }
    for (std::size_t j = 0; j < b; ++j) {
        if (!finite_positive(bath_kT_J[j])) {
            if (!std::isfinite(bath_kT_J[j])) {
                return PB11_STATUS_INVALID_ARGUMENT;
            }
            return PB11_STATUS_OUT_OF_RANGE;
        }
    }
    const std::size_t diffusion_count =
        (cells > 1) ? b * static_cast<std::size_t>(cells - 1) : 0;
    for (std::size_t k = 0; k < diffusion_count; ++k) {
        if (!finite_nonnegative(diffusion_J2_s[k])) {
            return PB11_STATUS_INVALID_ARGUMENT;
        }
    }

    try {
        std::vector<double> effective_escape(n);
        std::vector<double> effective_birth_t(n);
        std::vector<double> s_trial(n);
        std::vector<double> t_trial(n);
        std::vector<double> s_heat(b);
        std::vector<double> t_heat(b);
        fusion_kinetic_ledger_v1 s_ledger{};
        fusion_kinetic_ledger_v1 t_ledger{};

        for (std::size_t i = 0; i < n; ++i) {
            const Real sum = static_cast<Real>(escape_s_inv[i]) +
                             static_cast<Real>(transfer_s_inv[i]);
            if (!put_double(sum, effective_escape[i])) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
        }

        const int s_status = fusion_c_energy_fp_trial(
            cells, baths, dt_s, edges_J, old_s_m3, bath_kT_J, diffusion_J2_s,
            birth_s_m3_s, effective_escape.data(), 0.0, s_trial.data(),
            baths > 0 ? s_heat.data() : nullptr, &s_ledger);
        if (s_status != PB11_STATUS_OK) {
            return s_status;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!finite_nonnegative(s_trial[i])) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
        }
        for (std::size_t j = 0; j < b; ++j) {
            if (!std::isfinite(s_heat[j])) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            const Real transferred_birth =
                static_cast<Real>(transfer_s_inv[i]) *
                static_cast<Real>(s_trial[i]);
            const Real sum = static_cast<Real>(birth_t_m3_s[i]) +
                             transferred_birth;
            if (!std::isfinite(sum) || sum < 0.0L ||
                sum > std::numeric_limits<double>::max()) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
            // A positive product in a negligible tail may lie below double
            // range even when the global transfer is well resolved. Round
            // this cell source as the FP kernel rounds positive tail states;
            // the recomputed TOTAL N/U residuals below still bound the loss.
            // This is IEEE underflow, not a population floor or negative clip.
            effective_birth_t[i] = static_cast<double>(sum);
        }

        const int t_status = fusion_c_energy_fp_trial(
            cells, baths, dt_s, edges_J, old_t_m3, bath_kT_J, diffusion_J2_s,
            effective_birth_t.data(), escape_s_inv, 0.0, t_trial.data(),
            baths > 0 ? t_heat.data() : nullptr, &t_ledger);
        if (t_status != PB11_STATUS_OK) {
            return t_status;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!finite_nonnegative(t_trial[i])) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
        }
        for (std::size_t j = 0; j < b; ++j) {
            if (!std::isfinite(t_heat[j])) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
        }

        fusion_two_component_ledger_v1 result{};
        Real initial_number = 0.0L;
        Real final_number = 0.0L;
        Real initial_energy = 0.0L;
        Real final_energy = 0.0L;
        Real born_number = 0.0L;
        Real born_energy = 0.0L;
        Real escaped_number = 0.0L;
        Real escaped_energy = 0.0L;
        Real transferred_number = 0.0L;
        Real transferred_energy = 0.0L;

        std::vector<double> heat_total(b);
        Real net_heat = 0.0L;
        for (std::size_t j = 0; j < b; ++j) {
            const Real sum = static_cast<Real>(s_heat[j]) +
                             static_cast<Real>(t_heat[j]);
            if (!put_double(sum, heat_total[j])) {
                return PB11_STATUS_NUMERICAL_FAILURE;
            }
            net_heat += static_cast<Real>(heat_total[j]);
        }

        for (std::size_t i = 0; i < n; ++i) {
            const Real lower = static_cast<Real>(edges_J[i]);
            const Real upper = static_cast<Real>(edges_J[i + 1]);
            const Real energy = (lower + upper) / 2.0L;
            const Real old_s = static_cast<Real>(old_s_m3[i]);
            const Real old_t = static_cast<Real>(old_t_m3[i]);
            const Real trial_s = static_cast<Real>(s_trial[i]);
            const Real trial_t = static_cast<Real>(t_trial[i]);
            const Real birth = static_cast<Real>(dt_s) *
                               (static_cast<Real>(birth_s_m3_s[i]) +
                                static_cast<Real>(birth_t_m3_s[i]));
            const Real escape = static_cast<Real>(dt_s) *
                                static_cast<Real>(escape_s_inv[i]) *
                                (trial_s + trial_t);
            const Real transfer = static_cast<Real>(dt_s) *
                                  static_cast<Real>(transfer_s_inv[i]) *
                                  trial_s;

            initial_number += old_s + old_t;
            final_number += trial_s + trial_t;
            initial_energy += energy * (old_s + old_t);
            final_energy += energy * (trial_s + trial_t);
            born_number += birth;
            born_energy += energy * birth;
            escaped_number += escape;
            escaped_energy += energy * escape;
            transferred_number += transfer;
            transferred_energy += energy * transfer;
        }

        const Real particle_error =
            final_number - initial_number - born_number + escaped_number;
        const Real energy_error = final_energy - initial_energy - born_energy +
                                  escaped_energy + net_heat;
        const Real particle_scale = initial_number + final_number +
                                    born_number + escaped_number;
        Real heat_scale = 0.0L;
        for (std::size_t j = 0; j < b; ++j) {
            heat_scale += std::abs(static_cast<Real>(heat_total[j]));
        }
        const Real energy_scale = initial_energy + final_energy + born_energy +
                                  escaped_energy + heat_scale;
        if (!std::isfinite(particle_error) ||
            !std::isfinite(energy_error) ||
            !std::isfinite(particle_scale) ||
            !std::isfinite(energy_scale) ||
            std::abs(particle_error) > 1.0e-10L * particle_scale ||
            std::abs(energy_error) > 1.0e-10L * energy_scale) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        if (!put_double(initial_number, result.total.initial_number_m3) ||
            !put_double(final_number, result.total.final_number_m3) ||
            !put_double(initial_energy, result.total.initial_energy_J_m3) ||
            !put_double(final_energy, result.total.final_energy_J_m3) ||
            !put_double(born_number, result.total.born_number_m3) ||
            !put_double(born_energy, result.total.born_energy_J_m3) ||
            !put_double(escaped_number, result.total.escaped_number_m3) ||
            !put_double(escaped_energy, result.total.escaped_energy_J_m3) ||
            !put_double(particle_error,
                        result.total.particle_balance_error_m3) ||
            !put_double(energy_error,
                        result.total.energy_balance_error_J_m3) ||
            !put_double(transferred_number, result.transferred_number_m3) ||
            !put_double(transferred_energy, result.transferred_energy_J_m3)) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }
        /* The decomposition has no host thermalization operator. */
        result.total.thermalized_number_m3 = 0.0;
        result.total.thermalized_energy_J_m3 = 0.0;

        std::copy(s_trial.begin(), s_trial.end(), trial_s_m3);
        std::copy(t_trial.begin(), t_trial.end(), trial_t_m3);
        if (b > 0) {
            std::copy(heat_total.begin(), heat_total.end(),
                      heat_to_bath_J_m3);
        }
        *ledger = result;
        return PB11_STATUS_OK;
    } catch (...) {
        clear_outputs(cells, baths, trial_s_m3, trial_t_m3,
                      heat_to_bath_J_m3, ledger);
        return PB11_STATUS_EXCEPTION;
    }
}
