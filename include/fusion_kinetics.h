#ifndef FUSION_KINETICS_H
#define FUSION_KINETICS_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fusion_kinetic_ledger_v1 {
    double initial_number_m3, final_number_m3;
    double initial_energy_J_m3, final_energy_J_m3;
    double born_number_m3, born_energy_J_m3;
    double escaped_number_m3, escaped_energy_J_m3;
    double thermalized_number_m3, thermalized_energy_J_m3;
    double particle_balance_error_m3, energy_balance_error_J_m3;
} fusion_kinetic_ledger_v1;

/* One implicit finite-volume trial step for an isotropic ENERGY distribution.
 * This numerical kernel accepts externally computed collision coefficients;
 * it does not itself assert that arbitrary supplied coefficients are Coulomb.
 *
 * cells >=1, baths>=0. Edges[cells+1] are increasing, nonnegative J;
 * cell energy is its arithmetic midpoint. old_number[cells] and returned
 * trial_number[cells] are integrated cell populations in m^-3, NOT density/J.
 * bath_kT[baths] are positive J. diffusion[baths*(cells-1)] stores D_E at
 * interior faces in J^2/s, bath-major order. Birth[cells] is m^-3 s^-1;
 * escape[cells] is s^-1. thermalization_s_inv removes first-cell particles
 * only, with their carried midpoint energy. It is an explicit caller policy.
 *
 * Each bath's flux is -D_E [dp/dE+(1/kT-1/(2E))*p], p=dN/dE.
 * Exponential fitting preserves the sampled Maxwellian sqrt(E)*exp(-E/kT).
 * Physical grid boundaries otherwise reflect. heat_to_bath[baths] is J/m^3,
 * positive for energy deposited in the bath, negative for fast-ion heating.
 * Thermalization/escape carry separate particle AND residual-energy ledgers.
 *
 * Inputs remain unchanged; all input/output arrays must be nonoverlapping.
 * Repeated calls from the same old state are identical; rejection means
 * discard outputs, acceptance means caller copies the trial state. No hidden
 * state, counters or file I/O. A host must preserve its accepted grid/state.
 * Bath arrays may be NULL when baths=0; diffusion may be NULL when cells=1.
 * For valid dimensions, non-NULL output arrays/ledger are zeroed on failure.
 * All inputs must be finite/nonnegative, dt>0, kT>0. Returns PB11_STATUS_*.
 */
int fusion_c_energy_fp_trial(int cells, int baths, double dt_s,
    const double *edges_J, const double *old_number_m3,
    const double *bath_kT_J, const double *diffusion_J2_s,
    const double *birth_m3_s, const double *escape_s_inv,
    double thermalization_s_inv, double *trial_number_m3,
    double *heat_to_bath_J_m3, fusion_kinetic_ledger_v1 *ledger);

#ifdef __cplusplus
}
#endif
#endif
