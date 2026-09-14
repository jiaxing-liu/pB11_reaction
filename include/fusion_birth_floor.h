#ifndef FUSION_BIRTH_FLOOR_H
#define FUSION_BIRTH_FLOOR_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit numerical projection controls; no defaults or physical heat model. */
typedef struct fusion_birth_floor_options_v1 {
 double first_center_J, ion_kT_J, ion_energy_J_m3;
 double max_center_over_ion_kT, max_ion_energy_fraction;
} fusion_birth_floor_options_v1;
typedef struct fusion_birth_floor_ledger_v1 {
 double born_number_m3[6], born_energy_J_m3[6];
 double mapped_number_m3[6], mapped_energy_J_m3[6];
 double ion_energy_correction_J_m3[6], energy_residual_J_m3[6];
 double remaining_ion_energy_J_m3;
} fusion_birth_floor_ledger_v1;
/* Pure, additive candidate utility, NOT installed in coupled trial by default.
 * Inputs are per-trial amounts, canonical charged species order. For each
 * N>=0 and 0<=U<=N*first_center, map ALL N to the first center. The returned
 * signed ion correction U-Umap must be applied exactly once, atomically with
 * the kinetic source. It is numerical remapping energy, NOT Coulomb heat,
 * thermalization or nuclear Q. No particles enter the fluid population.
 * Caller keeps physical birth U distinct from mapped kinetic U in its ledger.
 * Explicit gates require center/kTi<=max_center_over_ion_kT (0..1) and total
 * borrowed energy/Uion<=max_ion_energy_fraction (0..1), and positive residual
 * ion energy. Empty input returns unchanged ion energy; geometry is validated.
 * This bounds an energy displacement, not distribution L1 or physical error.
 * The approximation must converge with grid refinement in the coupled host.
 * No handling of upper spill, no hidden state, no automatic reservoir choice.
 * Invalid or overflowing outputs reject with all output bytes zero. Positive
 * energy amounts can underflow to zero; signed correction uses the RETURNED
 * mapped energy, with residual checked against only measured IEEE rounding.
 * Particle numbers are copied exactly, including representable subnormals.
 * Inputs/outputs must not overlap. Existing coupled source ABI is unchanged.
 */
int fusion_c_birth_floor_project(const fusion_birth_floor_options_v1 *options,
 const double born_number_m3[6], const double born_energy_J_m3[6],
 fusion_birth_floor_ledger_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
