#ifndef FUSION_TARGET_BURN_H
#define FUSION_TARGET_BURN_H
#include "fusion_beam.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_target_burn_v1 {
 double initial_fast_number_m3, final_fast_number_m3;
 double initial_fast_energy_J_m3, final_fast_energy_J_m3;
 double initial_target_number_m3, final_target_number_m3;
 double initial_target_energy_J_m3, final_target_energy_J_m3;
 double reactions_m3, removed_fast_energy_J_m3, removed_target_energy_J_m3;
 double fast_number_residual_m3, target_number_residual_m3, energy_residual_J_m3;
} fusion_target_burn_v1;
/* Backward-Euler reaction-only trial for distinct fast and thermal pools.
 * Each event consumes one fast particle and one shared target. At fixed
 * supplied coefficients K_i (m^3/s): N_i'=N_i/(1+dt*K_i*n_target'), with
 * n_target'=n_target-sum(N_i-N_i'). No identical-pair 1/2: pools are distinct.
 * edges[cells+1] J, arithmetic centers; old/trial/loss[cells] m^-3.
 * target_energy_reactivity[cells] is integral sigma*v*E_target (J m^3/s),
 * NOT 3kT/2 times K unless independently justified. Removed target energy
 * uses this reaction-conditioned moment. Coefficients remain frozen during
 * this first-order trial. Negative target energy rejects the trial with
 * NUMERICAL_FAILURE; caller reduces dt or iterates the background model.
 * No product deposition, Q release, collisions, fast-fast reactions or ash
 * source is performed here. Reactions/removed energy are handed to the birth
 * source owner once after acceptance, never immediately deposited as heat.
 * All inputs finite/nonnegative, dt>0, cells>=1; K=0 requires moment=0.
 * Inputs unchanged; arrays must not overlap. Explicit trial has no hidden
 * state/counters: discard means rollback, caller copies outputs to accept.
 * Non-null outputs clear on error for valid dimensions.
 */
int fusion_c_target_burn_trial(int cells,double dt_s,const double *edges_J,
 const double *old_fast_m3,double target_number_m3,double target_energy_J_m3,
 const double *reactivity_m3_s,const double *target_energy_reactivity_J_m3_s,
 double *trial_fast_m3,double *reaction_loss_m3,fusion_target_burn_v1 *out);
/* Physical coefficient composition: evaluate the documented finite nuclear
 * window against a Maxwellian target, then call the trial above. target_kT
 * in J>=0; initial target energy is 3*n*kT/2. Masses are explicit SI inputs
 * and must match the channel. All window diagnostics[cells] remain visible;
 * OK means a WINDOW burn trial, not a complete physical burn prediction.
 * Unknown-window probabilities are NOT missing-reaction error bounds.
 * This wrapper does not select species slots or write to a transport host.
 */
int fusion_c_beam_target_burn_window_trial(int channel,int cells,double dt_s,
 double projectile_mass_kg,double target_mass_kg,const double *edges_J,
 const double *old_fast_m3,double target_number_m3,double target_kT_J,
 double *trial_fast_m3,double *reaction_loss_m3,fusion_beam_window_v1 *windows,
 fusion_target_burn_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
