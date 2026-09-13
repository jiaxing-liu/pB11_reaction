#ifndef FUSION_HANDOFF_H
#define FUSION_HANDOFF_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fusion_maxwellian_grid_v1 {
    double below_probability, above_probability;
    double represented_mean_energy_J;
    double probability_balance_error;
} fusion_maxwellian_grid_v1;

typedef struct fusion_handoff_ledger_v1 {
    double initial_number_m3, initial_energy_J_m3;
    double remaining_number_m3, remaining_energy_J_m3;
    double fluid_number_m3, fluid_energy_J_m3;
    double bath_energy_correction_J_m3;
    double distribution_L1, relative_mean_energy_error;
    double outside_grid_probability, represented_Maxwellian_mean_energy_J;
    double particle_balance_error_m3, energy_balance_error_J_m3;
} fusion_handoff_ledger_v1;

/* Normalized nonrelativistic 3D Maxwellian energy probability integrated
 * over each [edge_i,edge_i+1]. Uses Gamma(3/2), not midpoint sampling.
 * Cells>=1, finite increasing nonnegative edges[J], kT>0[J]. q[cells]
 * is probability per FULL-distribution particle (not renormalized in-grid).
 * Below/above tails complete normalization. represented_mean=sum(q*center)
 * uses arithmetic cell centers; it generally differs from 1.5*kT and omits
 * outside-grid energy. Positive probability underflow may round to zero.
 * Output arrays/ledger zero on failure for valid dimensions. No I/O/state.
 */
int fusion_c_maxwellian_energy_grid(int cells, double kT_J,
    const double *edges_J, double *probability,
    fusion_maxwellian_grid_v1 *out);

/* Controlled WHOLE-candidate kinetic->fluid projection, not a time advance.
 * Candidate old[cells] is cell-integrated number in m^-3. Measure complete
 * BIN-probability L1=sum|old_i/N-q_i|+q_outside and independently
 * abs(U/(1.5*kT*N)-1). Both must be <= explicit tolerances before projection.
 * 0<=max_L1<=2, 0<=max_relative_mean_error<=1. Empty candidates do not project.
 * On OK with projected=0, trial=old and fluid/correction sources are zero.
 * On OK with projected=1, trial=0, fluid N=Nold, fluid U=1.5*kT*Nold;
 * correction=Uold-fluidU (signed energy TO the target thermal bath).
 * Returned source quantities are amounts per trial, NOT rates. Caller must
 * account for fluid energy plus correction exactly once, validate the bath
 * can accept the signed correction, and commit/discard atomically.
 * This function does not evolve the bath or assert a Maxwellian stays valid
 * under later unequal/evolving backgrounds, collisions or transport. It
 * measures grid-resolved error, not unknown within-cell distribution shape.
 * Changing the target kT requires a new trial from the accepted old state.
 * Inputs and outputs must not overlap. No clipping/hidden counters. On any
 * failure, nonnull outputs are zero for valid dimensions; no exception crosses
 * C. Numerical overflow or invalid probability normalization rejects.
 */
int fusion_c_maxwellian_handoff_trial(int cells, double target_kT_J,
    double max_L1, double max_relative_mean_error,
    const double *edges_J, const double *old_number_m3,
    double *trial_number_m3, int *projected, fusion_handoff_ledger_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
