#ifndef FUSION_THERMAL_BIRTH_H
#define FUSION_THERMAL_BIRTH_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit closure of the effective pB population remainders. None is a
 * calibrated ghost/continuum decomposition:0 uses l2 for narrow remainder
 * and broad_mode for other;1 uses l2 for both;2 broad_mode for both.
 * Alpha0 peak fraction always comes from fusion_pb_population. */
enum fusion_pb_remainder_source_policy {
 FUSION_PB_REMAINDER_ENTRANCE_PROXY=0,
 FUSION_PB_REMAINDER_ALL_LOW=1,
 FUSION_PB_REMAINDER_ALL_BROAD=2
};
typedef struct fusion_thermal_birth_options_v1 {
 double relative_max_J, cm_max_kT, ground_state_q_J, cutoff_J;
 double l1_fraction, relative_phase, narrow_peak_fraction, continuum_peak_scale;
 int continuation, pb_low, remainder_policy, broad_mode, fsci_policy;
 int relative_order, cm_order, nq, ncos;
} fusion_thermal_birth_options_v1;
typedef struct fusion_thermal_birth_v1 {
 double reactivity_m3_s, reference_reactivity_m3_s;
 double reactant_energy_moment_J_m3_s[2];
 double reference_reactant_energy_moment_J_m3_s[2];
 double product_energy_moment_J_m3_s;
 double number_residual_m3_s, energy_residual_J_m3_s;
 double relative_rate_discrepancy, relative_reactant_energy_discrepancy;
 double cm_retained_probability, cm_tail_probability, cm_tail_energy_moment_J;
 double max_cm_energy_shift_fraction, max_shell_remap_fraction;
 double below_number_m3_s[7], below_energy_J_m3_s[7];
 double above_number_m3_s[7], above_energy_J_m3_s[7];
} fusion_thermal_birth_v1;
/* Product energy-grid source COEFFICIENTS for two isotropic, zero-drift
 * Maxwellians at the SAME kT[J]>0. Canonical masses and Q; channels0..4.
 * cells must be 1..100000. birth[7*cells], species-major, IDs0..5 plus neutron6; edges[cells+1]J,
 * arithmetic centers. Units number*m^3/s, not density/time: multiply by
 * n_a*n_b/(1+delta_ab) exactly once. Neutrons are explicit; no automatic
 * thermalization, escape or treatment of grid spill.
 *
 * Integrates reacting relative energy AND the full radial CM Maxwellian
 * within explicit finite cutoffs; NO evaluation at mean reacting energy.
 * The outgoing CM event orientation is uniform over rotations (isotropic
 * unpolarized closure). This API supplies scalar-energy marginals, NOT
 * pitch/current, beam or unequal-temperature angular correlations.
 * Classical reactant budgets match the existing NR cross sections. Products
 * are individually on shell: alpha1 NR amplitude events retain weights at
 * A0=Q+Erel, then a COMMON momentum rescaling enforces actual parent A and
 * zero CM momentum before exact boosts. This is an explicit kinematic
 * mapping of an NR amplitude model, not a relativistic nuclear amplitude.
 * Alpha0 uses exact narrow-width sequential kinematics, not equal energies.
 *
 * Relative Gauss order4..64 per segmented energy interval; CM order4..32
 * per speed interval; nq/ncos4..1024; cm_max_kT8..80. relative_max_J>0;
 * for pB the largest parent A must be<=12MeV. All controls validated even
 * when inactive. Ground-state q>=0; numerical cutoff1..10keV; k/peak0..1;
 * continuum scale0..2; broad mode1/3/13; policies as above/rate-model APIs.
 *
 * Finite tails are NOT renormalized. Returned rate/debits correspond to
 * the actual quadrature measure; reference values integrate the full
 * declared cross-section model. Signed rate and max absolute debit relative
 * discrepancies combine finite support and numerical error, NOT certified
 * nuclear uncertainty bounds. Convergence and tail checks are REQUIRED
 * before using these coefficients to advance a source. cm_retained_probability is the actual quadrature sum,
 * not the analytic retained probability; rounding may put it just above 1.
 * Success only means
 * finite, conservative evaluation. All non-NULL outputs clear on failure;
 * birth clears only for a valid cell count (its extent is otherwise unknown);
 * no global mutable state or hidden sampling. Arrays must not overlap. */
int fusion_c_thermal_birth_grid(int channel,double kT_J,
 const fusion_thermal_birth_options_v1 *options,int cells,const double *edges_J,
 double *birth,fusion_thermal_birth_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
