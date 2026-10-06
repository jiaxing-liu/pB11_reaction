#ifndef FUSION_BEAM_BIRTH_H
#define FUSION_BEAM_BIRTH_H
#include "fusion_thermal_birth.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_beam_birth_options_v1 {
 double relative_max_J, angular_max_exponent, ground_state_q_J, cutoff_J;
 double l1_fraction, relative_phase, narrow_peak_fraction, continuum_peak_scale;
 int continuation, pb_low, remainder_policy, broad_mode, fsci_policy;
 int relative_order, angular_order, nq, ncos;
} fusion_beam_birth_options_v1;
typedef struct fusion_beam_birth_v1 {
 fusion_thermal_birth_v1 spectrum;
 double relative_retained_probability;
 double retained_pair_probability;
 double angular_omitted_pair_probability;
} fusion_beam_birth_v1;
/* Explicit per-projectile support builder for coupled/cache callers.
 * Copy base unchanged except relative_max_J=max(base.relative_max_J,E[J]).
 * Validates only finite E>=0 and finite base.relative_max_J>0; remaining
 * source controls and retained-domain accuracy are checked by grid/table APIs.
 * It does not extend the pB event domain or certify any discrepancy gate.
 * BASE==OUT is allowed. Non-null OUT clears on failure; no hidden state.
 * Existing grid/table APIs continue to honor their explicit source cutoff.
 */
int fusion_c_beam_birth_resolve_support(double projectile_energy_J,
 const fusion_beam_birth_options_v1 *base,fusion_beam_birth_options_v1 *out);
/* Monoenergetic projectile with E[J]>=0 against a stationary Maxwellian
 * target kT[J]>=0; projectile_slot0/1 refers to canonical channel reactants.
 * Nuclear models, isotropic outgoing CM-event closure, product mapping,
 * arrays and units follow thermal_birth. Reactant debit arrays always use
 * CANONICAL order, including when projectile_slot=1. No density/self-pair
 * factor, burn/state update, heat deposition, or beam pitch is supplied.
 *
 * Relative energy is integrated over[0,relative_max_J], with Gaussian and
 * nuclear knots. At fixed relative speed w, target thermal speed u and beam
 * speed v, integrate c=cos(w,v) with density proportional to exp(kappa*c),
 * kappa=2*v*w/u^2. For kappa>=1 use z=kappa*(1-c), truncated explicitly at
 * angular_max_exponent8..80 (or2*kappa if smaller). For kappa<1 integrate
 * the full cosine interval. Orders: relative4..64, angular4..32; remaining
 * event controls are as in thermal_birth. Cutoffs are not renormalized.
 *
 * Probabilities describe the retained relative-energy domain and angular
 * truncation; numerical quadrature sums may slightly exceed one and are not
 * clipped. Omitted probability is NOT a missing-reaction bound. The three
 * spectrum.cm_* fields are zero/not applicable. Reference rate/debits use
 * the full matching beam cross-section continuation model. Discrepancies
 * expose truncation/quadrature errors; OK is finite/conservative evaluation,
 * not accuracy certification. The existing beam speed-ratio domain applies.
 * Target kT=0 uses the exact stationary target with the same relative cutoff.
 * pB parent energy must stay within the existing12MeV event domain.
 * Non-null outputs clear on failure; birth clears for valid cell counts.
 * No hidden mutable state; arrays must not overlap.
 */
int fusion_c_beam_birth_grid(int channel,int projectile_slot,
 double projectile_energy_J,double target_kT_J,
 const fusion_beam_birth_options_v1 *options,int cells,const double *edges_J,
 double *birth,fusion_beam_birth_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
