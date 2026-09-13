#ifndef FUSION_TWO_COMPONENT_H
#define FUSION_TWO_COMPONENT_H
#include "fusion_coulomb.h"
#include "fusion_kinetics.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fusion_two_component_ledger_v1 {
    fusion_kinetic_ledger_v1 total;
    double transferred_number_m3;
    double transferred_energy_J_m3;
} fusion_two_component_ledger_v1;

/* Localized finite-temperature ion-friction-divergence rate, s^-1.
 * lambda_i(v)=4 nu_i/(sqrt(pi)*vti^3)*exp(-v^2/vti^2),
 * nu_i=4 pi ni Za^2 <Zi^2> (e^2/(4 pi eps0))^2 lnLambda/(ma mi).
 * Caller selects an ION bath explicitly. Same dilute/classical/isotropic
 * assumptions as fusion_c_coulomb_energy. This defines an internal kinetic
 * component transfer, NOT physical escape or a fluid thermalization rate.
 * energy>=0, mass>0, bath temperature>0, log>0. Zero charge/density -> zero.
 * A positive exponential tail below double range rounds to zero; overflow
 * rejects. Output is zero on failure, no exception crosses the ABI.
 */
int fusion_c_coulomb_transfer_rate(double energy_J, double mass_kg,
    double charge_number, const fusion_maxwellian_bath_v1 *ion_bath,
    double *rate_s_inv);

/* Shared-grid exact decomposition of the SAME frozen linear FP operator L:
 * (I-dt L+dt E+dt Lambda) S_new = S_old+dt birth_S
 * (I-dt L+dt E) T_new = T_old+dt birth_T+dt Lambda S_new.
 * Thus S_new+T_new satisfies the original backward Euler FP equation with
 * birth_S+birth_T and physical escape E, up to floating-point roundoff.
 * Both components are kinetic distributions over the WHOLE grid. T is a
 * thermal-SCALE component, not a Maxwellian or the host's thermal ash.
 * Lambda>=0 is explicit, cellwise s^-1, e.g. sum of selected ion rates above.
 * Units/layout/domain match fusion_c_energy_fp_trial; identical bath/escape
 * coefficients act on both components. There is no first-cell removal.
 * Ledger total contains ONLY external births/physical escape; its thermalized
 * fields are exactly zero. Transfer N and carried energy are internal; do not
 * add them to a bath or host source. Per-bath heat is S+T collision heat.
 * Inputs/output arrays must not overlap. Explicit old/trial state, no hidden
 * counters. Positive per-cell transfer-source tails may round below double
 * range; the combined N/U residual check still bounds this rounding loss.
 * For valid dimensions nonnull output arrays and ledger zero on
 * any failure. Bath arrays may be NULL at baths=0, diffusion at cells=1.
 */
int fusion_c_two_component_trial(int cells, int baths, double dt_s,
    const double *edges_J, const double *old_s_m3, const double *old_t_m3,
    const double *bath_kT_J, const double *diffusion_J2_s,
    const double *birth_s_m3_s, const double *birth_t_m3_s,
    const double *escape_s_inv, const double *transfer_s_inv,
    double *trial_s_m3, double *trial_t_m3, double *heat_to_bath_J_m3,
    fusion_two_component_ledger_v1 *ledger);
#ifdef __cplusplus
}
#endif
#endif
