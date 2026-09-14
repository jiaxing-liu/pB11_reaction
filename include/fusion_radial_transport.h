#ifndef FUSION_RADIAL_TRANSPORT_H
#define FUSION_RADIAL_TRANSPORT_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_radial_ledger_v1 {
 double initial_number,final_number,inner_inward_number,outer_outward_number;
 double balance_error;
} fusion_radial_ledger_v1;
/* Stateless implicit upwind advection/diffusion of independent number-density
 * components across radial cells. Components can be species/energy/S/T bins;
 * this operator does not change energy-bin identity or compute adiabatic work.
 * volumes_old/new[zones]: positive m3. dt: positive seconds. Every face has a
 * signed ADVECTIVE VOLUME rate a [m3/s], positive from lower to higher cell.
 * Magnetic grid motion alone supplies a=-swept-volume-rate, not +rate.
 * conductance[components*(zones+1)]: nonnegative diffusive m3/s, explicitly
 * provided by caller (no inferred diffusivity or orbit-loss model).
 * old/trial densities[components*zones]: m^-3; boundary[components*2]: m^-3.
 * Arrays are component-major. Face0 is inner, face zones outer. For left/right
 * densities, outward face flux = (max(a,0)+k)*n_left+(min(a,0)-k)*n_right.
 * Boundary density is used for both prescribed influx and diffusive exchange;
 * reflecting boundaries require a=k=0. Volume update plus mesh flux must obey
 * caller's discrete geometry law to preserve a uniform physical density.
 * Returns integrated signed face amounts[components*(zones+1)] and one ledger
 * per component: Nnew-Nold=inner_inward-outer_outward. Negative boundary amounts
 * are allowed. Old/new density normalization uses the supplied actual volumes.
 *
 * Finite nonnegative densities, positive volumes, finite signed rates, no array
 * overlap. zones/components>=1 and total state+face entries<=50000000. Invalid
 * dimensions do not dereference outputs; otherwise outputs clear on failure.
 * IEEE underflow may round finite trace values to zero; extended-precision
 * zone/global balance checks include measured conversion errors, not floors.
 * No host state, grid indexing convention, hidden advancement or file I/O.
 */
int fusion_c_radial_transport_trial(int zones,int components,double dt_s,
 const double*volume_old_m3,const double*volume_new_m3,
 const double*advection_volume_m3_s,const double*conductance_m3_s,
 const double*old_density_m3,const double*boundary_density_m3,
 double*trial_density_m3,double*face_amount,fusion_radial_ledger_v1*ledger);
#ifdef __cplusplus
}
#endif
#endif
