#ifndef FUSION_KINETIC_GEOMETRY_H
#define FUSION_KINETIC_GEOMETRY_H
#include "fusion_source_state.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Stateless first-order composition: radial FV transport from Va to Vb, then
 * energy drift/work at Vb, for both six-species kinetic components S and T.
 * Existing radial/work discretizations and their model limitations apply.
 * No nuclear source, thermal feedback, state commit or host dependency here.
 *
 * Densities: [zone][species][energy cell] (Fortran cells,6,zones), m^-3.
 * edges[cells+1] J; Va/Vb[zones] m3; dt>0 s; compression[zones] s^-1.
 * advection[zones+1] outward-positive m3/s; pure mesh motion supplies -G.
 * conductance[species][cell][face] >=0 m3/s, common to S/T at the same
 * species/energy; caller supplies any physical diffusion explicitly.
 * Each boundary array [species][cell][2] contains inner/outer densities.
 * Vacuum, reflecting, or imposed influx must be expressed explicitly through
 * advection/conductance/boundary inputs, never inferred from scalar moments.
 *
 * Returns physical Vb-normalized S/T and one EXTENSIVE transport/work/domain
 * ledger per zone, directly usable by stage_volume. Shared interior face
 * amounts have opposite signs; energy-work density ledgers are multiplied by
 * Vb. Lower/upper energy-domain outflow is neither orbit loss nor thermal ash.
 * Source evolution, if composed before this call at Va, must retain its own
 * Va-normalized source ledger for the later state/host source stage.
 *
 * zones 1..300, cells 1..100000; 12*cells*(2*zones+1)<=50000000.
 * Inputs/outputs must not overlap. For admissible dimensions all supplied
 * outputs clear on any failure; no partial zone or component is published.
 * Invalid dimensions do not dereference outputs. Finite trace rounding follows
 * the underlying kernels; no hidden population floor or empirical loss.
 */
int fusion_c_kinetic_geometry_trial(int zones,int cells,double dt_s,
 const double*edges_J,const double*volume_old_m3,const double*volume_new_m3,
 const double*advection_m3_s,const double*compression_s_inv,
 const double*conductance_m3_s,const double*boundary_s_m3,const double*boundary_t_m3,
 const double*old_s_m3,const double*old_t_m3,double*trial_s_m3,double*trial_t_m3,
 fusion_transport_ledger_v1*transport_step);
#ifdef __cplusplus
}
#endif
#endif
