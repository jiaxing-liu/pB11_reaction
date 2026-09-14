#ifndef FUSION_TARGET_NETWORK_H
#define FUSION_TARGET_NETWORK_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
/* A reaction edge consumes one particle from fast_index and one from target_index.
 * Parallel edges represent competing channels and share BOTH inventories.
 * K is m^3/s; M is reaction-conditioned target energy integral J m^3/s.
 * K=0 requires M=0. No identical-pair half factor for distinct pools. */
typedef struct fusion_target_network_edge_v1 {
 int fast_index, target_index;
 double reactivity_m3_s, target_energy_reactivity_J_m3_s;
} fusion_target_network_edge_v1;
typedef struct fusion_target_network_v1 {
 double reactions_m3, removed_fast_energy_J_m3, removed_target_energy_J_m3;
 double max_fast_number_relative_residual, max_target_number_relative_residual;
 double max_target_energy_relative_residual;
 int iterations;
} fusion_target_network_v1;
/* Simultaneous frozen-coefficient backward-Euler depletion:
 * F_i'=F_i/(1+dt*sum_e(i) K_e B_target(e)'),
 * B_j'=B_j-sum_e(j) dt*K_e*F_fast(e)'*B_j'.
 * Per-edge losses remain explicit, preserving weak events when F-F' rounds away.
 * Fast representative energies are supplied per bin [J]; target U[J/m^3].
 * Target debit=sum loss*M/K. Any energy exhaustion fails, without clipping.
 * nfast1..600000, ntarget1..6, nedge0..2000000; zero edges pointer may be null.
 * No products/Q/heat/FP/state mutation. Caller must compose births using these
 * SAME coefficients and accept complete fuel/kinetic/energy state atomically.
 * All input and output arrays disjoint. On error nonnull outputs clear for
 * valid dimensions; no exceptions cross C. dt positive finite. */
int fusion_c_target_network_trial(int nfast,int ntarget,int nedge,double dt_s,
 const double *fast_energy_J,const double *old_fast_m3,
 const double *target_number_m3,const double *target_energy_J_m3,
 const fusion_target_network_edge_v1 *edges,double *trial_fast_m3,
 double *trial_target_m3,double *trial_target_energy_J_m3,
 double *edge_reactions_m3,fusion_target_network_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
