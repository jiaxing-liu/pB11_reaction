#ifndef FUSION_FLUX_VOLUME_DIAGNOSTIC_H
#define FUSION_FLUX_VOLUME_DIAGNOSTIC_H
#include "fusion_flux_core.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_flux_volume_diagnostic_v1 {
 int version, angular_nodes, boundary_orientation, reserved;
 double signed_volume_m3, absolute_volume_m3;
 double min_abs_boundary_jacobian_m3;
} fusion_flux_volume_diagnostic_v1;
/* Enclosed toroidal volume of the supplied symmetric Fourier/core snapshot.
 * Layout, units, regular-core model and domain are fusion_flux_core.h's.
 * V = pi integral R(theta)^2 dZ/dtheta dtheta. R-independent constant terms
 * integrate to zero and are removed before summation to reduce cancellation.
 * The zero-mean radius is recovered from its angular derivative by Fourier
 * projection, avoiding subtraction of almost equal absolute radii near axis.
 * The uniform periodic rule uses a power-of-two node count >3*harmonics;
 * it integrates this finite Fourier polynomial exactly in exact arithmetic.
 * Floating-point field evaluation and very small/thin volumes remain subject
 * to conditioning. This is not an interval-arithmetic error bound.
 * Axis volume is exactly zero. Otherwise all sampled boundary Jacobians must
 * have one nonzero orientation consistent with the signed enclosed volume.
 * This checks the requested boundary only: no global nesting/self-intersection,
 * interior Jacobian, equilibrium, field epoch, source volume, orbit or loss
 * validation follows. Shell volume is a caller's difference of enclosed values,
 * with its own conditioning and monotonicity checks. No rescaling is performed.
 * Read-only/stateless/reentrant; inputs must remain live and nonoverlapping
 * with output, as in the underlying field API. Errors clear output. No I/O,
 * host indexing, transport mutation, reaction-cache or other physics change.
 */
int fusion_c_flux_core_enclosed_volume_diagnostic_v1(int nodes,int harmonics,
 const double*xi,const double*coeff_m,const double*slopes_m,const double*F_T_m,
 const double*C_T_m2,int match_node,double point_xi,
 fusion_flux_volume_diagnostic_v1*out);
#ifdef __cplusplus
}
#endif
#endif
