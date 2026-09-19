#ifndef FUSION_FLUX_FIELD_H
#define FUSION_FLUX_FIELD_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Additional status: coordinate singularity, NOT physical loss/outside plasma. */
enum { FUSION_FLUX_FIELD_MODEL_REVISION = 1,
       FUSION_FIELD_AXIS_COORDINATE_SINGULAR = 1001 };
typedef struct fusion_flux_field_value_v1 {
 double R_m,Z_m,R_xi_m,Z_xi_m,R_theta_m,Z_theta_m;
 double signed_jacobian_m3,BR_T,BZ_T,Bphi_T;
} fusion_flux_field_value_v1;
/* Stateless up/down symmetric Fourier snapshot evaluator, SI units.
 * nodes 2..10000; harmonics 1..64; xi[0]=0, strictly increasing to xi[last].
 * coefficient arrays are node-major, stride=1+2*harmonics:
 * [R0,Rcos(1..harmonics),Zsin(1..harmonics)]. Values and d/dxi slopes in m.
 * Geometry uses cubic Hermite; its radial derivative is differentiated from
 * that SAME polynomial. F=R*Bphi [T m] and C [T m2] are flux functions,
 * explicitly linearly interpolated; no separate B-component interpolation.
 * R=R0+sum Rcos*cos(m theta); Z=sum Zsin*sin(m theta).
 * J=R*(Ztheta*Rxi-Rtheta*Zxi); BR=C*Rtheta/J, BZ=C*Ztheta/J.
 * Theta is analytic radians, independent of any host FFT storage ordering.
 * Signed F/C and geometry handedness are caller model conventions; no inferred
 * device direction. Negative J is retained, never changed to abs(J).
 * Revision 1 fixes Hermite geometry / linear F,C / analytic Fourier policy.
 * OUT_OF_RANGE is grid domain, not a certified LCFS/wall crossing.
 * No extrapolation. xi==0 returns the dedicated singular status. Full-domain
 * Cartesian evaluation/axis regularization/inverse mapping are NOT provided
 * by this stage; failed evaluation must never be counted as particle escape.
 * Local conditioning rejects |det|<=64*DBL_EPSILON*|tangent|*|radial|.
 * Does not certify global nesting, nonintersection or physical equilibrium.
 * All arrays finite; arrays/output nonoverlapping; input lengths as above.
 * Errors clear output. No state, allocations, file I/O or host dependency.
 */
int fusion_c_flux_field_eval(int nodes,int harmonics,const double*xi,
 const double*coeff_m,const double*slopes_m,const double*F_T_m,
 const double*C_T_m2,double point_xi,double theta_rad,
 fusion_flux_field_value_v1*out);
#ifdef __cplusplus
}
#endif
#endif
