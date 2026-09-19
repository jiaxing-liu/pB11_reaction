#ifndef FUSION_FLUX_CORE_H
#define FUSION_FLUX_CORE_H
#include "fusion_flux_field.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_FLUX_CORE_MODEL_REVISION = 1 };
/* Explicit regular-core surrogate; does NOT replace the v1 interpolation.
 * Same snapshot layout/units as fusion_flux_field.h. match_node is zero-based,
 * 1..nodes-2, xc=xi[match_node]. q>=xc delegates exactly to v1.
 * q<xc uses t=q/xc: R0,F=axis+alpha*t^2+beta*t^4;
 * harmonic m>=2 and C(m=1): t^m*(alpha+beta*t^2);
 * geometry m=1: xc*t*(axis_slope+b*t^2+d*t^4).
 * Coefficients match node values/slopes at xc. F/C slopes are the RIGHT v1
 * segment slopes; geometry m=1 retains supplied axis slopes. Higher geometry
 * axis slopes are explicitly replaced by this new model, not silently clipped.
 * Requires zero axis harmonic values and C(0)=0, Raxis>0 and nondegenerate
 * axis m=1 mapping (min abs slope/max abs slope >64*DBL_EPSILON).
 * Axis physical field is (BR,BZ,Bphi)=(0,0,Faxis/Raxis), independent of theta.
 * Axis J and theta derivatives are zero; xi derivatives describe the polar
 * coordinate direction and can depend on theta. No inverse-map claim follows.
 * Geometry is C1 and field is continuous at xc, NOT generally C1 field gradient.
 * Core fold against axis orientation is a numerical failure; global topology
 * still requires validation. This is not a unique physical reconstruction;
 * matching-radius sensitivity is REQUIRED before scientific use.
 * No orbit, collision, loss, host-state or I/O. Errors clear output.
 */
int fusion_c_flux_core_eval(int nodes,int harmonics,const double*xi,
 const double*coeff_m,const double*slopes_m,const double*F_T_m,
 const double*C_T_m2,int match_node,double point_xi,double theta_rad,
 fusion_flux_field_value_v1*out);
#ifdef __cplusplus
}
#endif
#endif
