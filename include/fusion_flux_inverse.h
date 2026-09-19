#ifndef FUSION_FLUX_INVERSE_H
#define FUSION_FLUX_INVERSE_H
#include "fusion_flux_core.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_FIELD_NOT_CONVERGED = 1002 };
typedef struct fusion_flux_inverse_value_v1 {
 double u,v,xi,theta_rad,residual_m;
 int iterations;
} fusion_flux_inverse_value_v1;
/* Regular coordinates u=xi*cos(theta), v=xi*sin(theta), using explicit core.
 * Same immutable snapshot/units as fusion_flux_core.h. caller tolerance_m>0,
 * max_iterations 1..100. seed_mode0 uses the axis linear map, mode1 uses the
 * finite caller seed_u/v; seeds and line-search trials are projected onto the
 * closed exported xi disk (never field-extrapolated). At most24backtracks/step.
 * Axis physical point is recognized by exact R==Raxis and Z==0, not a core hole.
 * Success means physical position residual<=tolerance_m; angle is canonical
 * atan2(v,u) (0 at exact axis). Caller controls tolerance/iteration budget.
 * Does not certify unique inverse or topology. Caller can compare finite seeds.
 * NOT_CONVERGED covers exhausted/stalled Newton and may include external points.
 * This API does NOT classify inside/outside: numerical failure is never escape.
 * OUT_OF_RANGE from validation concerns arguments/snapshot domain only. A future
 * independently validated boundary classifier must own actual exit events.
 * Errors clear output, no hidden state/allocation/I/O, inputs/output nonoverlap.
 */
int fusion_c_flux_inverse(int nodes,int harmonics,const double*xi,
 const double*coeff_m,const double*slopes_m,const double*F_T_m,const double*C_T_m2,
 int match_node,double R_m,double Z_m,double tolerance_m,int max_iterations,
 int seed_mode,double seed_u,double seed_v,fusion_flux_inverse_value_v1*out);
#ifdef __cplusplus
}
#endif
#endif
