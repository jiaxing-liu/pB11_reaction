#ifndef FUSION_FLUX_BOUNDARY_H
#define FUSION_FLUX_BOUNDARY_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY = 1003 };
enum { FUSION_BOUNDARY_AMBIGUOUS = 0, FUSION_BOUNDARY_INSIDE = 1,
       FUSION_BOUNDARY_OUTSIDE = 2 };
typedef struct fusion_flux_boundary_value_v1 {
 double polygon_distance_m, uncertainty_m, angular_margin_m2;
 int classification;
} fusion_flux_boundary_value_v1;
/* Model LCFS only, NOT a material wall or physical particle-loss decision.
 * coeff=[R0,Rcos1..m,Zsin1..m] in metres, m1..64; segments16..65536.
 * Explicit aggregate error allowances must be finite and positive:
 * geometry_allowance_m covers curve/polygon evaluation, distances, bound
 * arithmetic and query roundoff; angular_allowance_m2 covers W and P*M*h/2
 * arithmetic/evaluation uncertainty. Caller must justify these allowances for
 * its platform/data. They do not represent device uncertainty. Computation is
 * long double, NOT outward-rounded interval arithmetic or a certified proof.
 * Supports only curves passing guarded positive angular derivative relative to
 * (R0,0), winding=1 and center-clearance checks. Unsupported is never outside.
 * Global second-derivative bound gives curve-to-chord deviation M*h*h/8.
 * Distance within this bound plus geometry allowance returns AMBIGUOUS.
 * Inverse mapping is neither used nor needed. All errors clear output;
 * no allocation, hidden state, I/O or host dependence. Inputs/output nonoverlap.
 */
int fusion_c_flux_boundary(int harmonics,const double*coeff_m,int segments,
 double R_m,double Z_m,double geometry_allowance_m,double angular_allowance_m2,
 fusion_flux_boundary_value_v1*out);
/* Immutable prepared geometry. Additive API; legacy point queries retain their
 * allocation-free implementation. Owns all derived data, never retains coeff.
 * Storage <=65536 vertices of four long doubles (platform-dependent, about4MiB
 * on x86-64). Query arithmetic/allowances are the legacy algorithm. No relaxed
 * geometry gate or spatial index. Queries are read-only and thread-safe; caller
 * must not destroy while any query runs. NULL destroy is safe; other handles
 * must be live handles returned by prepare (no arbitrary/freed-pointer probe).
 * Prepare sets *out=NULL before validation and publishes only on success.
 * Queries clear output on error. Allocation failure returns EXCEPTION.
 */
typedef struct fusion_flux_boundary_prepared fusion_flux_boundary_prepared;
int fusion_c_flux_boundary_prepare(int harmonics,const double*coeff_m,int segments,
 double geometry_allowance_m,double angular_allowance_m2,fusion_flux_boundary_prepared**out);
int fusion_c_flux_boundary_prepared_point(const fusion_flux_boundary_prepared*context,
 double R_m,double Z_m,fusion_flux_boundary_value_v1*out);
void fusion_c_flux_boundary_destroy(fusion_flux_boundary_prepared*context);
#ifdef __cplusplus
}
#endif
#endif
