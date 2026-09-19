#ifndef FUSION_BOUNDARY_SEGMENT_H
#define FUSION_BOUNDARY_SEGMENT_H
#include "fusion_flux_boundary.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_SEGMENT_UNRESOLVED=0, FUSION_SEGMENT_CLEAR=1,
       FUSION_SEGMENT_GEOMETRIC_CANDIDATE=2 };
enum { FUSION_SEGMENT_FINISHED=0, FUSION_SEGMENT_FRACTION_LIMIT=1,
       FUSION_SEGMENT_BUDGET_LIMIT=2, FUSION_SEGMENT_INITIAL_AMBIGUOUS=3,
       FUSION_SEGMENT_RESOLUTION_LIMIT=4 };
typedef struct fusion_boundary_segment_value_v1 {
 double fraction_lo,fraction_hi,boundary_band_m;
 int result,tolerance_met,termination_reason,intervals_visited;
} fusion_boundary_segment_value_v1;
/* Cartesian straight segment x(s)=start+s*(end-start), s in[0,1]. This is a
 * GEOMETRIC query only: no field, trajectory integration or physical loss.
 * Boundary inputs/allowances as fusion_flux_boundary.h; geometry allowance
 * must additionally cover Cartesian interpolation/length arithmetic.
 * fraction_tolerance in(0,1], max_intervals1..1000000; fixed stack depth<=52.
 * Start must be INSIDE; AMBIGUOUS gives UNRESOLVED, OUTSIDE gives OUT_OF_RANGE.
 * Left-first subdivision certifies entire inside subsegments by endpoint
 * clearance and Cartesian length (RZ projection is nonexpansive). Endpoints
 * alone never certify no contact. CLEAR means the whole segment passed.
 * CANDIDATE has safe INSIDE/OUTSIDE endpoints and all earlier fractions covered
 * by inside subsegments; its interval contains the first possible contact,
 * but need not have only one crossing. It is NOT a physical first-exit event.
 * A retained ancestor candidate may be wider than requested: tolerance_met=0
 * explicitly reports this at ambiguity/budget/resolution limits. A bounded
 * outside-endpoint refinement may shrink hi; sampled inside/ambiguous points
 * never advance the proven prefix. intervals_visited counts traversal nodes
 * plus outside-refinement queries, all charged to the same budget. UNRESOLVED
 * has no known safe outside endpoint. No guess resolves an ambiguous query.
 * All errors clear outputs. Stateless/no allocation/I/O, no input/output overlap.
 */
int fusion_c_boundary_segment(int harmonics,const double*coeff_m,int segments,
 const double*start_m,const double*end_m,double geometry_allowance_m,
 double angular_allowance_m2,double fraction_tolerance,int max_intervals,
 fusion_boundary_segment_value_v1*out);
#ifdef __cplusplus
}
#endif
#endif
