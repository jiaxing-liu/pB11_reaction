#ifndef FUSION_PROMPT_ORBIT_H
#define FUSION_PROMPT_ORBIT_H
#include "fusion_magnetic_push.h"
#include "fusion_boundary_segment.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_PROMPT_UNRESOLVED=0, FUSION_PROMPT_RETAINED=1,
       FUSION_PROMPT_EVENT=2, FUSION_PROMPT_CALLBACK_CONTRACT=1004 };
enum { FUSION_PROMPT_FINISHED=0, FUSION_PROMPT_BOUNDARY_UNRESOLVED=1,
       FUSION_PROMPT_STEP_BUDGET=2, FUSION_PROMPT_REFINEMENT_LIMIT=3,
       FUSION_PROMPT_TIME_RESOLUTION=4 };
typedef int (*fusion_orbit_point_callback_v1)(void*,const double*,int*);
typedef int (*fusion_orbit_segment_callback_v1)(void*,const double*,const double*,
 double fraction_tolerance,fusion_boundary_segment_value_v1*);
typedef struct fusion_prompt_options_v1 {
 double horizon_s,initial_dt_s,geometry_time_budget_s;
 double event_time_tolerance_s,event_position_tolerance_m;
 double final_position_tolerance_m,final_u_relative_tolerance;
 int min_levels,max_levels,max_steps_per_level;
} fusion_prompt_options_v1;
typedef struct fusion_prompt_result_v1 {
 int outcome,reason,levels_used,steps_last_level,phase;
 double final_time_s;
 fusion_orbit_state_v1 final_state;
 double event_time_lo_s,event_time_hi_s,event_x_lo_m[3],event_x_hi_m[3];
 double estimated_time_spread_s,estimated_position_spread_m,estimated_u_relative_spread;
} fusion_prompt_result_v1;
/* Static E=0 collisionless numerical model-boundary driver, NOT physical wall
 * loss. Independent uniform runs use ceil(horizon/initial_dt)*2^level steps.
 * At least min_levels (2..20), at most max_levels<=20; per-level step cap<=1e6.
 * All tolerances positive finite, explicitly supplied before running.
 * Adjacent complete retained paths or localized event candidates must agree.
 * Event time spread includes both geometric half-widths; position spread is
 * maximum distance over all four endpoint pairs. No Richardson extrapolation.
 * Agreement is an EMPIRICAL numerical test, NOT a rigorous true-error bound.
 * EVENT contains only the fine split-path geometric bracket plus spreads;
 * it is not an enclosure of the continuous orbit. No exit momentum is provided.
 * RETAINED is finite-horizon numerical retention, not permanent confinement.
 * UNRESOLVED carries only reason/counters; no valid final state or loss weight.
 * Callbacks must be deterministic for an immutable static physical context and
 * must not throw. Counters may change; physical state may not. point uses
 * FUSION_BOUNDARY_{INSIDE,OUTSIDE,AMBIGUOUS}. segment has D220 semantics and may
 * refine geometry; its requested internal tolerance is separate from the
 * external geometry_time_budget. Endpoint classifications are cross-checked,
 * but the driver cannot independently prove the callback's whole-segment claim.
 * Field is queried only after the first half-segment and endpoint are INSIDE.
 * Callback errors propagate unchanged; malformed claims return1004. No outside
 * field evaluation, clipping, extrapolation, hidden state, allocation or I/O.
 * Nonoverlapping inputs/output; caller initial state never modified. Nonzero
 * API status clears output (outcome zero is UNRESOLVED, never RETAINED).
 */
int fusion_c_prompt_orbit(double mass_kg,double charge_C,
 const fusion_orbit_state_v1*initial,const fusion_prompt_options_v1*options,
 fusion_magnetic_field_callback_v1 field,fusion_orbit_point_callback_v1 point,
 fusion_orbit_segment_callback_v1 segment,void*context,fusion_prompt_result_v1*out);
#ifdef __cplusplus
}
#endif
#endif
