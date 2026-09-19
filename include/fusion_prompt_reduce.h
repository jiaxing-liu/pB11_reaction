#ifndef FUSION_PROMPT_REDUCE_H
#define FUSION_PROMPT_REDUCE_H
#include "fusion_prompt_orbit.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_prompt_weight_v1 {
 double number_weight,kinetic_energy_J;
 int api_status,outcome;
} fusion_prompt_weight_v1;
typedef struct fusion_prompt_totals_v1 {
 double number[3],energy_J[3]; /* indexed by UNRESOLVED0/RETAINED1/EVENT2 */
 double total_number,total_energy_J;
 double number_fraction_low,number_fraction_high;
 double energy_fraction_low,energy_fraction_high;
 int number_fraction_defined,energy_fraction_defined;
} fusion_prompt_totals_v1;
/* Pure extensive one-particle ensemble reduction. N weight>=0, single-particle
 * K[J]>=0; energy weight is derived internally as N*K, never independently
 * supplied. Caller groups compatible species/channels/epochs/horizons and must
 * ensure K matches the marker's initial proper velocity and nuclear mass.
 * Does not create markers, run orbits, infer omitted weights or validate sampler
 * physics. Zero count accepts null items. Count0..1e6, inputs/output disjoint.
 * Any nonzero marker API status propagates unchanged and clears ALL output,
 * even for zero weight: errors cannot become physical UNRESOLVED. Invalid
 * outcome or negative/nonfinite weights reject. Input order is deterministic.
 * Category sums and total use compensated long-double accumulation. Nonzero
 * values unrepresentable as double fail instead of disappearing or clipping.
 * Bounds are event/total and (event+unresolved)/total for number and energy.
 * These describe only classification uncertainty in the SUPPLIED measure,
 * subject to floating rounding. No quadrature/model/confidence bound is implied.
 * Zero totals return fraction_defined=0 and numeric placeholders0, not0%loss.
 * Positive number with zero kinetic energy may have only number fraction defined.
 * Missing spectral/spatial tails must remain separately accounted by caller.
 * No FP/source mutation, hidden state, allocation or I/O; errors clear output.
 */
int fusion_c_prompt_reduce(int count,const fusion_prompt_weight_v1*items,
 fusion_prompt_totals_v1*out);
#ifdef __cplusplus
}
#endif
#endif
