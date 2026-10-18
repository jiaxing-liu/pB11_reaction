#ifndef FUSION_THERMAL_BIRTH_BUDGET_H
#define FUSION_THERMAL_BIRTH_BUDGET_H
#include "fusion_thermal_birth.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit numerical accuracy policy. version=1, reserved=0. Tolerances
 * are normalized L1 in [0,2], against the SAME finite full-order source.
 * They do not bound nuclear uncertainty or subsequent plasma response.
 * Both coarse orders are4..corresponding full order; all fields validated.
 * Either zero tolerance retains every parent at full order. */
typedef struct fusion_thermal_birth_budget_policy_v1 {
 int version,coarse_nq,coarse_ncos,reserved;
 double number_L1_tolerance,energy_L1_tolerance;
} fusion_thermal_birth_budget_policy_v1;
typedef struct fusion_thermal_birth_budget_report_v1 {
 int version,full_nq,full_ncos,coarse_nq,coarse_ncos,mixed;
 uint64_t relative_nodes,cm_nodes,joint_nodes,high_pairs,low_pairs;
 uint64_t high_relative_nodes,low_relative_nodes;
 double number_L1_bound,energy_L1_bound;
 double number_roundoff_allowance,energy_roundoff_allowance;
 double reaction_weight_m3_s,product_energy_weight_J_m3_s;
} fusion_thermal_birth_budget_report_v1;
/* pB only, two equal-temperature isotropic zero-drift Maxwellians. Same
 * physical options/grid/units/finite-cutoff rules as fusion_thermal_birth_grid.
 * No positive representable parent from that reference is dropped: selected
 * parents receive a cheaper positive conservative decay rule, including alpha0.
 * Rates, reactant debits and Q use the unchanged full relative/CM quadrature.
 *
 * Reported bounds are ideal positive-mixture bounds (cells AND both tails).
 * Roundoff allowances are separately declared engineering allowances1e-12,
 * not interval-arithmetic certificates; check observed errors and convergence.
 * Zero budgets reproduce the original full-order evaluator. A zero source
 * succeeds with zero weights/bounds/counts apart from quadrature metadata.
 * No library cache/table/BALDUR policy is silently changed by this call.
 *
 * All three outputs are required; nonnull report/out clear on any failure.
 * Nonnull birth clears only if cells1..100000. No overlapping arrays/structs.
 * Null policy/edges is INVALID_ARGUMENT; null required output is NULL_OUTPUT.
 * Unsupported policy version/reserved/order/nonfinite controls reject. All
 * options keep the original domain guards, including12MeV parent support.
 * Stateless, no hidden globals/environment/sampling. */
int fusion_c_pb_thermal_birth_budget_grid(double kT_J,
 const fusion_thermal_birth_options_v1*options,
 const fusion_thermal_birth_budget_policy_v1*policy,int cells,
 const double*edges_J,double*birth,fusion_thermal_birth_v1*out,
 fusion_thermal_birth_budget_report_v1*report);
#ifdef __cplusplus
}
#endif
#endif
