#ifndef FUSION_BIRTH_TABLE_BUDGET_H
#define FUSION_BIRTH_TABLE_BUDGET_H
#include "fusion_birth_table.h"
#include "fusion_thermal_birth_budget.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit pB-only common-T source-accuracy table. Same immutable handle,
 * grid/evaluate/destroy/pack/digest interfaces as full-source tables.
 * Old request/coupled matching MUST reject budget-enabled tables, even with
 * zero budgets; use the explicit matching entry below. Full-source constructors
 * and format1 stay unchanged. Budget tables serialize as format2 with policy,
 * sampled error accounting and stored-knot source reports; no source reruns
 * occur on import. These are sampled numerical checks, not supremum/physical
 * uncertainty certificates. */
typedef struct fusion_birth_table_budget_info_v1 {
 fusion_birth_table_info_v1 base;
 fusion_thermal_birth_budget_policy_v1 policy;
 /* Effective sampled source bounds include the separately reported
  * engineering allowance for nonzero sources (zero sources have zero bound). */
 double max_sampled_source_number_bound,max_sampled_source_energy_bound;
 double max_sampled_number_roundoff_allowance,max_sampled_energy_roundoff_allowance;
 double max_validated_interpolation_number_L1,max_validated_interpolation_energy_L1;
} fusion_birth_table_budget_info_v1;
/* base.max_validated_number/energy_L1 include sampled interpolation error
 * normalized to expected source N/E PLUS the direct-reference mixture bound
 * and its declared engineering roundoff allowance. Endpoint source budgets
 * must also fit controls. Rate/debit gates unchanged. Interpolation-only
 * maxima are retained separately. No off-sample guarantee is implied.
 * Clears*out before all validation; no partial table escapes. */
int fusion_c_pb_birth_table_create_budget(double lower_kT_J,double upper_kT_J,
 const fusion_thermal_birth_options_v1*source,
 const fusion_thermal_birth_budget_policy_v1*policy,
 const fusion_birth_table_control_v1*control,int cells,const double*edges_J,
 fusion_birth_table_v1**out);
/* Budget-only introspection. Nonnull out clears on failure; full-source tables
 * reject rather than acquiring an implicit policy. Generic old info returns
 * the base information with the combined bounds described above. */
int fusion_c_pb_birth_table_budget_info(const fusion_birth_table_v1*table,
 fusion_birth_table_budget_info_v1*out);
int fusion_c_pb_birth_table_budget_matches_request(const fusion_birth_table_v1*table,
 double lower_kT_J,double upper_kT_J,const fusion_thermal_birth_options_v1*source,
 const fusion_thermal_birth_budget_policy_v1*policy,
 const fusion_birth_table_control_v1*control,int cells,const double*edges_J,int*matches);
/* Zero-based stored knot index. Both outputs required and clear on failure. */
int fusion_c_pb_birth_table_budget_knot_report(const fusion_birth_table_v1*table,
 int index,double*kT_J,fusion_thermal_birth_budget_report_v1*out);
#ifdef __cplusplus
}
#endif
#endif
