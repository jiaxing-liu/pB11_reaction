#ifndef FUSION_COUPLED_SOURCES_BUDGET_H
#define FUSION_COUPLED_SOURCES_BUDGET_H
#include "fusion_coupled_sources.h"
#include "fusion_birth_table_budget.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit pB thermal-source numerical policy on the complete covered/packet
 * trial. Other thermal channels and every fast-target source remain full.
 * No change to existing options/results/entries or implicit activation.
 * These estimates concern the SAME finite source reference, not nuclear data,
 * continuous-source accuracy or subsequent plasma response.
 * mode=0: no pB thermal source evaluated (disabled or absent reactants).
 * mode=1: direct pB call; direct contains its actual report, table is zero.
 * mode=2: immutable pB budget table; table contains its sampled envelopes,
 * direct is zero (no fictitious per-call quadrature report).
 * number/energy_L1_estimate are direct ideal bound + declared engineering
 * allowance, or the table's sampled COMPOSITE interpolation/source estimate.
 * For a zero-weight direct source both estimates are zero.
 * version=1 on whole-trial success, including inactive mode0; all bytes clear
 * on failure. Observations are not accepted state; host acceptance is required.
 */
typedef struct fusion_coupled_source_accuracy_v1 {
 int version, mode;
 double number_L1_estimate, energy_L1_estimate;
 fusion_thermal_birth_budget_report_v1 direct;
 fusion_birth_table_budget_info_v1 table;
} fusion_coupled_source_accuracy_v1;
/* pb_policy and accuracy are mandatory, distinct from all other inputs/outputs.
 * The policy is validated even when pB is inactive. With thermal_tables null,
 * pB uses the direct budget API. Otherwise an active pB needs an exactly matching
 * BUDGET table; full pB tables and mismatched policies reject. Non-pB entries
 * retain the original full-table matching and pb_low normalization. Out-of-table
 * thermal temperatures reject; beam coverage policy is unchanged. All floor,
 * packet, handoff, conservation, weak-source and atomic failure rules from the
 * original packet entry apply. No extra source solve/cache mutation occurs.
 * Old entries continue rejecting budget tables, including zero-budget ones.
 */
/* Covered source accuracy is independent of packet extraction. This entry has
 * the covered trial's physical behavior, including its legacy weak-source
 * rounding policy; it does not require positive packets to be representable.
 * pb_policy and accuracy are mandatory and use the policy/table rules above.
 */
int fusion_c_coupled_sources_covered_budget_trial(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,
 const fusion_birth_table_v1 *const *thermal_tables,
 int beam_table_count, const fusion_beam_table_entry_v1 *beam_tables,
 int effective_charge, int cells, const double *edges_J,
 const double thermal_number_m3[6], double electron_energy_J_m3,
 double ion_energy_J_m3, double electron_density_m3,
 const double thermal_charge_squared[6], int inert_count,
 const fusion_inert_ion_v1 *inert, const double *coulomb_logs,
 const double *old_s_m3, const double *old_t_m3,
 const double *external_birth_m3_s, const double *escape_s_inv,
 double new_thermal_m3[6], double *new_s_m3, double *new_t_m3,
 fusion_coupled_thermal_v1 *out, fusion_handoff_diagnostics_v1 *diagnostics,
 fusion_beam_table_usage_v1 *usage,
 const fusion_coupled_floor_limits_v1 *floor_limits,
 fusion_birth_floor_ledger_v1 *floor_ledger,
 int table_domain_policy, uint64_t *outside_domain_direct_evaluations,
 const fusion_thermal_birth_budget_policy_v1 *pb_policy,
 fusion_coupled_source_accuracy_v1 *accuracy);
int fusion_c_coupled_sources_packets_budget_trial(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,
 const fusion_birth_table_v1 *const *thermal_tables,
 int beam_table_count, const fusion_beam_table_entry_v1 *beam_tables,
 int effective_charge, int cells, const double *edges_J,
 const double thermal_number_m3[6], double electron_energy_J_m3,
 double ion_energy_J_m3, double electron_density_m3,
 const double thermal_charge_squared[6], int inert_count,
 const fusion_inert_ion_v1 *inert, const double *coulomb_logs,
 const double *old_s_m3, const double *old_t_m3,
 const double *external_birth_m3_s, const double *escape_s_inv,
 double new_thermal_m3[6], double *new_s_m3, double *new_t_m3,
 fusion_coupled_thermal_v1 *out, fusion_handoff_diagnostics_v1 *diagnostics,
 fusion_beam_table_usage_v1 *usage,
 const fusion_coupled_floor_limits_v1 *floor_limits,
 fusion_birth_floor_ledger_v1 *floor_ledger,
 int table_domain_policy, uint64_t *outside_domain_direct_evaluations,
 double *mapped_packets_m3, fusion_birth_packets_v1 *packets,
 const fusion_thermal_birth_budget_policy_v1 *pb_policy,
 fusion_coupled_source_accuracy_v1 *accuracy);
#ifdef __cplusplus
}
#endif
#endif
