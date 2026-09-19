#ifndef FUSION_COUPLED_SOURCES_H
#define FUSION_COUPLED_SOURCES_H
#include "fusion_coupled_fast.h"
#include "fusion_birth_floor.h"
#include "fusion_beam_birth_table.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Borrowed immutable table for one canonical channel/projectile orientation
 * and zero-based energy cell. The caller retains ownership throughout trial. */
typedef struct fusion_beam_table_entry_v1 {
 int channel, projectile_slot, energy_cell;
 const fusion_beam_birth_table_v1 *table;
} fusion_beam_table_entry_v1;
typedef struct fusion_beam_table_usage_v1 {
 int direct_evaluations, table_evaluations;
 double max_validated_rate_error, max_validated_debit_error;
 double max_validated_number_L1, max_validated_energy_L1;
} fusion_beam_table_usage_v1;
/* Additive full-source coupled trial. Null thermal_tables selects direct thermal
 * births; otherwise pass the existing five thermal table handles. Beam entries
 * form an explicit sparse set (0..10*cells entries). Missing entries use direct
 * full beam integration. Present entries MUST exactly match channel, slot,
 * energy-cell center, complete beam options and all grid edges. Invalid, null,
 * duplicate, disabled-channel, or redundant identical-reactant slot1 entries
 * reject, even when the corresponding fast population is zero. Used tables
 * reject temperatures outside their domain; no fallback masks a bad table.
 *
 * The interpolated full source supplies BOTH target-network coefficients and
 * accepted product births. Existing burn/FP/handoff splitting, spill checks,
 * conservation and output atomicity are unchanged. No table construction or
 * cache mutation occurs, and no physical state is advanced. Count0 reproduces
 * the corresponding direct-fast arithmetic. effective_charge is exactly0/1.
 *
 * usage counts actual source calls (including uncached second passes) and
 * reports sampled interpolation envelopes of USED tables. These are distinct
 * from result.max_source_* (sampled direct quadrature discrepancies). Tables'
 * explicit constructor controls govern interpolation; no new default tolerance.
 * Neither envelope is a uniform error guarantee. Coupled convergence is needed.
 * All outputs, diagnostics and usage clear on failure. No existing ABI changes.
 */
int fusion_c_coupled_sources_trial(double dt_s,
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
 fusion_beam_table_usage_v1 *usage);
/* Provisional opt-in numerical boundary. Limits apply to post-reactant-debit
 * ion reservoir BEFORE numerical correction; first center and ion kT are derived internally. All below-grid
 * nuclear births stay kinetic. Nuclear ledger energy stays physical; floor
 * ledger correction is separate from heat/Q. All outputs clear on failure.
 * Above-grid spill retains strict rejection. Existing entry points unchanged.
 * Non-floor per-cell source rounding admits only both-unrepresentable packet/rate
 * tails, with upward-bounded trial aggregate missing N/E below half a binary64
 * subnormal quantum in SI. Mapped floor packets retain strict positive-rate checks.
 * See docs/BIRTH_FLOOR.md; coupled physical convergence remains required. */
typedef struct fusion_coupled_floor_limits_v1 {
 double max_center_over_ion_kT, max_ion_energy_fraction;
} fusion_coupled_floor_limits_v1;
int fusion_c_coupled_sources_floor_trial(double dt_s,
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
 fusion_birth_floor_ledger_v1 *floor_ledger);

/* Explicit domain coverage policy. STRICT preserves the existing used-table
 * out-of-range error. DIRECT_OUTSIDE uses the original direct full source only
 * when the actual post-thermal-burn ion kT is outside a valid matching table.
 * Invalid/mismatched/duplicate entries and table evaluation errors still reject.
 * No temperature clamping, extrapolation or relaxed numerical gate.
 * Both floor pointers null disables the provisional floor; otherwise both are
 * required. All ordinary outputs and the mandatory outside-domain counter clear
 * on failure. The counter counts actual direct source calls caused by a present
 * out-of-domain table, and is a subset of usage.direct_evaluations. Missing
 * entries remain ordinary direct calls. Existing APIs and structs unchanged. */
enum { FUSION_BEAM_TABLE_STRICT=0, FUSION_BEAM_TABLE_DIRECT_OUTSIDE=1 };
int fusion_c_coupled_sources_covered_trial(double dt_s,
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
 int table_domain_policy, uint64_t *outside_domain_direct_evaluations);
/* Additive source-only TRIAL packet, not an accepted snapshot or new owner.
 * source0=thermal, source1=fast-target; channel0..4, species0..6(neutron).
 * mapped array layout [2][5][7][cells], AMOUNTS in m^-3, arithmetic centers.
 * Below/above retain ORIGINAL physical N/E before optional floor remapping.
 * No external birth, FP evolution, handoff or numerical floor correction is
 * folded into this packet. Fast angles/pitch are NOT recovered by this API.
 * Events[2][5] distinguish thermal from fast-target reaction amounts.
 * Packet is derived from the exact quadratures/tables and accepted thermal/
 * target-network subsolves already used by this trial; no second source solve.
 * Caller must publish ONLY when its enclosing host attempt is accepted, discard
 * on rejection, and multiply by the explicit source volume exactly once.
 * No callback or persistent side effect. All outputs clear on failure; mapped
 * clears only for valid cells1..100000. Nonoverlapping output arrays required.
 * Opt-in extraction may reject unrepresentable positive per-packet values even
 * when legacy summed quantities are representable. No tiny packet clipping.
 * Existing entry points/ABI/arithmetic remain unchanged. Extra workspace is
 * O(70*cells) long doubles plus O(70*cells) doubles, allocated only on opt-in.
 */
typedef struct fusion_birth_packets_v1 {
 double events_m3[2][5];
 double below_number_m3[2][5][7],below_energy_J_m3[2][5][7];
 double above_number_m3[2][5][7],above_energy_J_m3[2][5][7];
} fusion_birth_packets_v1;
int fusion_c_coupled_sources_packets_trial(double dt_s,
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
 double *mapped_packets_m3, fusion_birth_packets_v1 *packets);
#ifdef __cplusplus
}
#endif
#endif
