#ifndef FUSION_COUPLED_SOURCES_H
#define FUSION_COUPLED_SOURCES_H
#include "fusion_coupled_fast.h"
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
#ifdef __cplusplus
}
#endif
#endif
