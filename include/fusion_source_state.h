#ifndef FUSION_SOURCE_STATE_H
#define FUSION_SOURCE_STATE_H
#include "fusion_network.h"
#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_source_state_v1 fusion_source_state_v1;
/* Step AMOUNTS: number m^-3, energy J/m^3. Species order is the six-species
 * network. Seven baths: electron, then the six thermal ion species. Heat is
 * fast-species-major [6][7], positive to bath, and is the only signed field.
 * Nuclear birth, external birth and removed fuel are separate. Products are
 * born fast; handoff is a later, separately accounted transfer to the host.
 * No implicit thermalization of grid spill or conversion of missing data.
 */
typedef struct fusion_source_ledger_v1 {
 double events_m3[5];
 double nuclear_born_number_m3[6], nuclear_born_energy_J_m3[6];
 double external_born_number_m3[6], external_born_energy_J_m3[6];
 double thermal_consumed_number_m3[6], thermal_consumed_energy_J_m3[6];
 double fast_consumed_number_m3[6], fast_consumed_energy_J_m3[6];
 double escaped_number_m3[6], escaped_energy_J_m3[6];
 double handed_off_number_m3[6], handed_off_energy_J_m3[6];
 double heat_to_bath_J_m3[42];
 double neutron_number_m3, neutron_energy_J_m3;
} fusion_source_ledger_v1;
/* Per-zone owner of two kinetic component populations [species][cell].
 * Both are kinetic, not fluid ash. Grid edges in J; arithmetic-center energy.
 * Host owns thermal state and must commit/restore it atomically with this
 * context. Do not access the same context concurrently; independent contexts
 * share no mutable state. Input/output buffers must not overlap. model_tag is a caller-defined configuration fingerprint checked
 * on restart, not a security token or automatically inferred model identity.
 * cells in[1,1000000], finite increasing nonnegative edges; time finite>=0.
 * Initial populations finite/nonnegative. Explicit initial inventory anchors
 * cumulative conservation checks. No file I/O, global state or physics advance.
 */
int fusion_c_source_state_create(int cells,const double *edges_J,
 const double *initial_s_m3,const double *initial_t_m3,double initial_time_s,
 uint64_t model_tag,fusion_source_state_v1 **out);
void fusion_c_source_state_destroy(fusion_source_state_v1 *state);
/* Immutable layout query, including after unpack. cells is zero on errors. */
int fusion_c_source_state_cells(const fusion_source_state_v1 *state,int *cells);
/* Snapshot accepted state only; arrays each6*cells, ledger cumulative amounts.
 * Outputs required; caller knows cells from its configuration. Returns epoch.
 */
int fusion_c_source_state_snapshot(const fusion_source_state_v1 *state,
 double *accepted_s_m3,double *accepted_t_m3,fusion_source_ledger_v1 *cumulative,
 double *accepted_time_s,uint64_t *epoch);
/* A new begin supersedes any prior trial and produces a fresh ticket. No
 * accepted state changes. Invalid replacement begin also discards an older
 * trial. Tickets are context-local. Positive dt must advance representable time.
 * Stage REPLACES a trial (never accumulates nonlinear re-evaluations). It
 * validates stoichiometry, coherent new nuclear mass/Q energy, each fast
 * species' particle/energy balance, and cumulative representability.
 * Thermal fuel availability remains the caller/burn solver's responsibility.
 * Relative ledger residual tolerance1e-10, scaled to terms of each identity;
 * no unit-size absolute floor. A rejected stage invalidates the staged trial.
 * Commit requires current ticket and a valid stage, increments epoch once.
 * Repeated/stale commit rejected; discard removes only pending state.
 */
int fusion_c_source_state_begin(fusion_source_state_v1 *state,double dt_s,uint64_t *ticket);
int fusion_c_source_state_stage(fusion_source_state_v1 *state,uint64_t ticket,
 const double *trial_s_m3,const double *trial_t_m3,const fusion_source_ledger_v1 *step);
int fusion_c_source_state_commit(fusion_source_state_v1 *state,uint64_t ticket);
int fusion_c_source_state_discard(fusion_source_state_v1 *state,uint64_t ticket);
/* Portable versioned little-endian IEEE754 restart, with accidental-corruption
 * checksum. Pack rejects pending trials; no unaccepted result enters restart.
 * Unpack creates a NEW context, checks tag, dimensions, all ledgers and total
 * balances. Host must checkpoint its thermal/geometry state at same epoch/time.
 * bytes_written/required are zero on errors. No struct padding serialized.
 * Buffer lengths use size_t; caller frees packed storage, destroy frees state.
 */
int fusion_c_source_state_pack_size(const fusion_source_state_v1 *state,size_t *required);
int fusion_c_source_state_pack(const fusion_source_state_v1 *state,
 unsigned char *buffer,size_t capacity,size_t *bytes_written);
int fusion_c_source_state_unpack(const unsigned char *buffer,size_t length,
 uint64_t expected_model_tag,fusion_source_state_v1 **out);
#ifdef __cplusplus
}
#endif
#endif
