#ifndef FUSION_CAPTURE_TRIAL_H
#define FUSION_CAPTURE_TRIAL_H
#include "fusion_capture_registry.h"
#include "fusion_coupled_sources.h"
#ifdef __cplusplus
extern "C" {
#endif
enum fusion_capture_entry_v1 { FUSION_CAPTURE_COVERED=1, FUSION_CAPTURE_PACKETS=2 };
enum fusion_capture_file_status_v1 {
 FUSION_CAPTURE_FORMAT=105, FUSION_CAPTURE_IO=106,
 FUSION_CAPTURE_IDENTITY=107, FUSION_CAPTURE_UNSUPPORTED=108
};
/* Borrowed exact inputs of one local trial, canonical species-major arrays.
 * This is a C ABI view only: padding, pointers and native integers are NEVER
 * serialized. All required pointers must be valid and live throughout capture.
 * The numerical API's non-overlapping buffer contract applies. Unsupported
 * aliased input layouts reject capture separately; pointer addresses are not
 * used as persistent identities. Output buffers must satisfy the original API.
 * cells1..100000, inert_count0..32, beam_count0..10*cells. Optional pointer
 * presence is retained, including count-zero inert/beam pointers and the
 * distinction between direct thermal and an all-null table-handle array.
 * Invalid dimensions/null required inputs are a separate capture failure;
 * the original numerical status is never replaced. No output state is recorded.
 */
typedef struct fusion_capture_trial_input_v1 {
 int entry_point, original_status, host_zone;
 double host_time_s, dt_s;
 const fusion_coupled_thermal_options_v1 *options;
 const fusion_fast_target_options_v1 *fast_options;
 const fusion_birth_table_v1 *const *thermal_tables;
 int beam_table_count;
 const fusion_beam_table_entry_v1 *beam_tables;
 int effective_charge, cells;
 const double *edges_J, *thermal_number_m3;
 double electron_energy_J_m3, ion_energy_J_m3, electron_density_m3;
 const double *thermal_charge_squared;
 int inert_count;
 const fusion_inert_ion_v1 *inert;
 const double *coulomb_logs, *old_s_m3, *old_t_m3, *external_birth_m3_s, *escape_s_inv;
 const fusion_coupled_floor_limits_v1 *floor_limits;
 int table_domain_policy;
} fusion_capture_trial_input_v1;
typedef struct fusion_capture_replay_report_v1 {
 int original_status, replay_status, status_matches, entry_point, host_zone, cells;
 double host_time_s, dt_s;
 char kernel_identity[65], output_sha256[65];
} fusion_capture_replay_report_v1;
/* Exclusive atomic POSIX publication. max_bytes1..128MiB, including checksum.
 * Nonnull tables must have live matching registrations in context. Capture
 * reads the loader identities, not large re-packed table payloads. Context can
 * be null if no nonnull tables occur. Existing evidence is never overwritten.
 * I/O/format/resource/registry failures are diagnostic-only statuses. No
 * physics calls, coefficient/cache changes or acceptance occur in this API.
 */
int fusion_capture_write_trial_v1(fusion_capture_context_v1 *context,
 const char *path,const fusion_capture_trial_input_v1 *input,uint64_t max_bytes);
/* Validate schema/checksum/current kernel and every referenced table file's
 * exact bytes before invoking the RECORDED entry point with the original bits.
 * Returns diagnostic status; successful replay does not imply physics success.
 * report.replay_status is the actual original API result, possibly nonzero.
 * report.output_sha256 covers all explicitly encoded output fields/arrays
 * including failed-cleared outputs; these are not measurements of the failed
 * operation. Replay has no accepted host state and never changes its input file.
 */
int fusion_capture_replay_file_v1(const char *path,uint64_t max_bytes,
 fusion_capture_replay_report_v1 *report);
#ifdef __cplusplus
}
#endif
#endif
