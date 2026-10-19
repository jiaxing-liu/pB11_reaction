#ifndef FUSION_CAPTURE_REGISTRY_H
#define FUSION_CAPTURE_REGISTRY_H
#include "fusion_birth_table.h"
#include "fusion_beam_birth_table.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Optional diagnostic ownership, independent of the numerical kernel.
 * Handles created here MUST be released with the matching diagnostic destroy
 * hook, never the original library destroy. Context destruction releases all
 * remaining handles. Do not destroy a context concurrently with any operation
 * or use any returned handle after removal/context destruction. Immutable
 * numerical evaluation can borrow registered handles while their owner lives.
 * Existing handles are supported only through the explicit ownership and
 * packed-content identity constrained adopt APIs below. No arbitrary-pointer
 * validity check is possible; these APIs do not change numerical physics.
 * Paths are absolute, <=4096 bytes, and identify the exact supplied packed
 * bytes. Replay must independently verify the referenced file against digest.
 * This registry does not certify table interpolation accuracy.
 */
typedef struct fusion_capture_context_v1 fusion_capture_context_v1;
enum fusion_capture_status_v1 {
 FUSION_CAPTURE_OK=0, FUSION_CAPTURE_INVALID=100,
 FUSION_CAPTURE_ALLOCATION=101, FUSION_CAPTURE_CAPACITY=102,
 FUSION_CAPTURE_UNKNOWN_HANDLE=103, FUSION_CAPTURE_UNPACK_REJECTED=104
};
enum fusion_capture_table_kind_v1 {
 FUSION_CAPTURE_THERMAL_TABLE=1, FUSION_CAPTURE_BEAM_TABLE=2
};
typedef struct fusion_capture_table_identity_v1 {
 int kind;
 uint64_t packed_bytes;
 char content_sha256[65], kernel_identity[65], qualified_path[4097];
} fusion_capture_table_identity_v1;
/* On failure outputs clear. Nonzero capacity, at most 1000000 entries. */
int fusion_capture_context_create_v1(uint64_t capacity,
 fusion_capture_context_v1 **out);
void fusion_capture_context_destroy_v1(fusion_capture_context_v1 *context);
/* Return a diagnostic status, not the numerical status. unpack_status is
 * separately set to the actual original unpack status if unpack was called,
 * otherwise -1. No handle is published until digest/registration succeed.
 * Unpublished handles are destroyed on registration failure. Inputs unchanged.
 */
int fusion_capture_thermal_unpack_v1(fusion_capture_context_v1 *context,
 const void *bytes,size_t length,const char *qualified_path,
 fusion_birth_table_v1 **out,int *unpack_status);
int fusion_capture_beam_unpack_v1(fusion_capture_context_v1 *context,
 const void *bytes,size_t length,const char *qualified_path,
 fusion_beam_birth_table_v1 **out,int *unpack_status);
/* Adopt a valid, sole-owned typed raw handle, not registered in another
 * context. Success transfers ownership without replacing/reconstructing the
 * handle: use the matching registry destroy hook thereafter. Failure never
 * destroys the raw handle, changes registered entries, or modifies inputs.
 * pack_status is mandatory: -1 until the core content_digest is called, then
 * its actual status (separate from the diagnostic return). Digest/size mismatch
 * or core rejection returns INVALID. Duplicate returns INVALID; capacity returns
 * CAPACITY. Caller must already have saved the exact packed bytes to the
 * absolute qualified_path (<=4096 bytes); no file is read here. Replay verifies
 * the file independently. Limits are thermal 512MiB and beam 256MiB.
 */
int fusion_capture_thermal_adopt_v1(fusion_capture_context_v1 *context,
 fusion_birth_table_v1 *table,const void *packed_bytes,size_t length,
 const char *qualified_path,int *pack_status);
int fusion_capture_beam_adopt_v1(fusion_capture_context_v1 *context,
 fusion_beam_birth_table_v1 *table,const void *packed_bytes,size_t length,
 const char *qualified_path,int *pack_status);
int fusion_capture_table_identity_v1_get(fusion_capture_context_v1 *context,
 int kind,const void *handle,fusion_capture_table_identity_v1 *out);
/* Wrong kind/unknown handle rejects without releasing or clearing the handle. */
int fusion_capture_thermal_destroy_v1(fusion_capture_context_v1 *context,
 fusion_birth_table_v1 **handle);
int fusion_capture_beam_destroy_v1(fusion_capture_context_v1 *context,
 fusion_beam_birth_table_v1 **handle);
#ifdef __cplusplus
}
#endif
#endif
