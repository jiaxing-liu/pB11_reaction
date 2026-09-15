#ifndef FUSION_BEAM_BIRTH_TABLE_H
#define FUSION_BEAM_BIRTH_TABLE_H
#include "fusion_birth_table.h"
#include "fusion_beam_birth.h"
#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_beam_birth_table_v1 fusion_beam_birth_table_v1;
typedef struct fusion_beam_birth_table_info_v1 {
 double projectile_energy_J,lower_kT_J,upper_kT_J;
 double max_validated_rate_error,max_validated_debit_error;
 double max_validated_number_L1,max_validated_energy_L1;
 double max_sampled_direct_rate_discrepancy,max_sampled_direct_debit_discrepancy;
 int channel,projectile_slot,cells,knots,direct_evaluations;
 uint64_t spectral_entries_evaluated,stored_spectral_entries;
 fusion_beam_birth_options_v1 source;
 fusion_birth_table_control_v1 control;
} fusion_beam_birth_table_info_v1;
/* Immutable full-source table at FIXED channel, projectile slot and kinetic
 * energy. The interpolation coordinate is target kT[J], strictly0<lower<upper.
 * Options/grid are copied. Source coefficient units and canonical debit order
 * are those of fusion_c_beam_birth_grid; multiply by projectile and target
 * densities outside this API, without an identical-pair factor for a distinct
 * beam/thermal pool. No host state, I/O, hidden cache or state advance.
 *
 * Same adaptive geometric quarter/mid/three-quarter sampled gates as thermal
 * birth tables: rate, both debits and maximum per-product number/energy L1,
 * including spills. One convex log-kT weight applies to ALL source quantities.
 * Gates are sampled evidence, not a supremum error guarantee. Independent
 * off-sample and coupled convergence validation remain required.
 *
 * Internal knots store every strictly nonzero grid value, including subnormals;
 * only exact zeros are omitted. Public output is dense[7*cells]. No shape or
 * moment surrogate, renormalization, weak-source cutoff or extrapolation.
 * Construction caps cumulative nonzero entries in all direct samples at
 *12,500,000, in addition to explicit control evaluation/knot/depth limits.
 * This conservatively bounds sparse payload work/storage; metadata, vector
 * capacity and dense temporary buffers are additional. No partial table escapes.
 * cells1..100000; control domains as fusion_birth_table_control_v1. Exact
 * stored knots reproduce direct coefficients. Create clears*out; error eval
 * clears outputs for valid cells. Arrays must not overlap. No ABI change to
 * existing tables. Concurrent immutable info/evaluate is allowed; destroy once.
 */
int fusion_c_beam_birth_table_create(int channel,int projectile_slot,
 double projectile_energy_J,double lower_kT_J,double upper_kT_J,
 const fusion_beam_birth_options_v1 *source,const fusion_birth_table_control_v1 *control,
 int cells,const double *edges_J,fusion_beam_birth_table_v1 **out);
void fusion_c_beam_birth_table_destroy(fusion_beam_birth_table_v1 *table);
int fusion_c_beam_birth_table_info(const fusion_beam_birth_table_v1 *table,
 fusion_beam_birth_table_info_v1 *out);
int fusion_c_beam_birth_table_evaluate(const fusion_beam_birth_table_v1 *table,
 double target_kT_J,int cells,double *birth,fusion_birth_coefficients_v1 *out);
/* Stable process-lifetime ASCII SHA256 string (64 hex characters plus NUL).
 * Conservatively identifies source/header content and Boost version. It is
 * not binary provenance, authentication, or an accuracy certificate. */
const char *fusion_c_beam_birth_table_kernel_identity(void);
/* Portable immutable-table bytes, not file I/O or host restart. Version1 uses
 * explicit little-endian binary64/integer fields, full sparse knots (including
 * subnormals), exact byte length, kernel identity and an FNV-1a checksum.
 * Source options, grid, controls and sampled validation metadata are retained.
 * Import requires the current kernel identity, validates structure/options and
 * per-knot particle/energy balance, and publishes only a complete table.
 * It does NOT repeat direct quadratures or certify unsampled interpolation;
 * the caller must retain trusted cache provenance and convergence evidence.
 * At most256MiB and the existing node/entry resource bounds are accepted.
 * required/written clear to0 on failure; unpack clears*out first. Short pack
 * buffers are untouched. Caller owns buffer/file I/O; arrays must not overlap.
 * See docs/BEAM_TABLE_BYTES.md for the format and compatibility contract. */
int fusion_c_beam_birth_table_pack_size(const fusion_beam_birth_table_v1 *table,
 size_t *required);
int fusion_c_beam_birth_table_pack(const fusion_beam_birth_table_v1 *table,
 void *buffer,size_t capacity,size_t *written);
int fusion_c_beam_birth_table_unpack(const void *buffer,size_t length,
 fusion_beam_birth_table_v1 **out);
#ifdef __cplusplus
}
#endif
#endif
