#ifndef FUSION_BIRTH_TABLE_H
#define FUSION_BIRTH_TABLE_H
#include "fusion_thermal_birth.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_birth_table_v1 fusion_birth_table_v1;
/* Raw source coefficients: K [m^3/s], reacting-energy moments [J m^3/s],
 * spill number [m^3/s] and spill energy [J m^3/s]. Species0..5, neutron6.
 * Products in the accompanying grid are coefficients, NOT rates or amounts.
 * Multiply by n_a*n_b/(1+delta_ab) exactly once outside this API. */
typedef struct fusion_birth_coefficients_v1 {
 double reactivity_m3_s, reactant_energy_moment_J_m3_s[2];
 double below_number_m3_s[7],below_energy_J_m3_s[7];
 double above_number_m3_s[7],above_energy_J_m3_s[7];
} fusion_birth_coefficients_v1;
typedef struct fusion_birth_table_control_v1 {
 /* All tolerances finite in [0,1]. No hidden absolute error floor. */
 double max_rate_error,max_debit_error,max_number_L1,max_energy_L1;
 double max_direct_rate_discrepancy,max_direct_debit_discrepancy;
 /* knots2..100000, evaluations5..1000000, depth0..24. Construction rejects
  * if these limits cannot meet the gates. No partially usable table escapes. */
 int max_knots,max_evaluations,max_depth;
} fusion_birth_table_control_v1;
typedef struct fusion_birth_table_info_v1 {
 double lower_kT_J,upper_kT_J;
 double max_validated_rate_error,max_validated_debit_error;
 double max_validated_number_L1,max_validated_energy_L1;
 double max_sampled_direct_rate_discrepancy,max_sampled_direct_debit_discrepancy;
 int channel,cells,knots,direct_evaluations;
 fusion_thermal_birth_options_v1 source;
 fusion_birth_table_control_v1 control;
} fusion_birth_table_info_v1;
/* Immutable single-channel common-T table. Constructor evaluates the existing
 * direct thermal-birth API at endpoints and geometric quarter/mid/three-quarter
 * points, adaptively bisecting failed intervals in log(kT). Every stored leaf
 * passes the four sampled interpolation gates AND direct rate/debit gates.
 * This is sampled numerical evidence, NOT a supremum error bound or nuclear
 * source-model validation. Independent off-sample refinement is still required.
 *
 * One convex weight interpolates K, both debit moments, EVERY grid coefficient
 * and spill moment; it preserves positivity and linear number/energy identities.
 * Rates/debits/shape are never independently renormalized. Number L1 includes
 * grid and spill-number differences, divided by reference number of that
 * product species; the gate uses the MAXIMUM over all seven species.
 * Energy L1 likewise uses center-weighted grid differences and spill-energy
 * moment differences divided by reference energy of that species, taking the
 * maximum across species. These latter moments
 * do not resolve spill spectral shapes. Zero reference requires exact zero
 * error; otherwise that check fails. No unit-sized tolerance or weak-source cut.
 *
 * cells1..100000, finite increasing nonnegative edges[cells+1]J; finite
 * 0<lower<upper, both in the direct source model domain. To bound memory,
 * 7*cells*(max_knots+3*max_depth+6) must not exceed 50000000 doubles.
 * Create clears *out before validation. All functions catch C++ exceptions.
 * No file I/O, host state, hidden global cache or temperature extrapolation.
 * Destroy once; immutable info/evaluate may be called concurrently. */
int fusion_c_birth_table_create(int channel,double lower_kT_J,double upper_kT_J,
 const fusion_thermal_birth_options_v1 *source,
 const fusion_birth_table_control_v1 *control,int cells,const double *edges_J,
 fusion_birth_table_v1 **out);
void fusion_c_birth_table_destroy(fusion_birth_table_v1 *table);
int fusion_c_birth_table_info(const fusion_birth_table_v1 *table,
 fusion_birth_table_info_v1 *out);
/* Output birth length7*cells. Explicit cells must match table before any
 * evaluation. On errors clear coefficient/info outputs and birth for valid
 * cells only. Edges/model fixed at construction, copied internally. Exact
 * knots reproduce direct coefficients. Strict out-of-range temperatures fail.
 * Arrays must not overlap. This API does not advance a plasma or kinetic state. */
int fusion_c_birth_table_evaluate(const fusion_birth_table_v1 *table,double kT_J,
 int cells,double *birth,fusion_birth_coefficients_v1 *out);
/* Exact constructor-request comparison for safe caller-owned cache reuse.
 * Returns OK with *matches=0 for a different request; OK with 1 only when
 * channel, full grid, both temperature limits, every source option and every
 * construction control match numerically (no tolerance). Null arguments or
 * cells outside1..100000 return INVALID_ARGUMENT, except null matches returns
 * NULL_OUTPUT. NaN request fields cannot match a valid table. No integration,
 * allocation or mutation occurs. */
int fusion_c_birth_table_matches_request(const fusion_birth_table_v1 *table,
 int channel,double lower_kT_J,double upper_kT_J,
 const fusion_thermal_birth_options_v1 *source,
 const fusion_birth_table_control_v1 *control,int cells,const double *edges_J,
 int *matches);
/* Stable process-lifetime ASCII SHA256 string (64 hex characters plus NUL).
 * It is the same generated source/data identity used by beam tables and is
 * not binary provenance, authentication or an accuracy certificate. */
const char *fusion_c_birth_table_kernel_identity(void);
/* Portable immutable-table bytes, not file I/O or host restart. Version1 uses
 * explicit little-endian binary64/integer fields, full dense knots (including
 * subnormals and zero bits), exact byte length, kernel identity and an FNV-1a
 * checksum. Source options, grid, controls and sampled validation metadata are
 * retained. Import requires the current kernel identity, validates structure,
 * source/control domains and per-knot particle/energy balance, and publishes
 * only a complete table. It does not repeat direct quadratures or certify
 * unsampled interpolation. At most512MiB and the original constructor
 * resource bounds are accepted. required/written clear to0 on failure;
 * unpack clears*out first. Short pack buffers are untouched. Caller owns
 * buffer/file I/O; arrays must not overlap. */
int fusion_c_birth_table_pack_size(const fusion_birth_table_v1 *table,
 size_t *required);
int fusion_c_birth_table_pack(const fusion_birth_table_v1 *table,
 void *buffer,size_t capacity,size_t *written);
int fusion_c_birth_table_unpack(const void *buffer,size_t length,
 fusion_birth_table_v1 **out);
#ifdef __cplusplus
}
#endif
#endif
