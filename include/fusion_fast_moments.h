#ifndef FUSION_FAST_MOMENTS_H
#define FUSION_FAST_MOMENTS_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_fast_moments_v1 {
 double number_m3[6], energy_J_m3[6];
 double charge_number_m3, charge_squared_number_m3, pressure_Pa;
} fusion_fast_moments_v1;
/* Explicit moments of the existing isotropic nonrelativistic energy-grid
 * representation. cells1..100000; edges[cells+1] strictly increasing finite
 * nonnegative J. S/T are cell-integrated m^-3, species-major [6*cells].
 * Both components remain kinetic until actual handoff. Arithmetic cell centers
 * match the FP/source-state energy ledger. Nuclear charges come from the
 * canonical nuclear-data API. Pressure=2/3 sum(U) is the SAME isotropic NR
 * closure as host legacy fast pressure, not a relativistic/tensor closure.
 * No Maxwellian fit, thermal-ion heat capacity, rest energy or neutron pressure.
 * All inputs required/finite/nonnegative; inputs/output must not overlap.
 * Stateless, no mutation of input; clears output on failure. Nonzero output
 * underflow or overflow rejects rather than silently removing a weak moment.
 */
int fusion_c_fast_moments(int cells,const double *edges_J,const double *s_m3,
 const double *t_m3,fusion_fast_moments_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
