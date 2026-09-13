#ifndef FUSION_THERMAL_BURN_H
#define FUSION_THERMAL_BURN_H
#include "fusion_network.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Six-species/five-channel order is defined in fusion_network.h.
 * All arrays are fixed size; energies J, densities m^-3, time s.
 */
typedef struct fusion_thermal_burn_v1 {
 double events_m3[FUSION_CHANNEL_COUNT];
 double reactant_removed_m3[FUSION_SPECIES_COUNT];
 double reactant_removed_energy_J_m3[FUSION_SPECIES_COUNT];
 double fast_product_birth_m3[FUSION_SPECIES_COUNT];
 double neutron_birth_m3;
 double number_residual_m3[FUSION_SPECIES_COUNT];
 double energy_residual_J_m3[FUSION_SPECIES_COUNT];
} fusion_thermal_burn_v1;
/* Reaction-only backward-Euler trial with FROZEN reactivities K[5], m^3/s.
 * Events = dt*K*n_a_new*n_b_new/(1+identical), hence each DD branch
 * includes 1/2 exactly once and consumes two D. Secondary DT/DHe3 burning
 * uses T/He3 present in the accepted old thermal pool. Products of this
 * trial are born FAST and do not burn as thermal fuel in the same stage.
 * Disjoint pB pair and coupled D/DD/DT/DHe3 depletion are solved jointly;
 * all species compete for the same thermal inventories. He4 is passive.
 * mean_a/mean_b[5] are REACTION-CONDITIONED reactant energies J/event,
 * e.g. moments/reactivity from fusion_c_thermal_pair_maxwellian_window.
 * They are not generally 1.5*kT. For DD both means debit the same D pool.
 * K=0 requires both means=0. Initial N=0 requires initial U=0. All input
 * values finite/nonnegative; dt>0. Negative remaining U rejects the whole
 * trial (NUMERICAL_FAILURE); caller reduces dt/updates frozen coefficients.
 * No Q or product kinetic energy is created here: birth owner must consume
 * events and removed reactant energies ONCE and use consistent masses/Q.
 * No collisions, external fueling, ash transfer, transport, state or I/O.
 * Amounts, not rates, are returned. Outputs and inputs must not overlap.
 * On error nonnull outputs clear; old state unchanged. Commit/discard is
 * explicit in caller. This solver does not certify supplied nuclear data.
 */
int fusion_c_thermal_burn_trial(double dt_s,
 const double old_number_m3[FUSION_SPECIES_COUNT],
 const double old_energy_J_m3[FUSION_SPECIES_COUNT],
 const double reactivity_m3_s[FUSION_CHANNEL_COUNT],
 const double mean_a_energy_J[FUSION_CHANNEL_COUNT],
 const double mean_b_energy_J[FUSION_CHANNEL_COUNT],
 double trial_number_m3[FUSION_SPECIES_COUNT],
 double trial_energy_J_m3[FUSION_SPECIES_COUNT],fusion_thermal_burn_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
