#ifndef FUSION_RATES_H
#define FUSION_RATES_H
#include "fusion_network.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Common Maxwellian kT in joules; output <sigma v> in m^3/s.
 * DT, DD(Tp), DD(He3n): Bosch-Hale 1992, 0.2 <= kT/keV <= 100.
 * DHe3: Bosch-Hale 1992, 0.5 <= kT/keV <= 190.
 * pB: existing Tentori-Belloni integral or FAST model (10..500 keV).
 * Positive finite kT is accepted by pB integral; this is numerical support,
 * not an assertion that nuclear data are validated at arbitrary temperature.
 * pb_method must always be a valid PB11_REACTIVITY_* ID; only pB uses it.
 * One ULP of unit-conversion rounding is admitted at fitted domain edges.
 */
int fusion_c_thermal_reactivity(int channel, double kT_J, int pb_method,
                                double *reactivity_m3_s);

/* Nuclear rest-energy release per event in joules; not product kinetic
 * energy including reactant motion, and not deposited heat. */
int fusion_c_channel_q(int channel, double *q_J);

typedef struct fusion_thermal_rates_v1 {
    double reactivity_m3_s[FUSION_CHANNEL_COUNT];
    double event_rate_m3_s[FUSION_CHANNEL_COUNT];
    fusion_particle_sources_v1 particles;
    double nuclear_power_W_m3[FUSION_CHANNEL_COUNT];
    double total_nuclear_power_W_m3;
} fusion_thermal_rates_v1;

/* Bits 0..4 select channels by fusion_channel_v1. All other bits invalid.
 * Every selected channel's temperature domain is checked even at zero fuel.
 * Disabled channel outputs are zero. All six densities must be finite >=0.
 * DD event rates include 1/2. Product births are NOT thermalization sources.
 * This stateless evaluation performs no integration in time or energy
 * deposition, and does not model non-Maxwellian reactants.
 */
int fusion_c_thermal_rates(double kT_J,
                           const double density_m3[FUSION_SPECIES_COUNT],
                           int channel_mask, int pb_method,
                           fusion_thermal_rates_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
