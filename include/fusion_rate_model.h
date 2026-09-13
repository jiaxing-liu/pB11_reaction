#ifndef FUSION_RATE_MODEL_H
#define FUSION_RATE_MODEL_H
#include "fusion_beam.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit model assumptions, not confidence bounds or measured cross sections.
 * All energies J, masses kg, sigma m^2. No density or identical-pair factor.
 * ENDPOINT_S: constant endpoint astrophysical S below/above fit domain with
 * the channel's fixed Gamow exponent. HIGH_FLAT changes only the high tail
 * to constant sigma. pB fit reaches zero: no extra low continuation there.
 * pB low variants replace only E<=400 keV; all higher fit pieces stay TB2023.
 * NS_LOW is NS as transcribed by TB2023, NOT the complete NS rate model.
 * C0 variations +/-12 MeV barn are one-parameter sensitivity, not covariance.
 * For non-pB channels pb_low must be TB; inappropriate options are rejected.
 */
enum fusion_continuation_policy {
    FUSION_ENDPOINT_S = 1, FUSION_HIGH_FLAT = 2
};
enum fusion_pb_low_policy {
    FUSION_PB_LOW_TB = 0, FUSION_PB_LOW_NS = 1,
    FUSION_PB_LOW_C0_MINUS12 = 2, FUSION_PB_LOW_C0_PLUS12 = 3
};
/* Each segment uses the existing moment layout; "resolved" means integrated
 * in that segment under the EXPLICIT model, not experimentally certified.
 * Fit owns both endpoints and the E=0 limit. Below/above are strictly outside.
 * In total, domain_incomplete=0 means this model defines all positive energies;
 * it does NOT certify nuclear-data accuracy. Total unresolved diagnostics=0.
 * Segment unresolved diagnostics describe the complement of that segment.
 * Same nonrelativistic Maxwellian and numerical speed-ratio domain as window
 * APIs. Gaussian tails below representable precision may round to zero.
 * No hidden state, clipping, probability renormalization, or silent fallback.
 * All non-NULL outputs cleared on failure; inputs unchanged.
 */
typedef struct fusion_rate_model_v1 {
    fusion_beam_window_v1 total, fit, below, above;
} fusion_rate_model_v1;
int fusion_c_cross_section_model(int channel, int continuation, int pb_low,
    double relative_energy_J, double *cross_section_m2);
int fusion_c_beam_maxwellian_model(int channel, int continuation, int pb_low,
    double projectile_mass_kg, double target_mass_kg,
    double projectile_energy_J, double target_kT_J, fusion_rate_model_v1 *out);
int fusion_c_thermal_pair_maxwellian_model(int channel, int continuation, int pb_low,
    double reactant_a_mass_kg, double reactant_b_mass_kg,
    double reactant_a_kT_J, double reactant_b_kT_J, fusion_rate_model_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
