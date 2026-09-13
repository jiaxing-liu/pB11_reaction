#ifndef FUSION_PB_POPULATION_H
#define FUSION_PB_POPULATION_H
#include "fusion_beam.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_pb_population_v1 {
 double total_cross_section_m2, narrow_fit_cross_section_m2;
 double narrow_fit_fraction, alpha0_peak_fraction;
 double narrow_remainder_fraction, other_remainder_fraction;
 double continuum_peak_fraction, continuum_extrapolated_fraction;
 double equivalent_proton_lab_energy_J;
} fusion_pb_population_v1;
/* EXPLICIT effective population model, not an amplitude decomposition or
 * unfolded channel fit. E is relative CM energy[J]>=0; continuation1/2 and
 * pb_low0..3 are fusion_rate_model choices. The positive narrow S-factor term
 * at E<=400keV sets r; r=0 above that FIT piece boundary. This does not claim
 * that the physical resonance vanishes there. E=0 uses the limiting S ratio.
 * caller narrow_peak in[0,1] is isolated narrow-parent ground-state PEAK
 * fraction (Laursen5.1% excludes its ghost); continuum_scale in[0,2] multiplies
 * the effective continuum peak table. The table uses Taskaev2024 ratios from
 * E_lab_target_avg247..2186keV, linear interpolation in equivalent stationary
 * proton lab energy, endpoint hold outside. This is an explicitly approximate
 * use of target-averaged data, NOT deconvolution. The75/134keV rows are excluded.
 * E_lab_equiv=E*(mp+mB)/mB uses canonical masses. continuum_extrapolated_fraction
 * is1-r outside the tabulated range,0 inside, not an uncertainty bound.
 * Outputs: alpha0_peak=r*narrow_peak+(1-r)*table_scaled;
 * narrow_remainder=r*(1-narrow_peak); other=(1-r)*(1-table_scaled).
 * The remainders are NOT declared pure alpha1 or specific angular states.
 * No spectrum/ghost/laboratory distribution is chosen. Outputs clear errors.
 */
int fusion_c_pb_population(int continuation,int pb_low,double E_J,
 double narrow_peak,double continuum_scale,fusion_pb_population_v1 *out);
typedef struct fusion_pb_population_rates_v1 {
 fusion_beam_window_v1 total, alpha0_peak, narrow_remainder, other_remainder;
 fusion_beam_window_v1 continuum_extrapolated;
} fusion_pb_population_rates_v1;
/* Source-weighted rates AND all selected reactant/CM/relative energy moments.
 * Same nonrelativistic distributions and explicit model/domain contracts as
 * fusion_*_maxwellian_model, with pB channel fixed. No density factor or
 * identical-pair factor. Caller supplies canonical proton/B masses in the
 * corresponding roles; projectile may be either p or B. Mass equality with
 * the canonical pair is checked, preventing table-energy role confusion.
 * Probability diagnostics also use the declared population weight in the
 * kinematic pair distribution (they are NOT normalized reaction fractions).
 * total equals the existing full model. Three populations sum to total;
 * continuum_extrapolated is a diagnostic subset, NOT a fourth population.
 * These do not yet compute product birth spectra or reactant angular samples.
 */
int fusion_c_pb_thermal_population_rates(int continuation,int pb_low,
 double mass_a_kg,double mass_b_kg,double kTa_J,double kTb_J,
 double narrow_peak,double continuum_scale,fusion_pb_population_rates_v1 *out);
int fusion_c_pb_beam_population_rates(int continuation,int pb_low,
 double mass_a_kg,double mass_b_kg,double projectile_energy_J,double target_kT_J,
 double narrow_peak,double continuum_scale,fusion_pb_population_rates_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
