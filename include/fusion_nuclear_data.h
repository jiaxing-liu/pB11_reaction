#ifndef FUSION_NUCLEAR_DATA_H
#define FUSION_NUCLEAR_DATA_H
#include "fusion_network.h"
#include "fusion_laboratory.h"
#ifdef __cplusplus
extern "C" {
#endif
/* IDs 0..5 retain fusion_species_v1; these two are data IDs, not additions
 * to the six thermal network slots. */
enum fusion_mass_data_id { FUSION_MASS_NEUTRON=6, FUSION_MASS_ELECTRON=7 };
typedef struct fusion_nuclear_mass_v1 {
 double mass_kg, rest_energy_J, reference_mass_u;
 double known_uncertainty_scale_kg;
 int nuclear_charge, mass_number;
} fusion_nuclear_mass_v1;
typedef struct fusion_nuclear_channel_v1 {
 double q_J, known_uncertainty_scale_J;
 int reactant_ids[2], product_ids[3], product_count;
} fusion_nuclear_channel_v1;
/* Fixed new-model mass set: CODATA2022 relative nuclear masses and common
 * atomic mass constant; B11 from AME2020 ATOMIC mass minus5electron masses
 * plus NIST ASD2024 total B electron binding (670.9838405eV)/c^2.
 * B binding is an element-table approximation, not isotope-specific data.
 * Known uncertainty scale is a linear propagation of published standard
 * uncertainties (conservative without covariance), NOT a confidence interval
 * or bound on this unquantified isotope/model systematic. See NUCLEAR_DATA.md.
 * Charge is bare nuclear Z (electron:-1), not coronal mean ionic charge.
 * kg is canonical rounded-double mass; rest J is calculated from that mass.
 * Output clears on error; invalid IDs return OUT_OF_RANGE; no state/I/O.
 */
int fusion_c_nuclear_mass(int particle_id,fusion_nuclear_mass_v1 *out);
/* Q is computed from the SAME returned kg masses with exact SI c, not an
 * independently rounded Q table. Missing third product ID is -1 on success.
 * pB product_count=3; other channels2; neutrons are explicit products.
 * Old fusion_c_channel_q / fusion_c_thermal_rates preserve their previous
 * compatibility convention, including pB8.68MeV; do not mix conventions.
 */
int fusion_c_nuclear_channel(int channel,fusion_nuclear_channel_v1 *out);
/* Deterministic two-product birth composition for DD branches, DT, DHe3:
 * available product CM kinetic energy = reactant_CM_kinetic_J + new-model Q.
 * Input CM kinetic energy[J]>=0; direction[3] norm1. Uses exact two-body
 * kinematics and masses above. No angular probability model is implied.
 * pB is OUT_OF_RANGE here; it requires the separate three-alpha source.
 * Output[2] and inputs must not overlap. Outputs clear on error. */
int fusion_c_nuclear_two_body_cm(int channel,double reactant_CM_kinetic_J,
 const double *direction,fusion_particle_four_vector_v1 *output);
#ifdef __cplusplus
}
#endif
#endif
