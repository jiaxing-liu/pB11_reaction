#ifndef FUSION_NUCLEAR_COULOMB_H
#define FUSION_NUCLEAR_COULOMB_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_ALPHA_BE8_L1=0, FUSION_ALPHA_BE8_L2=1, FUSION_ALPHA_BE8_L3=2,
       FUSION_ALPHA_ALPHA_L0=3, FUSION_ALPHA_ALPHA_L2=4, FUSION_NUCLEAR_CHANNEL_COUNT=5 };
/* Fixed nuclear channel conventions, unrelated to plasma Coulomb slowing:
 * IDs0,1,2: alpha-Be8 l=1,2,3, radius5.1fm, reduced mass2*m_alpha/3.
 * IDs3,4: alpha-alpha l=0,2, radius4.5fm, reduced massm_alpha/2.
 * m_alpha*c^2=3727.3794118MeV; alpha^-1=137.035999084;
 * hbar*c=197.3269804MeVfm. No hidden radius or boundary selection. */
typedef struct fusion_nuclear_coulomb_v1 {
    double log_penetrability, shift, phase_real, phase_imag, rho;
} fusion_nuclear_coulomb_v1;
/* Energy SI J, supported relative energy 0.001..12MeV inclusive.
 * Returns log(P_l), S_l, exp(i[omega_l-phi_l]), rho=k*a.
 * Out-of-domain energies fail, NOT silently zero/extrapolated. This is a
 * generated numerical approximation, not a nuclear cross-section model.
 * Stateless; output zero on error. Validation/generator are versioned. */
int fusion_c_nuclear_coulomb(int channel, double relative_energy_J,
    fusion_nuclear_coulomb_v1 *out);
/* Same channel IDs, masses and energy domain as above, but ALL channel
 * radii are explicitly16fm. Separately generated/validated tables; used in
 * the Refsgaard2018 Model-II FSCI factor. Does not replace ordinary radii.
 * Same errors/clearing contract. */
int fusion_c_nuclear_coulomb_radius16(int channel,double relative_energy_J,
 fusion_nuclear_coulomb_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
