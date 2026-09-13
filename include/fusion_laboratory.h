#ifndef FUSION_LABORATORY_H
#define FUSION_LABORATORY_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
/* SI units. Momentum is already oriented in the caller's Cartesian frame. */
typedef struct fusion_particle_four_vector_v1 {
    double mass_kg, kinetic_energy_J, momentum_kg_m_s[3];
} fusion_particle_four_vector_v1;
typedef struct fusion_boost_ledger_v1 {
    double expected_kinetic_energy_J, output_kinetic_energy_J;
    double energy_residual_J, momentum_residual_kg_m_s;
} fusion_boost_ledger_v1;
/* Two massive products in their CM frame. Available kinetic energy >=0.
 * direction[3] is the direction of product 0; norm must agree with unity
 * to 1e-12 and is normalized internally. No nuclear masses/Q/probabilities
 * inferred. output[2] must not overlap inputs. Errors clear outputs. */
int fusion_c_two_body_cm(double mass_a_kg, double mass_b_kg,
    double available_energy_J, const double *direction,
    fusion_particle_four_vector_v1 *output);
/* Active Lorentz boost: a particle initially at rest acquires velocity V.
 * count>=1, |velocity_m_s[3]|<c; massive on-shell input particles required.
 * No orientation randomization. Nonoverlapping arrays input[count],output[count].
 * Output retains masses. Ledger compares sums against a separate boost of
 * summed input four-momentum. Rest energies are excluded from ledger energies.
 * Stateless; errors clear available outputs, never clip invalid states. */
int fusion_c_boost_particles(int count, const double *velocity_m_s,
    const fusion_particle_four_vector_v1 *input,
    fusion_particle_four_vector_v1 *output, fusion_boost_ledger_v1 *ledger);
#ifdef __cplusplus
}
#endif
#endif
