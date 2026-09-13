#ifndef FUSION_PRODUCTS_H
#define FUSION_PRODUCTS_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fusion_birth_mapping_v1 {
    double input_number_m3_s, input_energy_W_m3;
    double mapped_number_m3_s, mapped_energy_W_m3;
    double below_number_m3_s, below_energy_W_m3;
    double above_number_m3_s, above_energy_W_m3;
    double number_residual_m3_s, energy_residual_W_m3;
} fusion_birth_mapping_v1;

/* Project supplied birth packets onto the arithmetic-center energy grid used
 * by fusion_c_energy_fp_trial. Each packet has energy J and nonnegative
 * particle birth rate m^-3 s^-1; it need not be a whole event or equal weight.
 * Adjacent-center linear weights conserve both number and kinetic energy.
 * Out-of-center-hull packets are returned in explicit below/above ledgers,
 * NOT clipped, deleted or automatically thermalized/escaped. OK means the
 * mapping and those ledgers were computed; caller must handle all spill.
 * cells>=1, packets>=0, increasing nonnegative edges[cells+1]. Empty packet
 * arrays may be NULL. Nonoverlapping input/output arrays; no hidden state.
 * Finite inputs, nonnegative rates/energies; outputs zero on error.
 */
int fusion_c_map_birth_packets(int cells, const double *edges_J, int packets,
    const double *packet_energy_J, const double *packet_rate_m3_s,
    double *cell_birth_m3_s, fusion_birth_mapping_v1 *out);

typedef struct fusion_three_body_cm_v1 {
    double kinetic_energy_J[3];
    double momentum_x_kg_m_s[3];
    double momentum_z_kg_m_s[3];
    double energy_residual_J;
    double momentum_residual_kg_m_s;
} fusion_three_body_cm_v1;

/* Exact relativistic SEQUENTIAL KINEMATICS for three equal-mass products.
 * Parent CM available kinetic energy is A>=0. Intermediate rest energy is
 * 2*m*c^2+q, with 0<=q<=A. It decays to products 1 and 2; product 0 is the
 * primary emitted along +z. cos_theta is the secondary direction relative to
 * the intermediate's motion (-z) in its rest frame; azimuth is chosen x-z.
 * Returns all three CM kinetic energies and momenta, with y=0.
 * This function selects NO q distribution, branching or angular probability.
 * It is NOT a pB spectrum model; caller must supply validated event weights
 * including any required coherent permutations and angular correlations.
 * Mass kg>0, energies J, -1<=cos_theta<=1. No nuclear Q inferred internally.
 * Outputs zero on errors, inputs unchanged, no random generator or state.
 */
int fusion_c_three_equal_sequential_cm(double product_mass_kg,
    double available_energy_J, double intermediate_relative_energy_J,
    double cos_theta, fusion_three_body_cm_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
