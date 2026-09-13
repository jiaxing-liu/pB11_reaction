#ifndef FUSION_REACTION_EVENT_H
#define FUSION_REACTION_EVENT_H
#include "fusion_nuclear_data.h"
#ifdef __cplusplus
extern "C" {
#endif
enum fusion_reactant_energy_convention {
 FUSION_REACTANT_CLASSICAL_BUDGET=0, FUSION_REACTANT_ON_SHELL=1
};
typedef struct fusion_reaction_parent_v1 {
 double reactant_kinetic_J[2], classical_kinetic_J[2], on_shell_kinetic_J[2];
 double relative_classical_energy_J, available_cm_energy_J;
 double boost_velocity_m_s[3], expected_product_lab_kinetic_J;
 double classical_minus_on_shell_J, momentum_sum_kg_m_s[3];
 int product_count, convention;
} fusion_reaction_parent_v1;
/* Channel0..4 and reactant momenta[3] in kg m/s, in channel reactant-ID order.
 * Canonical nuclear masses/Q only. Convention0 assigns p²/(2m) to the host
 * kinetic budget; its individual reactants are NOT on-shell four-vectors.
 * Convention1 uses exact on-shell kinetic energy for the SAME momenta.
 * In both cases the aggregate parent has input total energy and momentum;
 * CM product energy and boost are derived from that same aggregate, not from
 * a separately imposed Q+mean relative energy. The classical-minus-on-shell
 * diagnostic quantifies this energy-convention difference; it is not heat.
 * Classical relative energy is retained for existing sigma/rate models.
 * No reactant probability distribution, angular model, source weight or
 * automatic validity threshold for the classical approximation is selected.
 * Outputs clear on errors; all arrays nonoverlapping, finite momenta required.
 */
int fusion_c_reaction_parent(int channel,int convention,const double *pa,
 const double *pb,fusion_reaction_parent_v1 *out);
/* Complete deterministic laboratory event: two products for DD/DT/DHe3,
 * three alphas for pB. output[3] is always required; unused third slot clears.
 * direction[3], norm1 within1e-12, sets product0 direction in parent CM.
 * For pB, q[J] in[0,parent.available_cm_energy_J] is the intermediate excess
 * energy, cosine in[-1,1], azimuth finite[radians]. Exact sequential equal-mass
 * kinematics is used and then rotated/boosted. For two-product channels q,
 * cosine,azimuth must all be0, so irrelevant input is never silently ignored.
 * The auxiliary axis uses cross(z,direction), or cross(x,direction) when
 * |direction.z|>=0.9; azimuth rotates it toward direction cross auxiliary.
 * All PRODUCT four-vectors are on shell, even under classical input budget.
 * Ledger compares product LAB kinetic energy to Q+chosen reactant kinetic
 * energies and momentum to pa+pb. No probability or ghost/spectrum model is
 * implied by supplying one q/direction. Pure marginals cannot reconstruct
 * the correlated laboratory source without additional distribution inputs.
 */
int fusion_c_reaction_lab_event(int channel,int convention,const double *pa,
 const double *pb,const double *direction,double q_J,double cosine,double azimuth,
 fusion_particle_four_vector_v1 *output,fusion_reaction_parent_v1 *parent,
 fusion_boost_ledger_v1 *ledger);
#ifdef __cplusplus
}
#endif
#endif
