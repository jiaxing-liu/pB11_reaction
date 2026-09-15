#ifndef FUSION_THERMAL_PARENT_INTERNAL_H
#define FUSION_THERMAL_PARENT_INTERNAL_H
namespace fusion_detail {
// Shared direct-source parent construction/domain gate. No quadrature or I/O.
int thermal_parent_preflight(int channel,long double relative_energy_J,
 long double cm_energy_J,double mass_a_kg,double mass_b_kg,
 double reaction_q_J,double ground_state_q_J);
}
#endif
