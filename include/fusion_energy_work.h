#ifndef FUSION_ENERGY_WORK_H
#define FUSION_ENERGY_WORK_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_energy_work_ledger_v1 {
 double initial_number_m3,final_number_m3;
 double initial_energy_J_m3,final_energy_J_m3;
 double lower_number_m3,lower_energy_J_m3;
 double upper_number_m3,upper_energy_J_m3;
 double work_on_particles_J_m3;
 double particle_balance_error_m3,energy_balance_error_J_m3;
} fusion_energy_work_ledger_v1;
/* Stateless backward-Euler upwind trial for d_t g+d_E(-C E g)=0.
 * C [s^-1] is signed: positive cools, negative heats. Isotropic nonrelativistic
 * adiabatic closure uses C=(2/3)div(u); caller owns model and geometry.
 * edges[cells+1] increasing nonnegative J; population[cells] is cell-integrated
 * m^-3, NOT density/J. Energy is represented by the arithmetic cell midpoint.
 * dt>0 seconds. Cells 1..1000000. All inputs finite, populations nonnegative.
 * Boundary inflow is zero; outgoing particles are returned in distinct lower/
 * upper ledgers with boundary-edge energy, NOT silently thermalized. At E=0
 * the drift flux is zero. No volume dilution or physical-space transport here.
 * face_amount[cells+1] is integrated signed flux toward higher energy [m^-3].
 * Discrete work: Unew-Uold+Elower+Eupper. It includes discretization error;
 * grid/time refinement is required to compare with continuum -C*U power.
 * Inputs/outputs must not overlap. Valid-dimension outputs clear on any error.
 * Finite trace values round by IEEE; measured conversion errors are included
 * in extended-precision conservation checks, without a population floor.
 */
int fusion_c_energy_work_trial(int cells,double dt_s,double compression_s_inv,
 const double*edges_J,const double*old_number_m3,double*trial_number_m3,
 double*face_amount_m3,fusion_energy_work_ledger_v1*ledger);
#ifdef __cplusplus
}
#endif
#endif
