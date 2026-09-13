#ifndef FUSION_COUPLED_THERMAL_H
#define FUSION_COUPLED_THERMAL_H
#include "fusion_thermal_birth.h"
#include "fusion_source_state.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_inert_ion_v1 {
 double density_m3, mass_kg, mean_charge_squared;
} fusion_inert_ion_v1;
typedef struct fusion_coupled_thermal_options_v1 {
 fusion_thermal_birth_options_v1 birth;
 double max_source_rate_error, max_source_debit_error;
 double handoff_max_L1, handoff_max_mean_error;
 int channels[5], handoff_enabled;
} fusion_coupled_thermal_options_v1;
typedef struct fusion_coupled_thermal_v1 {
 fusion_source_ledger_v1 ledger;
 /* Extra, non-network ion baths have their own heat account, summed over
  * inert species for each fast species. Never label these as H or e heat.
  * ledger ALONE is insufficient for source_state_v1_stage if these are nonzero.
  * This complete trial needs an atomic thermal/kinetic owner before host use. */
 double inert_ion_heat_J_m3[6];
 double electron_energy_J_m3, ion_energy_J_m3;
 double particle_residual_m3[6], energy_residual_J_m3;
 double max_source_rate_discrepancy, max_source_debit_discrepancy;
 double handoff_L1[6], handoff_mean_error[6];
 int handoff_projected[6];
} fusion_coupled_thermal_v1;
/* Stateless local trial: common thermal-ion temperature, separate electrons.
 * Thermal ion energy includes the six network species AND explicit inert ions.
 * Inert densities do not evolve here. Electron density is supplied/frozen for
 * the trial; host owns charge balance, ionization and boundary electron fluxes.
 * Bare nuclear rates, nuclear charges for fast ions, supplied <Z_i^2> for baths.
 * No implicit ionization energy, geometry, transport, radiation or fueling.
 *
 * Split: freeze nuclear coefficients at old Ti; simultaneous backward-Euler
 * fuel depletion; integrate birth packets; freeze Coulomb baths at depleted
 * common Ti; implicit two-component FP; feed signed collision heat back to
 * Ue/Ui; optional measured handoff into a self-consistent mixed ion pool.
 * First-order splitting, NOT a converged implicit bath iteration. Time/source/
 * grid convergence and dilute NR collision validity remain caller obligations.
 *
 * cells 1..100000; inert_count 0..32. Grid in J; populations cell-integrated
 * m^-3, arrays old/trial_s,t and external_birth,escape are [6*cells], species
 * major. external_birth m^-3 s^-1 enters S; escape s^-1 is a physical caller
 * policy carrying particle/energy explicitly. Coulomb logs are positive,
 * [6*(7+inert_count)] fast-major, baths electron, network0..5, then inert ions.
 * thermal_charge_squared[6] is <Z^2>, not <Z> squared for a charge mixture.
 * birth.pb_low selects only pB; other channels use their native low model.
 * Every array required except inert when count=0. Arrays must not overlap.
 *
 * Charged nuclear grid spill rejects the entire trial; enlarge the grid.
 * Neutrons are returned as separate full birth number/energy, including their
 * grid spills. They are not an automatic local heat source. Source rate/debit
 * discrepancies must satisfy explicit options before burn; success still does
 * not certify pB angular/source-proxy convergence or differential data.
 *
 * Handoff processes species0..5 sequentially. For each complete T candidate,
 * targetTi=2*(Ui+Ucandidate)/(3*(Nthermal+Ninert+Ncandidate)); the measured
 * shape/energy gates apply at that target. On pass, fluid energy and signed
 * correction sum to removed kinetic energy and final common Ti equals target.
 * The correction is accounted to that species' thermal-ion ledger column.
 * Species-order/splitting sensitivity must be checked. S->T is internal, not ash.
 *
 * Output amounts are per step. Caller must accept/discard thermal N/U and both
 * kinetic arrays atomically; repeated trials have no accumulated state. All
 * Residual checks include explicit double-roundoff allowance from the actual
 * stored thermal and kinetic inventories, not a unit-sized absolute floor.
 * outputs clear on failure for valid dimensions; accepted inputs unchanged.
 */
int fusion_c_coupled_thermal_trial(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
