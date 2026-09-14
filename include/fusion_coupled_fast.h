#ifndef FUSION_COUPLED_FAST_H
#define FUSION_COUPLED_FAST_H
#include "fusion_coupled_thermal.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Fast-target channels are independent of the thermal channel selectors.
 * Nuclear model, relative cutoff/order and outgoing-event controls come from
 * thermal_options.birth; only incoming angular quadrature is additional. */
typedef struct fusion_fast_target_options_v1 {
 int channels[5];
 int angular_order;
 double angular_max_exponent;
} fusion_fast_target_options_v1;
/* Additive coupled trial: existing thermal burn, simultaneous old S/T fast
 * fuel vs remaining thermal targets, combined product birth, FP and handoff.
 * First-order split, no fast-fast or within-step re-burning of new products.
 * S/T in one species/energy cell share the same nuclear survival factor.
 * Full beam spectra are evaluated twice with bounded O(cells+reactionedges)
 * storage: once for consistent coefficients and once for actual accepted births.
 * This direct fast-source implementation can be expensive for dense pB spectra.
 * Charged spill rejects; neutrons remain explicit. Rate/debit gates use the
 * existing thermal options for both thermal and fast sources. Fast coefficients
 * freeze at the mixed thermal Ti after thermal burn. No new heat partition.
 * All existing thermal-trial input, clearing and atomic ownership rules apply.
 * fast_options is required even with all channels disabled. All-zero channels
 * preserve the old arithmetic path. No host state is advanced by these calls.
 */
int fusion_c_coupled_fast_trial(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out);
int fusion_c_coupled_fast_table_trial(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,
 const fusion_birth_table_v1 *const *tables,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out);
int fusion_c_coupled_fast_table_trial_effective_charge(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,
 const fusion_birth_table_v1 *const *tables,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out);
/* Diagnosed variants preserve the existing candidate/S-to-T observations.
 * diagnostics is required and clears on any failure. */
int fusion_c_coupled_fast_trial_diagnosed(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out, fusion_handoff_diagnostics_v1 *diagnostics);
int fusion_c_coupled_fast_table_trial_diagnosed(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,
 const fusion_birth_table_v1 *const *tables,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out, fusion_handoff_diagnostics_v1 *diagnostics);
int fusion_c_coupled_fast_table_trial_effective_charge_diagnosed(double dt_s,
 const fusion_coupled_thermal_options_v1 *options,
 const fusion_fast_target_options_v1 *fast_options,
 const fusion_birth_table_v1 *const *tables,int cells,const double *edges_J,
 const double thermal_number_m3[6],double electron_energy_J_m3,double ion_energy_J_m3,
 double electron_density_m3,const double thermal_charge_squared[6],
 int inert_count,const fusion_inert_ion_v1 *inert,const double *coulomb_logs,
 const double *old_s_m3,const double *old_t_m3,const double *external_birth_m3_s,
 const double *escape_s_inv,double trial_thermal_number_m3[6],
 double *trial_s_m3,double *trial_t_m3,fusion_coupled_thermal_v1 *out, fusion_handoff_diagnostics_v1 *diagnostics);
#ifdef __cplusplus
}
#endif
#endif
