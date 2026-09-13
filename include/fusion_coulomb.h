#ifndef FUSION_COULOMB_H
#define FUSION_COULOMB_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fusion_maxwellian_bath_v1 {
    double density_m3;
    double mass_kg;
    double mean_charge_squared;
    double kT_J;
    double coulomb_log;
} fusion_maxwellian_bath_v1;

typedef struct fusion_coulomb_energy_v1 {
    double diffusion_J2_s;       /* half of the energy variance per unit time */
    double mean_energy_rate_J_s; /* signed mean dE/dt of the test particle */
} fusion_coulomb_energy_v1;

/* Classical, nonrelativistic, dilute test particle in a Maxwellian bath.
 * energy_J>=0 is its kinetic energy, mass_kg>0, charge_number is in units e.
 * Bath uses <Z_b^2>, NOT necessarily <Z_b>^2 for a charge-state mixture.
 * coulomb_log is supplied explicitly and must be positive. This coefficient
 * API does not select a screening model, incorporate collective/magnetized
 * scattering, or evolve the Maxwellian bath. The caller must respect the
 * weakly coupled/nonrelativistic assumptions and account for bath heating.
 * Zero density or charge gives zero coefficients. Outputs zero on error.
 * Normalization: NRL Formulary exact relaxation rates (2019 p31), tracing
 * Trubnikov1965; D_E=(m_a v)^2*(nu_parallel*v^2)/2, A_E=-E*nu_energy.
 */
int fusion_c_coulomb_energy(double energy_J, double mass_kg,
    double charge_number, const fusion_maxwellian_bath_v1 *bath,
    fusion_coulomb_energy_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
