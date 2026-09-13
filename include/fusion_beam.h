#ifndef FUSION_BEAM_H
#define FUSION_BEAM_H
#include "fusion_cross_sections.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fusion_beam_window_v1 {
    double resolved_reactivity_m3_s;
    double projectile_energy_reactivity_J_m3_s;
    double target_energy_reactivity_J_m3_s;
    double relative_energy_reactivity_J_m3_s;
    double cm_energy_reactivity_J_m3_s;
    double resolved_pair_probability;
    double unresolved_pair_probability;
    double unresolved_relative_speed_m_s;
    double quadrature_error_m3_s;
    double energy_identity_error_J_m3_s;
    int domain_incomplete;
} fusion_beam_window_v1;

/* Monoenergetic projectile against a stationary isotropic Maxwellian target.
 * Classical nonrelativistic kinematics. Masses are explicit caller inputs in
 * kg, not inferred from the channel; the caller must supply consistent masses.
 * Projectile kinetic energy and target kT are J, both >=0; kT=0 is an exact
 * stationary-target limit. No density factors or identical-pair 1/2 included.
 *
 * Only the documented cross-section window is integrated, without extension
 * or renormalizing the retained relative-velocity probability. OK means the
 * WINDOW result was computed, not that it is the full physical reactivity.
 * domain_incomplete=1 for any nonzero-temperature target (infinite tail), or
 * for a cold target outside the cross-section window. The zero-energy
 * cross-section exception follows fusion_c_cross_section.
 *
 * The unresolved probability and mean relative speed are diagnostics, NOT a
 * reaction-rate error estimate. If an independently justified upper bound S
 * on cross section outside the window exists, S*unresolved_relative_speed is
 * the corresponding missing-rate bound in exact arithmetic. The reported
 * moment is a numerical estimate, not a certified interval; quadrature and
 * rounding uncertainty also apply. This API does not invent the physical S.
 * Gaussian tails below representable numerical precision may round to zero;
 * domain_incomplete still reports the mathematically incomplete domain.
 * quadrature_error is only the resolved integral's numerical error estimate.
 *
 * Energy quantities are integrals of sigma*w times the indicated reacting
 * particle/pair energy (J m^3/s), not ordinary Maxwellian mean energies.
 * Their division by resolved reactivity gives reaction-conditioned energies
 * when that reactivity is nonzero. Projectile+target=CM+relative energy.
 * These are source primitives, not a complete product/ash closure.
 * Inputs unchanged, output cleared on failure; no hidden state or I/O.
 */
int fusion_c_beam_maxwellian_window(int channel, double projectile_mass_kg,
    double target_mass_kg, double projectile_energy_J, double target_kT_J,
    fusion_beam_window_v1 *out);
/* Two independent, zero-drift Maxwellian reactants, possibly unequal kT.
 * The same output basis is used: "projectile" means reactant a, "target" b.
 * The relative distribution has Trel=(mb*Ta+ma*Tb)/(ma+mb), but the two
 * reactant energy moments retain their CM/relative correlations. They are
 * NOT obtained by assigning 3Ta/2 and 3Tb/2 to every reacting pair.
 * Same finite-window, density-free and identical-pair counting contract.
 */
int fusion_c_thermal_pair_maxwellian_window(int channel,
    double reactant_a_mass_kg, double reactant_b_mass_kg,
    double reactant_a_kT_J, double reactant_b_kT_J,
    fusion_beam_window_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
