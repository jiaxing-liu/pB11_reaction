#ifndef FUSION_ALPHA_AMPLITUDES_H
#define FUSION_ALPHA_AMPLITUDES_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_alpha_amplitudes_v1 {
    double unsym_real[5], unsym_imag[5];
    double sym_real[5], sym_imag[5];
    double phase_space_J;
} fusion_alpha_amplitudes_v1;
/* Sequential nonrelativistic spin-J=2 -> Be8(J=2) -> 3alpha amplitudes.
 * primary_l=1,2,3; secondary_l=2. Supplied A and q in SI J satisfy 0<q<A.
 * cos_theta is measured relative to intermediate recoil in its rest frame.
 * Nonrelativistic momenta are constructed internally with m_alpha c²=
 * 3727.3794118MeV; this is the convention of the R-matrix amplitude.
 * Output spin projections m=-2..2; 'sym' is coherent cyclic permutation sum.
 * Primary reduced width squared fixed at1MeV; intermediate E0=3.129MeV,
 * gamma²=1.075MeV, radius4.5fm and B=S(E0); primary radius5.1fm.
 * Phase-space factor sqrt(q*(A-q)) multiplies squared amplitudes with dq dcos.
 * Overall normalization is arbitrary: these are NOT normalized probabilities.
 * Angular-momentum algebra follows the explicit coherent amplitude, not a
 * separate empirical angle fit. No l mixture/branching/incident-energy model
 * or final-state Coulomb correction selected. Nuclear Coulomb domain applies
 * to every primary/secondary permutation; errors never omit permutations.
 * Outputs zero on error; no hidden state or random sampling. */
int fusion_c_alpha_amplitudes(int primary_l,double available_energy_J,
    double intermediate_energy_J,double cos_theta,fusion_alpha_amplitudes_v1 *out);
/* Explicit small-energy numerical cutoff variant. Same amplitude convention,
 * but each permutation with primary or pair energy < cutoff_J is set to zero.
 * Other permutations remain. cutoff_J in[0.001,0.01]MeV; q endpoints allowed.
 * Returns number of suppressed permutations (0..3), never hides suppression.
 * This approximation requires cutoff-sensitivity/reference validation for the
 * intended spectrum; it is not a physical branching or nuclear-data cutoff.
 * Both outputs required; cleared on error. No spectrum normalization occurs. */
int fusion_c_alpha_amplitudes_cutoff(int primary_l,double available_energy_J,
    double intermediate_energy_J,double cos_theta,double cutoff_J,
    fusion_alpha_amplitudes_v1 *out,int *pruned_permutations);
#ifdef __cplusplus
}
#endif
#endif
