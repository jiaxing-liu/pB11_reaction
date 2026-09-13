#ifndef FUSION_ALPHA_SPECTRUM_H
#define FUSION_ALPHA_SPECTRUM_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_alpha_spectrum_v1 {
 double mapped_number, mapped_energy_J, below_number, below_energy_J;
 double above_number, above_energy_J, number_residual, energy_residual_J;
 double normalization_J2, l1_normalization_J2, l3_normalization_J2;
 int quadrature_events, pruned_events;
} fusion_alpha_spectrum_v1;
/* Normalized, rotationally averaged three-alpha CM births PER REACTION.
 * mode=1,2,3 selects one primary l; mode=13 selects coherent l1/l3 mixture.
 * In mode13, k in[0,1] is the fraction between separately unit-normalized
 * amplitude bases BEFORE interference; phase in radians. Interference is
 * followed by overall normalization. This k is not a final branching ratio.
 * A is total CM kinetic energy in J, 0<A<=12 MeV; normalization must be
 * finite and positive. This is a numerical, not nuclear-data, interval.
 * cutoff in [1,10] keV is explicit numerical policy
 * from fusion_c_alpha_amplitudes_cutoff; errors do not drop whole events.
 * nq/ncos in[4,1024]; cells>=1. edges[cells+1] inJ, birth[cells] number/event.
 * Gauss-Legendre q=A*sin(theta)^2 with theta in[0,pi/2], and cos in[-1,1].
 * Arithmetic-center projection retains explicit below/above number and energy
 * per event. The total must close to3 particles and A, INCLUDING spill.
 * Pure source shape: no nuclear cross-section, resonance mixture, alpha0
 * branching, laboratory boost or device deposition is selected here.
 * Arrays must not overlap. Outputs clear on error; no hidden state.
 */
int fusion_c_alpha_spectrum_grid(int mode,double available_energy_J,
 double cutoff_J,double l1_fraction,double relative_phase,int nq,int ncos,
 int cells,const double *edges_J,double *birth_per_event,fusion_alpha_spectrum_v1 *out);
/* Same normalized CM-grid contract with explicit FSCI policy0(NONE) or
 * 1(Refsgaard2018 Model-II at16fm). Mixture basis normalization is recomputed
 * WITH the selected correction before coherent interference. Still no
 * incident-resonance, alpha0 branching, lab boost or deposition selection.
 * The no-policy entry point above is permanently equivalent to policy0. */
int fusion_c_alpha_spectrum_model_grid(int mode,int fsci_policy,
 double available_energy_J,double cutoff_J,double l1_fraction,
 double relative_phase,int nq,int ncos,int cells,const double *edges_J,
 double *birth_per_event,fusion_alpha_spectrum_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
