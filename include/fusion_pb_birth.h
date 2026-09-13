#ifndef FUSION_PB_BIRTH_H
#define FUSION_PB_BIRTH_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_pb_birth_v1 {
 double mapped_number, mapped_energy_J, below_number, below_energy_J;
 double above_number, above_energy_J, number_residual, energy_residual_J;
 double alpha0_fraction, low_alpha1_fraction, broad_alpha1_fraction;
 double primary_alpha0_energy_J, secondary_min_J, secondary_max_J;
} fusion_pb_birth_v1;
/* Narrow-width alpha0 CM marginal: one primary delta and two uniform
 * secondary boxes from exact sequential relativistic kinematics, using the
 * canonical alpha mass. A in(0,12MeV], q in[0,A], both J. Physical Be8(gs)
 * choice q=91.84keV is caller explicit; other q values are sensitivity or
 * kinematic tests, not measured ground-state energies. J=0 isotropic decay.
 * Per-event source on arithmetic centers; grid spill retained explicitly.
 * No laboratory boost or reactant/branch probabilities are inferred.
 * cells>=1; nonnegative increasing edges[cells+1], birth[cells].
 * All arrays nonoverlapping. Outputs clear on errors; no hidden state. */
int fusion_c_pb_alpha0_grid(double A_J,double q_J,int cells,
 const double *edges_J,double *birth,fusion_pb_birth_v1 *out);
/* Incoherent conditional CM source mixture. f0 and flow are TOTAL EVENT
 * fractions for alpha0 and low-parent l2 alpha1; fbroad=1-f0-flow. Both>=0,
 * sum<=1. These fractions are caller-supplied physics inputs, NOT selected
 * from energy thresholds or claimed to be calibrated by this routine.
 * Broad mode1/3/13, FSCI0/1, cutoff1..10keV, k0..1, finite phase,
 * nq/ncos4..1024. All controls validated even for zero-weight branches.
 * Alpha1 models are nonrelativistic R-matrix marginals; alpha0 uses exact
 * kinematics. Both preserve N=3 and CM kinetic energy A including spill;
 * this does not supply common angular event samples for laboratory boosts.
 * Returned fractions sum to1; alpha0 endpoint diagnostics always supplied.
 */
int fusion_c_pb_cm_source_grid(double A_J,double q_J,double f0,double flow,
 int broad_mode,int fsci_policy,double cutoff_J,double k,double phase,
 int nq,int ncos,int cells,const double *edges_J,double *birth,
 fusion_pb_birth_v1 *out);
#ifdef __cplusplus
}
#endif
#endif
