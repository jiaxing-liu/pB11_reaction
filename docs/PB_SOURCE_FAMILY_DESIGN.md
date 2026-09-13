# Physical pB birth-source composition: selected model family

This is the next source-composition design, following source-state3945bd9.
It does not claim that the full host source is already implemented.

## Decisions required by the evidence

A single broad-resonance l1/l3 spectrum must not be applied to all pB events.
The low 16.11-MeV 2+ parent and broad 2- region have different allowed primary
partial waves, and alpha0 (Be8 ground-state) events require a separate source.
The existing scalar total fusion cross section supplies event rates and must
not be multiplied by the parent's approximately 50% fusion-vs-elastic decay
probability a second time.

The physically motivated source family will retain:

1. Low-parent alpha1: coherent primary-l2, secondary-l2 amplitudes. Add the
   Refsgaard2018 Model-II final-state Coulomb factor, separately within every
   permutation before coherent summation. Use an explicit fixed-r16fm policy
   as a published-radius model, alongside the existing no-extra-FSCI policy.
   This is a phenomenological nuclear breakup model, not a full three-body
   scattering calculation. It does not resolve the earlier inconsistent
   standalone Laursen angular formula; the explicit CG amplitude stays intact.
2. Broad-region alpha1: explicit l1/l3 unit-basis mixture. The published k=.76
   is not relabeled a certified coefficient in our normalization convention.
   The previously quantified .76 versus raw-reinterpreted .81916 and pure-wave
   alternatives remain separate model sensitivities until source-weighted
   impacts are quantified. FSCI policy is explicit, never silently enabled.
3. Alpha0: sequential primary alpha plus Be8(gs) -> two alphas; secondary
   s-wave has no directional memory. Use the measured Be8 relative energy
   91.84(4) keV (Refsgaard Table2), distinguish it from independently rounded
   total nuclear Q, and retain three-particle kinematic conservation.
4. Event-energy/branch weighting: use the existing rate-model energy weights,
   distinguish the narrow entrance-resonance contribution from the remaining
   continuum, and retain alpha0/alpha1 branch evidence with its target-energy
   averaging assumptions. Taskaev table ratios must not masquerade as an
   unfolded resonance decomposition. This weighting is still to be calibrated
   and validated; no off-domain branch is silently set to zero.
5. Laboratory source: account for reactant-selected center-of-mass motion,
   not just add Q to a mean CM spectrum. For common-temperature Maxwellian
   reactants the center-of-mass and relative velocities are independent;
   unequal-temperature and beam cases need their actual correlations and
   declared angular model. A marginal energy grid alone cannot reconstruct
   every angular correlation. This remains required composition work.

## FSCI formula and implementation boundary

Refsgaard et al., Physics Letters B779 (2018)414–419, printed p415 Eq4:

`P_l1/rho1 -> (P_l1/rho1) * C`,

`C = (rho1_tilde/P_l1_tilde) * (P_l2_tilde(E12)/rho12_tilde) * (P_l2_tilde(E13)/rho13_tilde)`.

All tilde quantities use the same enlarged separation radius; the original
primary radius5.1fm and secondary radius4.5fm remain in the ordinary amplitude.
Thus multiply each unsymmetrized amplitude by sqrt(C), not C, and only then
sum permutations. Keep the existing Coulomb phases. For the alpha1 model,
l2=2 in both additional pair factors. Use logarithms to avoid multiplying an
underflowed ordinary amplitude by an overflowing correction.

The follow-up paper explicitly identifies a calculation error behind the
older reported10fm value; it describes approximately15fm as reasonable for
the short-lived excited Be8 decay and uses16fm in its tabulated Model-II
study. A fixed16fm choice is a declared model, not an accurately measured
universal separation. A future energy-dependent lifetime Model III is a
separate model and must not be inferred from a fixed-radius result.

The numerical Coulomb tables must be generated and checked at the new radius;
changing rho alone in old fixed-radius tables is invalid. Cache signatures
must contain charge product, reduced mass, radius and partial waves. Ordinary
channel tables and compatibility APIs remain unchanged. Any numerical cutoff
must apply explicitly to the new pair factors too, with suppression counters;
never replace failed Coulomb evaluations by physical zero without policy.

## Acceptance work

- Check new-radius Coulomb values independently with mpmath, including
  Wronskian, interpolation residual and metadata checks.
- Check amplitude correction against direct finite-radius Coulomb references,
  per-permutation coherence, no-FSCI compatibility and cutoff behavior.
- Check normalized source number/energy and spectral-bin convergence for
  the actual source grids; global moments alone were inadequate previously.
- Quantify corrected/uncorrected shape differences and compare the relevant
  published low-parent model/data with its selection definitions intact.
- Compose physical branch weighting and lab births with thermal/fast burn,
  collisions and measured handoff in the accepted-state manager, then validate
  evolving backgrounds. This work is not complete after the table/API step.

Primary local evidence: /tmp/m3-spectrum-sources/refsgaard-2018.pdf, pp415–417;
Laursen2016 and Kuhlwein2021 PDFs in the same directory. The root visually
checked Refsgaard Eq4; the completed bounded audit is in PB_FINAL_STATE_COULOMB_EVIDENCE.md.

## Alpha0 CM marginal: analytic implementation path

For the narrow Be8(gs) branch, q_gs=91.84keV and mass M*=2m_alpha+q_gs/c^2.
First use exact two-body kinematics for parent rest energy
3m_alpha*c^2+A -> alpha + Be8*. The primary has fixed kinetic energy K0;
the intermediate has gamma and beta from its recoil. In its rest frame each
secondary has total energy epsilon*=m_alpha*c^2+q_gs/2 and momentum p*.
Because the ground state has J=0, its secondary direction cosine is uniform.
The CM kinetic energy of either secondary is therefore uniform between

`K_minus = gamma*(epsilon* - beta*c*p*) - m_alpha*c^2`
and
`K_plus  = gamma*(epsilon* + beta*c*p*) - m_alpha*c^2`.

The per-event marginal is one primary delta plus two identical uniform boxes.
The three particles are correlated; this expression only specifies the
one-particle energy marginal. It gives exactly N=3 and
`K0 + K_minus + K_plus = A`. Integrate the piecewise-linear grid projection
analytically over the box, with explicit number/energy outside the center
hull. This avoids inventing three equal energies and avoids quadrature
aliasing of the narrow branch. Ground-state natural width5.57(25)eV is tiny
relative to keV bins; treating q_gs as fixed is an explicit narrow-width
approximation, whose use must be recorded rather than confused with exact
nuclear data. Laboratory rotation/boost correlations remain a separate step.

The ordinary amplitude/grid remains a nonrelativistic R-matrix convention.
Do not pass its nonrelativistic momenta directly to the on-shell Lorentz
boost API. Laboratory event composition must use one explicit consistent
kinematic mapping, retain the declared amplitude approximation, and quantify
the induced source-shape difference while preserving complete event budgets.
The source-state manager cannot repair inconsistent birth kinematics.

The fixed-r16 FSCI policy and1024-node numerical control are now implemented;
see FSCI_SOURCE.md. Incident weights, laboratory composition and
full coupled source validation remain subsequent work. In particular, a
continuum event below the narrow resonance is not automatically an l2 event;
weight entrance components rather than choosing a parent only by energy.

The analytic alpha0 marginal and explicit conditional three-family CM mixture
are now implemented in fusion_pb_birth.h; see PB_BIRTH.md. This supplies no
automatic incident-energy branching or laboratory angular samples. The
remaining composition steps above are unchanged.
