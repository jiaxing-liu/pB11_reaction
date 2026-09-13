# Nuclear Coulomb functions and coherent alpha amplitudes

Work in progress: this layer implements an explicitly stated amplitude and
numerical evaluation, not a validated full-energy pB birth probability model.
The initial 16.11/16.62MeV states, alpha0 branching, l mixing and final-state
corrections require separate model selection and experimental validation.

## Nuclear Coulomb evaluator

The channel IDs and fixed masses/radii/constants are in
fusion_nuclear_coulomb.h. All inputs are SI J. The numerical domain is
0.001..12MeV relative energy, distinct from a measured nuclear-data domain.
Outputs are log(P), shift S, exp(i(omega-phi)), and rho=k*a. These Coulomb
functions describe nuclear barrier penetration, not plasma collision slowing.

Tables are generated with pinned mpmath1.3.0 at40digits. In each log-energy
interval a17coefficient Chebyshev expansion represents logP,S and real/imaginary
phase. Independent off-node checks trigger adaptive subdivision. The generator
checks F'G-FG'=1 and stores constants, segment errors and reference evaluations.
The C++ evaluator uses Clenshaw recurrence, preserves the interpolated phase
error, and fails outside the domain. It performs no fitted nuclear extrapolation
or production-time Python calls. Numerical functions alone do not establish
any cross-section or spectral model.

## Amplitude conventions

fusion_c_alpha_amplitudes takes l_primary in1,2,3, total CM kinetic energy A,
intermediate pair energy q, and a secondary recoil-frame cosine. It constructs
nonrelativistic Jacobi momenta internally, as used by the stated R-matrix
amplitude. The earlier relativistic kinematic generator is a separate API;
its momenta are not silently substituted in this nonrelativistic model.
The parent and intermediate spins are2, and the secondary orbital l is2.

For each cyclic primary choice, the spin-projection amplitude is the sum over
intermediate projections of a Clebsch-Gordan coefficient times primary and
secondary spherical harmonics, i^(l_primary+2), Coulomb/hard-sphere phases,
and radial factors. The intermediate denominator is
E0-q-gamma_squared*(S(q)-S(E0))-i*gamma_squared*P(q),
with E0=3.129MeV, gamma_squared=1.075MeV and a_secondary=4.5fm.
The primary radius is5.1fm; its reduced width squared is fixed at1MeV as an
arbitrary overall amplitude convention. Condon-Shortley spherical harmonics
and the standard CG phase convention are implemented explicitly.

The API returns both a single unsymmetrized amplitude and the coherent sum
of all three cyclic permutations, for m=-2..2. Squaring each term before
summing permutations would discard physical interference. To obtain a density
in dq*dcos multiply the summed spin norm by sqrt(q*(A-q)); the returned
phase-space factor uses J. Overall normalization must be computed by the
caller/model; these are not normalized probabilities. Averaging over initial
spin gives a constant1/5 for this rotationally averaged quantity and does not
supply beam-aligned laboratory angular distributions.

No permutation is silently omitted if a Coulomb evaluation is out of range.
Thus the present finite numerical window does not yet cover every endpoint
of three-body phase space. A full spectrum integrator must resolve this
endpoint policy and demonstrate normalization/convergence before production.
Nor does this API choose the broad-resonance l1/l3 mixing coefficients or the
low-resonance ground-state branch.

## Source-algebra discrepancy under audit

The explicit J=2,l=2,Jb=2 single-permutation spin sum produces an angular
ratio W(theta)/W(0)=1-(9/16)sin^2(2theta). Laursen2016 Eq2 prints a different
standalone angle expression. The library follows the explicit coupled
amplitude; independent exact CG/source review confirms this discrepancy. Further experimental comparisons are required before labeling
its experimental spectrum validated. Do not silently change a CG coefficient
to reproduce that standalone expression. Root tests check cyclic relabeling,
rotational spin-norm invariance, identical-secondary exchange and the
independently derived unsymmetrized polynomial.

Primary formula sources: Laursen et al., arXiv1604.01244, p3 Eqs3-4;
Kuhlwein et al., arXiv2109.07886, p2 Eqs1-3. Parameter/boundary reference
and width reconstruction are in validation/rmatrix_boundary_reference.py.
These references and the open angular discrepancy are documented in the
companion program's RMATRIX_PARAMETER_SOURCES.md and this repository's ANGULAR_CORRELATION_AUDIT.md.
