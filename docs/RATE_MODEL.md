# Explicit cross-section continuation and rate model

This interface supplies a declared closure to the previously finite-window
integrals; it does not turn unmeasured cross sections into certified data.
It is independent of BALDUR. All inputs are SI, density-free, and have no
identical-pair factor. Existing finite-window interfaces keep their contract.

## Physics choice

The central selectable model uses the published fits internally and constant
endpoint astrophysical S outside their windows, with the channel's existing
fixed Gamow exponent. This is an explicit extrapolation assumption for BH.
For pB it matches the high-energy thermal-tail assumption of TB2023, whose
analytic low piece already extends to zero. See NUCLEAR_WINDOW_STUDY.md and
its primary-source audits for domains, formulas and evidence limitations.
The API takes a required policy; it never silently chooses a fallback.

A second policy replaces only the high continuation by constant endpoint
cross section. This is a sensitivity alternative, not an upper bound. For pB
there are separate low-piece alternatives: NS coefficients as transcribed by
TB2023, and TB C0+/-12MeV barn. They replace only E<=400keV; none is claimed
to be the complete NS rate model or an experimental confidence interval.
Non-pB channels reject nondefault pB options to expose configuration mistakes.

Model acceptance in a reactor requires both absolute budgets and changes
under declared alternatives. A probability-tail fraction is not an error
bound. A small common-Maxwellian continuation contribution cannot certify
an arbitrary fast-particle distribution. The several-percent low-energy
parameter sensitivity is not removed by precise numerical integration.

## API and accounting

`fusion_c_cross_section_model` evaluates the explicit cross-section model.
`fusion_c_beam_maxwellian_model` integrates a monoenergetic projectile against
a stationary isotropic Maxwellian target; zero target temperature has an
exact cold-target path. `fusion_c_thermal_pair_maxwellian_model` handles two
zero-drift Maxwellians with unequal temperatures and reaction-conditioned
reactant energy moments, preserving their CM/relative correlations.

Each returns fit, below, above and total records using the existing moment
layout. For every segment, the rate and projectile, target, CM and relative
energy moments are integrated without renormalization. Total is their sum;
projectile+target = CM+relative. Only fit owns endpoints and the E=0 limit;
this matters for exact cold targets. A pB below-fit segment has zero measure.

The total domain_incomplete=0 flag means the DECLARED MODEL covers the
positive energy axis, not that the nuclear-data domain has expanded. Total
unresolved probability/speed are zero. Segment unresolved diagnostics refer
to each segment's complement and must not be added: those complements overlap.
`quadrature_error_m3_s` is a numerical rate-error estimate, not model uncertainty
or a rigorous interval. Energy identities are independently checked.

The implementation shares the established relative-speed/Langevin quadrature
with the window API. Gaussian integration extends40 thermal speed units from
the beam drift; omitted tails are below representable scale for these models.
The preexisting numerical restriction v_beam/v_thermal<=1e6 is retained;
zero target temperature is handled analytically. Classical kinematics and
stationary isotropic Maxwellians are explicit assumptions. High-energy
numerical diagnostics do not validate nonrelativistic physics there.

All non-NULL outputs are cleared on error. No shared mutable state, host
indices, I/O, fuel depletion, product generation, ash transfer, or time advance
is hidden in evaluation. Burn/source composition must use consistent nuclear
masses and Q and apply event counting separately.

## Verification

GNU/Intel-default/Intel-r8 each pass38 tests; noFortran passes24. New regressions
cover frozen independent thermal references and model variants, exact cold
segment ownership, unequal-temperature interchange, energy identities and
error clearing. Root's65-case grid independently matches the energy-space
study within1.37e-13 in rates and1.36e-15 in conditional relative energy.
Four direct target-velocity/angle integrations agree within7.25e-12. Installed
C11 and Intel-r8 Fortran consumers pass. See docs/validation/rate-model/ for
raw outputs, test logs, source and reproducible consumer examples.
