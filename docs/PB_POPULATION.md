# Explicit pB entrance-energy population weighting

The new C/Fortran APIs combine an explicit effective population model with
thermal/beam reaction weighting and selected-energy moments. This is a step
toward physical source composition, not an unfolded channel fit or a complete
laboratory birth source. `fusion_pb_population.h` records the full contract.

## Evidence and deliberate model boundaries

[Laursen2016](https://arxiv.org/abs/1604.01244) reports a5.1(5)percent narrow
Be8 ground-state peak for the16.11MeV parent. Its section5.2 explicitly
excludes the ground-state ghost contribution and discusses an approximately
20percent ghost contribution from its reference17. Do not silently turn
5.1percent into all ground-state-mediated events, or convert the stated ghost
contribution into an independently calibrated narrow delta component.

[Taskaev2024](https://bnct.inp.nsk.su/publics/2024/NIMB_2024_3a.pdf) tables5/8
supply alpha0/alpha1 angle-integrated cross sections at target-depth-averaged
proton energies. Its broad energy averaging and two-angle reconstruction
remain part of the evidence. See TASKAEV_2024_CHANNEL_EVIDENCE.md and the
versioned CSV. The point values are not unfolded resonance-separated data.

We therefore expose three named event populations: ground-state PEAK,
narrow-fit REMAINDER, and OTHER REMAINDER. Neither remainder is asserted to
be a pure alpha1 branch or a measured angular state. In particular the
narrow remainder may contain the ghost. Mapping it wholesale onto the
existing l2/secondary-l2 model would be an additional approximation requiring
source-weighted sensitivity evidence, not something this API does secretly.

## Model equations

Within the existing E<=400keV fit piece, let

`S_n=18200/[(E_keV-148)^2+2.35^2]`,
`r=S_n/(C0+C1 E+C2 E^2+S_n)`.

Use the selected TB/NS/C0±12 coefficients already exposed by the rate model.
For E>400keV set r=0 because this term is absent from that fit piece, not
because a physical pole vanishes there. The diagnostic narrow cross section
is r times the unchanged total model cross section. At E=0, cross section is
zero and r has its limiting S-factor meaning.

The effective continuum peak proxy g is a linearly interpolated table of
`alpha0/(alpha0+alpha1)` from Taskaev's20 rows247–2186keV. Convert relative
CM E to equivalent stationary-proton lab energy with canonical masses:
`Ep_equiv=E*(mp+mB)/mB`. Rows75/134keV are excluded to avoid using the most
strongly smeared low-resonance rows as continuum anchors. The remaining rows
are still target averages and not pure continuum measurements: treating their
ratios as this continuum proxy is an explicit modeling assumption. No
resonance subtraction/deconvolution is claimed. Endpoints are held outside
the table interval; that coverage loss is reported independently.

Caller parameters are isolated narrow-peak fraction b in[0,1] and continuum
scale s in[0,2]; neither is hidden as a universal default. Define z=s*g:

`f_peak = r*b + (1-r)*z`,
`f_narrow_remainder = r*(1-b)`,
`f_other_remainder = (1-r)*(1-z)`.

The tabulated maximum g is below0.1, so the allowed scale does not require
clipping. The three fractions are nonnegative and sum to1. Scale0/2 and
variation of b are model sensitivities, not statistical confidence bounds.
The continuum extrapolation diagnostic is `(1-r)` outside the table and0
inside. It is a subset coverage indicator, not a fourth event population or
an error bar. No total fusion/elastic branching factor is applied a second
time to the existing fusion cross section.

## Weighted reaction rates and energies

For each fraction f, calculate `<sigma v_rel f(E)>` and the same weighted
integrals of projectile, target, relative and CM kinetic energy. The existing
nonrelativistic speed-space GK61 kernel is reused with a default-zero optional
population selector, retaining old API behavior. Gaussian, nuclear and table
knots split the integral; all declared continuation pieces are included.
Probability diagnostics carry the same fraction in the kinematic pair
distribution; they are not normalized reaction fractions. Sum the three
populations to recover total rate, energy moments and probability. The
extrapolation result overlaps the populations and must not be added again.

For unequal thermal temperatures, `Tr=(mb Ta+ma Tb)/(ma+mb)`, and
`C=1.5 Ta Tb/Tr`. At relative energy E, conditional energies are

`Ea=ma/M*C + mb/M*(Ta/Tr)^2 E`,
`Eb=mb/M*C + ma/M*(Tb/Tr)^2 E`,
`Ecm=C + mu/M*((Ta-Tb)/Tr)^2 E`.

The C term is energy, not an extra source, and these quantities satisfy
Ea+Eb=E+Ecm. Thus a branch need not consume the ordinary Maxwellian mean
energy per reacting particle. Beam cases retain relative-direction
conditioning and the exact cold-target limit of the existing classical
model. Proton and boron roles may be exchanged, but the returned canonical
mass pair is required explicitly; the equivalent proton energy is defined
independently of which species is the beam.

## Validation and observed model coverage

GNU/Intel default/Intel-r8 each46/46 tests; C++-only28/28; installed C11/r8
consumers and BALDUR Intel-r8 build pass. Old-model totals remain exactly
identical in the parity tests. Point fractions, all low variants, cold beams,
unequal thermal baths, reversed projectile roles and invalid inputs are tested.

An independent mpmath40-digit ENERGY-space calculation verifies40 component
cases (six unequal-temperature thermal pairs and two finite-temperature beam
cases, each with five outputs). All five reaction/energy moments agree with
maximum relative difference1.6741e-13. It uses an independent cross-section
transcription and energy-space probability density, not the production
speed-space implementation. IEEE-underflow references use a1e-300 denominator
floor in the recorded metric; this is a comparison convention, not a physical
source floor. The initial unscaled reference under-resolved the very small
1000keV beam narrow component (K~2.8e-209); factoring out its Gaussian
exponent recovered high relative precision. That failed reference log is
retained; no production threshold or physics was altered to pass it.

In the explicit demonstration b=.051/s=1/TB and Tb=.7Tp, narrow-remainder
reaction fractions are about14percent at Tp10/50keV,0.66percent at150keV,
and0.17percent at300keV. These are results of this effective model. The
continuum table-extension event fraction is about99percent at2keV,
85percent at10keV and24percent at50keV. Those large fractions forbid treating
low-temperature source calibration as established by the tabulated data.
The actual absolute source and transport impact remains a case-level task.

## Remaining required work

Source-weighted spectrum assumptions, including the narrow-parent ghost and
continuum angular/source alternatives, need explicit choice and sensitivity
accounting. Full birth composition must retain the distribution of reacting
energies and center-of-mass motion; replacing it by a source evaluated at its
mean E is not generally valid. Laboratory boosts/angular correlations,
self-consistent burn/collision/handoff state, host atomic acceptance/I007 and
EXL four-fuel0.6nG full-window validation remain incomplete.
