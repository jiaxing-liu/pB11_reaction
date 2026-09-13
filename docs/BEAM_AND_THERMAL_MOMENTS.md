# Beam/thermal rate windows and reaction-conditioned energy

M3 source-physics increment. All units are SI, with kT in joules. These APIs
are independent of BALDUR and complement the thermal network and kinetic
trial primitives. No population is advanced and no density factors are applied.
The caller supplies consistent reactant masses in kg; the kernel does not
infer or silently substitute masses from a channel ID.

## Data domain is part of the result

`fusion_c_cross_section_domain` returns the existing channel's fitted relative
energy interval. `fusion_c_beam_maxwellian_window` and
`fusion_c_thermal_pair_maxwellian_window` integrate only that interval. They
never extend a nuclear fit, renormalize the retained probability or interpret
missing data as a physically zero cross section. Existing cross sections and
legacy thermal-fit APIs retain their semantics and values.

The windows come from the already documented Tentori–Belloni pB model and
[Bosch–Hale 1992](https://doi.org/10.1088/0029-5515/32/4/I07): see
BOSCH_HALE_DATA.md and THERMAL_NETWORK.md. pB covers 0–9.76 MeV; the positive
DD/DT/DHe3 data windows have nonzero lower endpoints. The isolated exact
zero-speed limit sigma*v=0 is retained for every channel.

`PB11_STATUS_OK` means a window integral was computed successfully. It does
not certify the full physical reactivity. The output contains:

- Resolved reactivity, m^3/s, and quadrature error estimate for that integral.
- Reactivity-weighted energies of reactant a, reactant b, relative motion and
  center-of-mass motion, J m^3/s. Divide by the nonzero resolved reactivity to
  obtain conditional mean energies of reacting pairs.
- Resolved/unresolved relative-velocity probability, plus the unresolved
  integral of relative speed, m/s.
- `domain_incomplete`: always 1 for nonzero relative Maxwellian temperature,
  even when a tiny unresolved numerical integral rounds to zero. For an exact
  cold monoenergetic pair, the flag follows whether its relative energy lies
  inside the fitted interval (or equals the isolated zero-speed limit).

A small unresolved probability does not by itself bound a reaction-rate
error. If independent evidence supplies sigma(E)<=S outside the data window,
then the missing rate is bounded mathematically by S*integral_outside(w P dw).
The reported integral is a numerical estimate, not a certified interval, and
no physical S is invented by this library. These diagnostics do not certify
missing reaction-energy moments either. A production model must explicitly
resolve the domain/uncertainty policy; this increment does not silently
approve a truncated rate as a complete source.

## Monoenergetic projectile / Maxwellian target

The target has zero bulk velocity and an isotropic Maxwellian distribution.
Let v=sqrt(2 Ea/ma), u=sqrt(2 Tb/mb), mu=ma*mb/(ma+mb). Integration of target
velocity directions gives the relative-speed density

```
P(w) = w/(sqrt(pi)*u*v) *
       [exp(-(w-v)^2/u^2) - exp(-(w+v)^2/u^2)]
K_window = integral_window sigma(mu*w^2/2) * w * P(w) dw.
```

At v=0, P(w)=4*w^2*exp(-w^2/u^2)/(sqrt(pi)*u^3). At Tb=0, the code evaluates
sigma at the relative energy mu/ma*Ea and multiplies by v exactly. A beam's
laboratory energy is not passed directly as the nuclear relative energy.
Direction anisotropy of the monoenergetic projectile does not affect this
scalar rate against an isotropic target; product directions and spatial
transport require separate models.

For fixed w, define a=2*v*w/u^2 and L(a)=coth(a)-1/a. The conditional target
and center-of-mass energies follow from the angular distribution exp(a*cos):

```
Eb(w) = mb/2 * [(v-w)^2 + 2*v*w*(1-L(a))]
Ecm(w) = (ma+mb)/2 * [(v-mb/(ma+mb)*w)^2
                     + 2*v*mb/(ma+mb)*w*(1-L(a))]
Er(w) = mu*w^2/2.
```

Each is integrated with sigma*w*P. The nonnegative forms avoid subtracting
nearly equal squared velocities in the cold-target limit. Small-a series and
large-a formulas avoid cancellation in L. The returned energies satisfy
Ia+Ib=Icm+Ir, where Ia=Ea*K for a monoenergetic projectile. Reacting target
particles generally do not carry the unweighted Maxwellian mean 3Tb/2.

## Different-temperature thermal reactants

For independent zero-drift Maxwellians at Ta and Tb, the relative distribution
has Trel=(mb*Ta+ma*Tb)/(ma+mb). The implementation uses a zero-speed projectile
with an effective target temperature Tb+(mb/ma)*Ta solely as an integration
substitution to obtain K and Ir with the existing beam kernel. It does not
change the physical input temperatures.

Gaussian conditioning gives the following conditional energies at fixed Er:

```
M = ma+mb, mu = ma*mb/M
Ea(Er) = (ma/M)*(3/2*Ta*Tb/Trel) + (mb/M)*(Ta/Trel)^2*Er
Eb(Er) = (mb/M)*(3/2*Ta*Tb/Trel) + (ma/M)*(Tb/Trel)^2*Er
Ecm(Er) = 3/2*Ta*Tb/Trel + (mu/M)*((Ta-Tb)/Trel)^2*Er.
```

Derivation: for V=(ma*va+mb*vb)/M and w=va-vb, per-component covariance
Cov(V,w)=(Ta-Tb)/M and Var(w)=Trel/mu. Thus E[V|w] equals
mu*(Ta-Tb)/(M*Trel)*w and Var(V|w)=Ta*Tb/(M*Trel). Substitution into va and
vb yields the expressions above. Equal-temperature CM and relative motion
are independent, so Icm=(3/2*T)*K; unequal-temperature motion is correlated.
Ignoring that correlation gives wrong fuel-energy debits despite a correct
scalar rate. This Gaussian reduction is a root derivation, verified against
independent averaging of the monoenergetic kernel and direct velocity angles.

## Numerics and verification

Integration uses the dimensionless relative speed x=w/u. Gaussian knots,
existing fit joins and resonance scales split the quadrature; an expm1 form
avoids cancellation at small projectile speed. Gauss–Kronrod61 uses relative
request 1e-11 and depth15. A resolved quadrature estimate above 1e-8 relative,
energy identity mismatch above 1e-9 relative, or probability normalization
error above 1e-9 rejects the call. These are numerical consistency gates,
not nuclear-model error bars.

The Gaussian computational span is max(0,s-40)..s+40, s=v/u. Omitted far tails
are below representable precision for the nonrelativistic physical regimes
intended here; numerical underflow is not proof of complete nuclear coverage.
At s>1e6 the finite-temperature kernel rejects loss of numerical resolution;
the exact Tb=0 API remains available. No automatic cold approximation is made.
The models assume classical nonrelativistic motion and uncorrelated isotropic
backgrounds; they do not validate a relativistic or drifting plasma.

`tests/test_fusion_beam.cpp` verifies:

- All-channel cold target kinematics and explicit domain gaps.
- Analytic Maxwellian mean relative speed for small, intermediate and large
  drift, using completely unresolved test windows.
- Continuous approach to the cold target: DT at 50 keV laboratory energy gives
  relative finite-T corrections 0.0140586, 0.00140831, 0.000140855, 0.0000140858
  at target kT=0.1,0.01,0.001,0.0001 keV, without relaxing the test tolerance.
- Independent target-speed/angle integration for pB and DT, including target
  and CM energy moments.
- Maxwellian averaging of the projectile reproduces a separate relative-
  energy thermal integral.
- Unequal-temperature rate and species energy moments against explicit
  projectile averaging; species-exchange and one-cold-reactant limits.
- C/Fortran status, field layout, error clearing and energy identity.

No DD identical-pair factor is inserted here. The event network must apply
1/2 once when integrating two identical populations. The current outputs are
building blocks for conservative reactant depletion and product births;
validated product spectra, grid birth mapping, reaction loss from fast bins,
thermal ash, host source ownership and EXL-50U application remain separate
acceptance work.

## Build evidence

GNU Fortran, Intel ifx default real and Intel `-r8` each pass 15/15 tests;
Intel uses `-check all -traceback`. The no-Fortran build passes 11/11. An
installed static-library C11 caller also passes. Commands follow the prior
kinetic build setup; archived outputs are `docs/validation/m3-beam-*`.
The C example uses explicitly illustrative 2u/3u masses and reports the
resolved window and its incompleteness, rather than claiming a complete
DT source from an unresolved data interval.
