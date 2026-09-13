# Distribution-integrated thermal product births

`fusion_thermal_birth.h` and `fusion_thermal_birth_fortran` integrate laboratory
product energy sources for all five canonical channels. This increment is
stateless and independent of BALDUR. It does not yet advance coupled thermal
baths, fast-particle collisions, escape, or host transport. The supported
closure is two zero-drift isotropic Maxwellians at the **same temperature**,
with isotropic global rotation of the outgoing CM event. Unequal-temperature,
beam, pitch and current sources are not implemented by this API.

## Source coefficients and accounting

Inputs use joules and SI masses. Output `birth[species*cells+cell]` is a number
coefficient in m³/s, for species 0..5 and neutron 6. Fortran uses `(cells,7)`.
Multiply by `n_a*n_b/(1+delta_ab)` exactly once, including the DD symmetry
factor, to obtain m⁻³ s⁻¹. The returned rate and selected-reactant energy
moments use the same finite quadrature as the births. Use those together;
using a full-reference reaction rate with a truncated product source would
break conservation. Neutron energy is explicit and must be assigned an
escape/deposition owner by the caller. Grid spills are explicit number and
energy coefficients; they are never silently removed or thermalized.

The sum over product energies, including spills, is
`Q*reactivity + reactant_energy_moment[0] + reactant_energy_moment[1]`.
Reaction-selected reactant energy is a debit from the thermal distribution;
it is generally different from removing 3kT/2 per consumed particle. The
canonical nuclear Q is used throughout; this API does not mix the historical
8.68 MeV compatibility constant into new kinetic budgets.

## Independent relative and center-of-mass integration

For common temperature T in energy units, relative energy E and CM velocity
are independent under the incoming Maxwellians. The reacting measure is

`dK = sqrt(8/(pi*mu))/T^(3/2) * E*sigma(E)*exp(-E/T) dE`.

With `C=M*V²/2=T*y²`, the radial CM probability is
`4/sqrt(pi)*y²*exp(-y²) dy`. Both distributions are integrated, not replaced
by their mean energies. Composite Gauss-Legendre integration splits relative
energy at thermal scales, cross-section boundaries, resonance scales and the
actual population interpolation knots. A separate speed quadrature covers
0 <= C/T <= cm_max_kT. Finite tails are not normalized back to one.

The angle-averaged selected reactant energies are
`Ea=(ma/M)*C+(mb/M)*E` and `Eb=(mb/M)*C+(ma/M)*E`.
A representative perpendicular relative/CM orientation is used only to
construct the aggregate parent: its total E and P do not depend on that
angle in the chosen classical budget. This angular reduction is valid for
this scalar-energy, globally isotropic closure; it does not sample the joint
individual reactant-angle distribution or implement anisotropic nuclear
emission. See REACTION_EVENT.md for the classical-budget convention.

## Product spectrum and kinematic mapping

Two-body channels use the existing exact on-shell CM kinematics at the
actual aggregate available energy A. For a product CM kinetic energy K,
rest energy r and boost beta, an isotropic rotation gives a uniform lab
energy interval with midpoint `gamma*K+(gamma-1)*r` and half-width
`gamma*beta*sqrt(K*(K+2*r))`. The lower endpoint is evaluated using the stable
identity `lo=(K-(gamma-1)*r)^2/hi`.

For pB, PB_POPULATION.md supplies the effective ground-state peak weight.
Its exact narrow-width sequential three-alpha events use caller-supplied
intermediate q. The two remainder populations require an **explicit source
proxy**: entrance proxy maps narrow remainder to l2 and other remainder to
the chosen broad model; all-low and all-broad alternatives enable sensitivity
studies. None of these options is a calibrated decomposition of the ghost
and continuum. No new empirical ghost probability is invented.

The existing NR coherent alpha amplitude quadrature is factored into a
private helper without changing its public spectrum outputs. Event weights
are evaluated at `A0=Q+E`. A common squared-momentum scale s then maps each
NR event to on-shell products satisfying the actual parent A:

`K_i(s)=sqrt(r²+2*r*s*e_i)-r`, `sum_i K_i(s)=A`.

The common scale preserves zero CM momentum. This is a declared kinematic
mapping of NR amplitudes, not a new relativistic nuclear-amplitude model.
Diagnostics report the maximum A versus A0 shift and the maximum deviation
from a simple energy-proportional allocation. They are sampled maxima over
the evaluated source support, not probability-weighted error bounds.

Uniform lab intervals are integrated analytically against the linear hats
on arithmetic cell centers. This preserves number and first energy moment,
with separately counted support below/above the extreme centers. Number and
energy closure alone do not establish spectrum accuracy: nq/ncos, energy
quadrature and output-grid convergence are separate requirements.

## Domain and numerical acceptance

The header freezes explicit parameter ranges, layouts and policies. The C
cell count is 1..100000. On error, the result clears; the birth array clears
when its cell count is valid. The Fortran wrapper additionally validates array
extents before C_LOC and clears its actual arrays on every failure.

`reference_*` values use the existing full rate-model integrator, independent
of the finite product-source quadrature. The signed rate discrepancy and
maximum relative reactant-debit discrepancy include truncation and numerical
error; they are not nuclear-data uncertainty estimates. A success status
only certifies finite conservative evaluation, not accurate tail resolution.
The caller must check convergence and discrepancies before advancing state.

`cm_retained_probability` is the actual numerical quadrature sum and may
exceed one by roundoff; it is not clamped. `cm_tail_probability` and
`cm_tail_energy_moment_J` are analytic omitted Maxwellian CM moments:
`tail=erfc(sqrt(X))+2/sqrt(pi)*sqrt(X)*exp(-X)` and
`T*(1.5*tail+2/sqrt(pi)*X^(3/2)*exp(-X))`. Neither is itself a complete
reaction-source energy error bound. Relative-energy truncation and pB model
coverage/proxy uncertainties remain separate.

The present implementation directly integrates events and is intended as a
reference source evaluator. Reusing tables/caches with verified interpolation
will be needed before calling expensive pB spectrum quadrature in every host
zone and trial. It has no hidden global cache or state.

## Validation

`test_fusion_thermal_birth_width.cpp` independently integrates relative energy
and CM speed with composite midpoint sums (512 and 1024 points per axis),
using direct aggregate invariants and analytical angular second moments. It
calls neither the production Gauss quadrature, parent, boost nor grid mapper.
For DD(Tp), DD(He3n), DT and DHe3 at 20 keV it checks full rate parity and each
product's spectral second moment. Positive linear-hat interpolation increases
the second moment by at most `N*h²/4`; the test checks this signed bound plus
independent integration error. The 4000-cell 0..20 MeV grid covers the full
chosen support. This catches missing CM broadening that moment-only ledger
tests would miss. It validates the implemented isotropic closure, not an
experimental differential cross section.

The Fortran test checks C layouts, DT and pB evaluation, conservation, invalid
shape and policy clearing. The separate C++ source regression covers numerical
rate/debit convergence, five-channel accounting and explicit spills. GNU, Intel default-real and Intel `-r8` each pass 51/51 tests;
a C++-only build passes 31/31. Installed C11 and Intel `-r8` consumers pass.
Compiler matrix logs and complete GNU test output are retained under
`docs/validation/thermal-birth`. Full pB source-spectrum refinement and self-consistent coupled
bath/host validation are still required before production activation.
