# Thermal-transfer model investigation

Status: numerical model study, not a production thermal-ash API. The complete
program goal remains active. Existing first-cell removal is a numerical
kernel option, not a physical definition of thermal helium.

## Resolved reference and first passage

A trace population of 3 MeV alpha particles relaxes in fixed Maxwellian
10 keV electron/D/T baths, ne=1e20 m^-3, nD=nT=5e19 m^-3, lnLambda=15.
The reference retains all alpha particles on the energy grid and evolves the
full previously validated Coulomb drift/diffusion operator. An alternative
removes particles on first arrival at a selected low-energy absorbing cell.
The alternative's removed number and carried energy remain separate ledgers.
If those particles are instantaneously mixed into a Ti Maxwellian, the ion
energy correction is carried_energy - 3Ti*removed_number/2. Omitting this
correction would change the energy budget.

Both calculations conserve total particles and energy to better than1.1e-14
and7.7e-15 in this study. Nevertheless, they diagnose different things:
first-passage absorption versus overlap of the full distribution with a
Maxwellian. Neither can simply be renamed the other's thermalized fraction.
At the refined setting, first-passage t90 is0.5540/0.5436/0.5312s for lower
edges2.5/10/40keV, whereas the full-reference Maxwellian-overlap t90 is0.5336s.
These differences are not all numerical error. It would be unjustified to
pick the40keV boundary merely because that one diagnostic happens to agree.

The original logarithmic grid smeared the transient more than the terminal
energy partition suggested. A hybrid grid resolves0.001..20keV logarithmically
and20..6000keV uniformly. For the reference,2000→4000cells at2000steps changes
t90 by0.0002933s;2000→4000steps at4000cells changes it by0.001649s. These are
separate resolution checks. Deposited electron fractions remain near0.80635.
The absorbing penalty scan1e4..1e8 per step reached the numerical boundary
limit; it does not establish a physical thermalization rate.

Reproducibility: tools/study_alpha_transfer.cpp and analyze_alpha_transfer.py;
raw and summary files under docs/validation/alpha-transfer-study/.
The source is a diagnostic executable, not the host production integrator.

## Literature-motivated two-component prototype

Peigney et al., arXiv1402.6191, sections3.2–3.3 and4.6, define suprathermal
and thermal-scale alpha distributions that both exist over velocity space.
Their sum is physical. The thermal-scale component continues kinetic
relaxation; it is not immediately a Maxwellian or the host's thermal ash.

For an isotropic fixed-Maxwellian-bath reduction, the cold-ion suprathermal
energy drift can be written using

```
nu_i = 4*pi*n_i*Za^2*<Zi^2>*(e^2/(4*pi*epsilon0))^2*lnLambda/(ma*mi)
dv/dt = -sum_i(nu_i)/v^2.
```

This is the high-v/Ti limit of the existing normalized Coulomb energy drift:
A_E=-C_i/(mi*v), with C_i=4*pi*n_i*Za^2*<Zi^2>*ke^2*lnLambda.
nu has units velocity^3/time. Finite-volume flux at v=0 gives a transfer
rate proportional to4*pi*sum(nu_i)*fST(0), rather than a chosen energy cutoff.
A positive implicit downward flux on velocity volumes
(v_hi^3-v_lo^3)/3 implements this part. Electron drift/diffusion remains in
ST; the thermal-scale component receives the transfer and evolves under the
full finite-temperature electron/ion Coulomb operator.

The transfer from ion i is distributed according to that ion's velocity
Maxwellian. Its alpha energy temperature is therefore (ma/mi)*Ti, with
continuous mean energy(3/2)*(ma/mi)*Ti. In a discrete source, its number must
normalize correctly and its represented energy must be measured. Assigning
that finite injection energy to the thermal component requires a matching
ion energy debit when using a cold-ion ST approximation. A number-conserving
transfer alone does not establish energy conservation.

Root's prototype uses a shared grid first, explicit per-ion transfer weights,
full TH kinetics, and a separate fully resolved reference. It must be checked
for number/energy closure, full-distribution differences, source quadrature
error and separate time/grid refinement before any public API is accepted.
This reduced, isotropic, fixed-bath prototype is not a reproduction of the
paper's multidimensional algorithm. Fluid thermal-He handoff, unequal bath
temperatures, changing backgrounds, escape and host transport remain further
acceptance requirements.

## Cold-ion prototype outcome

The cold-ion prescription above was retained as a research comparison, not
promoted to the production transfer model. Its final distribution converges
to the same equilibrium but its transient does not reproduce the resolved
finite-temperature operator. Sampled maximum normalized distribution half-L1
errors for cells/steps1000/1000,2000/2000,4000/2000,4000/4000 are respectively
0.05216,0.06541,0.07383,0.07447. At4000/4000 the final electron heat differs
by-0.001449 of the initial alpha energy. Final half-L1~6.6e-14 would conceal this
transient discrepancy. Spatial/time refinement and closed energy budgets
are necessary but do not eliminate model error. This is an isotropic
cold-ion prototype, not a verdict on Peigney's multidimensional algorithm.

## Finite-temperature shared-grid decomposition

We therefore retain the full finite-temperature Coulomb operator L for BOTH
kinetic components and subtract/add a localized internal transfer:

```
dS/dt = L S - Lambda S
dT/dt = L T + Lambda S.
```

S+T obeys the original full linear FP equation exactly. This shared-grid
version deliberately does not use Peigney's cold-ion delta-function
replacement or a coarse/fine grid interpolation. For each explicitly selected
ion bath, use the divergence of its finite-temperature radial friction:

```
vti = sqrt(2*Ti/mi), y = v/vti
H(y) = erf(y) - 2*y*exp(-y*y)/sqrt(pi)
F_i(v) = nu_i*H(y)/v^2
lambda_i(v) = (1/v^2)*d[v^2*F_i(v)]/dv
            = 4*nu_i/(sqrt(pi)*vti^3)*exp(-y*y)
Lambda = sum_ions(lambda_i).
```

Here F is the inward friction magnitude, not the full stochastic mean
velocity/energy drift. Starting with the existing energy diffusion D_E,
F_i=D_E/(ma*v*Ti). For p(E)=4*pi*v*f(v)/ma, the source lambda*f in velocity
space becomes exactly lambda*p in energy space. The transfer therefore moves
particles at the SAME cell energy; no bath energy correction accompanies
this internal relabeling. Lambda(0) is finite and nonnegative. Its velocity
integral is4*pi*nu_i, recovering the localized finite-temperature source
normalization. Electrons remain in L but are not selected as ion-transfer
baths. nu_i uses the independently validated dimensional NRL normalization;
we do not transcribe a potentially inconsistent dimensionless paper prefactor.

The coefficient is independently tested by reconstructing v^2*F from the
existing Coulomb diffusion, applying a fourth-order numerical derivative,
and integrating lambda over velocity space. For D/T baths and0.001<=v/vti<=3,
maximum derivative relative disagreement is1.21e-9. Simpson integration
against the independently measured high-speed friction flux gives the
4*pi*nu normalization within1e-12. Tests also cover zero charge, the finite
zero-energy limit, tails, invalid arguments and numerical overflow.

For a frozen background and common physical escape E, use:

```
(I - dt*L + dt*E + dt*Lambda) S_new = S_old + dt*birth_S
(I - dt*L + dt*E) T_new = T_old + dt*birth_T + dt*Lambda*S_new.
```

Adding these equations recovers the SAME backward-Euler update for S+T.
The positive FP matrices preserve nonnegative states without clipping.
Per-bath collision heat is linear, so summed S/T heat must also reproduce
the single-component reference. Internal transfer number/energy is recorded
separately from physical escape, external births and bath heat. Trial state
is caller-owned; no rejected call changes an accepted population or counter.

This decomposition is a mathematically exact representation of the stated
isotropic test-particle equation, NOT a unique observable definition of
thermalized particles. T remains kinetic and must not be injected directly
into BALDUR's thermal-He transport equation. Variable backgrounds must use
the same frozen coefficients in both solves at each stage. Spatial transport,
fast pressure/charge, bath evolution and the host integrator remain separate
coupling responsibilities.

## Next physical handoff acceptance

Ricketson et al., arXiv1301.5678 / JCP273(2014)77–99, use distribution
proximity to a Maxwellian to control a hybrid kinetic/fluid approximation
(Eqs22–31), rather than the instantaneous speed of a particle. Their
relative-entropy bound motivates an independently measured distribution
error criterion here; their Monte Carlo algorithm and its tuned numerical
constant are not being copied into this deterministic energy-grid solver.

The next fluid handoff must measure both distribution and energy errors.
The thermal-scale label alone is insufficient. A prospective grid-resolved
criterion compares a normalized candidate distribution to Maxwellian bin
probabilities, includes finite-grid tails, and records the removed energy.
If the host puts removed N into a common-Ti fluid with energy3*Ti*N/2,
the separate bath correction is U_removed-3*Ti*N/2. It must not be silently
lost or counted twice. Any projection tolerance must be explicit and refined
against the full kinetic reference. This handoff, evolving unequal-temperature
backgrounds and host transport are not certified by the exact split identity.

Primary references:
- https://arxiv.org/abs/1402.6191 (Peigney et al.,2014).
- https://arxiv.org/abs/1301.5678 (Ricketson et al.,2014).

## Exact-split study result

The1000-cell/1000-step independent composition has maximum over ALL steps
normalized half-L1(S+T,FULL)=8.47e-17, maximum normalized accumulated
electron/D/T heat differences3.35e-15/2.60e-16/4.05e-16, and maximum combined
number/energy residuals1.23e-15/1.56e-15. This establishes the algebraic
composition in a physical slowing case; it does not assert that1000cells
resolve all continuous transient physics. Separate refinement is documented
in the preceding resolved-reference study.

Diagnostic half-L1 is (1/2)*sum(abs(N_i-Nref_i))/N0; for normalized
distributions it is total variation. The raw prototype originally called
this L1 despite the explicit half factor. Source/CSV/log field names were
corrected on import without changing values or rerunning calculations.
The thermal-scale Maxwellian diagnostic uses the same half-L1 convention.
An initially negligible TH population may still have order-unity normalized
shape error; density and shape must be interpreted together.

## Controlled fluid-handoff research prototype

A subsequent root-owned fixed-bath prototype tests actual Maxwellian
projection, separately from the public exact-transfer API. For candidate TH
cell populations N_i, let p_i=N_i/NTH. Maxwellian q_i is the analytic integral
of the3D energy Maxwellian over each cell (regularized gamma3/2), not a
midpoint sample. Exact tails outside the grid are included in the complete
bin-probability L1 test. Project all TH only when BOTH

```
sum_i |p_i-q_i| + q_outside <= epsilon
|U_TH/(1.5*Ti*N_TH)-1| <= epsilon.
```

The second condition is essential: a small particle-number norm alone does
not bound the energy of a rare high-energy tail. Projection records NTH and
its ACTUAL midpoint-represented energy. Fluid energy is1.5*Ti*NTH; the
difference is assigned to the common-Ti ion bath (D/T split by equal number
in this particular fixed-bath example). There is no selected particle-energy
cutoff and no negative-state repair. No projection occurs if either criterion
fails. The criterion is grid-resolved, not a bound on unknown within-cell
continuous structure; finite-grid representation still needs convergence.

The five one-second runs below preserve N/E to8e-15/2e-16. Maximum half-L1
compares S+TH+reconstructed fluid Maxwellian to a separately advanced full
kinetic reference throughout the window:

| cells/steps | epsilon | first projection(s) | max half-L1 |
|---|---:|---:|---:|
|2000/2000|.001|none|2.38e-16|
|4000/2000|.001|.6290|1.982e-4|
|4000/4000|.001|.6255|1.982e-4|
|8000/2000|.001|.6260|2.702e-4|
|8000/2000|.00025|.6575|4.954e-5|

The2000-cell grid never meets the shape criterion; this is a retained
insufficient-resolution result, not a thermal-ash success. The fixed-tolerance
maximum error need not decrease monotonically under refinement because the
projection event changes; lowering epsilon at the SAME8000/2000 resolution
reduces the discrepancy. No threshold was loosened to make the coarse case
project. At4000/2000, energy ledger closure is~2e-16 even though reconstructing
the fluid by midpoint energies has a separate grid-energy bias2.40e-6 of
initial energy. At8000cells this bias is~6.0e-7. These quantities must not be
conflated. The sum of the local projection shape errors alone also omits the
fact that exact bin-integrated Maxwellians are not precisely stationary under
the sampled-Maxwellian FP discretization. The full reference comparison
includes that discretization effect.

This prototype is NOT yet a public fluid-handoff API or a BALDUR thermal-He
source. Its equal, fixed temperatures make the projected Maxwellian
stationary in the continuous collision problem. Different/evolving bath
conditions, ongoing births, host bath-energy acceptance and spatial transport
must be checked before coupled production use. Next implementation should
factor the measured Maxwellian projection into a reusable old/trial library
operation, keeping its tolerance, actual energy, grid tails and eligibility
explicit. It must not turn failure to meet the criterion into forced ash.

## Public API end-to-end regression

The same study also compiles with FUSION_STUDY_PUBLIC_API=1 to use the public
wrapper, with an independent full FP reference and explicit per-step
half-L1/per-bath heat tolerance1e-11. CTest runs this1000/1000 physical case.
The first wrapper run rejected step818 because a negligible per-cell product
lambda*S_new underflowed double precision. The source now permits this
positive cell-source rounding, like the existing FP positive state tails,
while recomputing and checking TOTAL N/U residuals. No negative clipping,
population floor or relaxed balance threshold was introduced. A small
analytic regression reproduces that specific underflow, and the full public
run now completes one second. The failed log is retained. Public mode's
internal-transfer subledger aliases are bookkeeping only; its independent
checks are total N/U, FULL distribution and per-bath heat.

## Build/ABI validation at this checkpoint

GNU and Intel default/r8 configurations each pass29/29 CTest cases, including
the public1000/1000 slowing test. No-Fortran mode passes19/19. Installed C11
consumer passes (examples/two_component_consumer.c); the explicit C ledger
layout is14doubles. Both Intel configurations use runtime checks. A new
Fortran test initially selected the C++ linker and failed a PIE relocation;
setting LINKER_LANGUAGE Fortran, consistently with existing Fortran tests,
resolved it. The failed build log is retained. No production compiler flag
or physics setting was changed to hide the error.

C++ callers include fusion_two_component.h; Fortran callers use
fusion_two_component_fortran and link pb11::fusion_two_component_fortran.
The independent trace study defaults to separate FP calls; CTest enables
FUSION_STUDY_PUBLIC_API=1 to exercise the wrapper. Source/CSV/logs are in this
repo, with numerical/reference distinction and the prototype limitations
above. BALDUR rebuild/link passes, but these operators are not yet activated
in its source assembly.
