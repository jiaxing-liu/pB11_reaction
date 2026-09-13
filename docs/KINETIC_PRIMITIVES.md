# Isotropic test-particle collision primitives (M3 increment)

These independent C/C++ primitives implement Coulomb coefficients and a trial
energy-distribution step. They do not yet implement a complete fusion-product
spectrum, reaction/kinetic network, thermal-ash closure, host integration or
EXL-50U prediction. Existing instant pB and thermal-network APIs are unchanged.

## Physics and conventions

`fusion_c_coulomb_energy` describes a dilute nonrelativistic test particle in
one fixed Maxwellian bath. All units are SI; temperatures are kT in joules.
The caller supplies positive lnLambda and the bath charge-state second moment
<Z_b^2>. Screening, collective effects, relativistic corrections and fast-fast
collisions are outside this primitive. The bath is not advanced internally.

The coefficient normalization follows the exact test-particle relaxation rates
on printed p31 of the [official NRL Plasma Formulary (2019)](https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Formulary_2019.pdf?ver=p9F4Uq9wAtB0MPBwKYL9lw%3D%3D), which cites Trubnikov (1965):

```
v = sqrt(2 E / ma), x = mb v^2 / (2 Tb)
psi(x) = erf(sqrt(x)) - 2 sqrt(x) exp(-x) / sqrt(pi)
C = 4 pi nb Za^2 <Zb^2> (e^2 / (4 pi epsilon0))^2 lnLambda
nu0 = C / (ma^2 v^3)
nu_parallel = psi(x) / x * nu0
nu_energy = 2 [(ma/mb) psi(x) - psi'(x)] nu0
D_E = ma^2 v^4 nu_parallel / 2
A_E = -E nu_energy
```

D_E is HALF the energy variance per unit time. The factor 1/2 is essential.
The implementation uses an algebraically equivalent, cancellation-safe formula
and a small-speed series. E=0 has D_E=0 and finite positive A_E. A separate test
computes the NRL frequency route rather than copying the implementation form.
e and epsilon0 use [CODATA 2022](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

The energy-density operator for p(E)=dN/dE is

```
J_b = -D_E,b [dp/dE + (1/Tb - 1/(2E)) p]
dp/dt = -d/dE sum_b J_b
A_E = dD_E/dE - D_E (1/Tb - 1/(2E))
```

This retains energy diffusion and has p_eq proportional sqrt(E) exp(-E/Tb).
A deterministic friction-only approximation is not used near equilibrium.

## Finite-volume trial and accounting

`fusion_c_energy_fp_trial` accepts arbitrary increasing nonnegative energy
edges, integrated cell populations N_i (m^-3), and bath-major interior-face
D_E coefficients (J^2/s). Cell energies are arithmetic midpoints; bath kT values
are independent. The flux uses exponential fitting of the sampled Maxwellian,
and implicit Euler gives a tridiagonal M-matrix. A nonpositive pivot or negative
population is an error, never repaired by clipping. Conversion of a positive
tail below double's representable range follows IEEE rounding; the full
particle/energy residual tests still bound any associated loss.

Birth is a cell population rate. Escape removes lambda_i*N_i; the optional
first-cell thermalization rate removes lambda_th*N_0. That rate is an explicit
caller policy, not a validated physical ash prescription. Reflecting boundaries
apply otherwise. Escape and thermalization carry their cell midpoint energy.
A continuous birth spectrum must be mapped preserving number AND energy by a
separate source routine; this primitive does not perform that mapping.

For each bath, returned heat is minus the discrete fast-energy change caused
by that bath's face flux. Positive heat deposits in the bath; negative heat
heats the fast distribution. With the returned state, both budgets must close:

```
Nfinal - Ninitial - Nborn + Nescaped + Nthermalized = residual_N
Ufinal - Uinitial - Uborn + Uescaped + Uthermalized + sum(heat_b) = residual_U
```

Each residual is checked against 1e-10 of its accumulated nonnegative budget
scale. Outputs clear on errors. No hidden state, cumulative counters or I/O:
reject a trial by discarding it, accept by copying its returned state. Host
rollback, multistage extrapolation and restart consistency remain host-level
requirements and have not been certified by these primitive tests.

## Validation and physical limits

- One-cell implicit birth/loss solution, first-order time refinement, reflecting
  particle conservation, common-bath equilibrium, opposite-sign two-bath heat,
  repeated trials and errors.
- Independent NRL coefficient normalization, derivative identity, small-speed
  limits and actual short-step energy-moment convergence.
- `test_fusion_slowing`: trace 3 MeV alpha population 1e12 m^-3 in fixed e/D/T
  baths at 10 keV, ne=1e20 m^-3, nD=nT=5e19 m^-3, supplied lnLambda=15. No
  reaction, escape or ash conversion. Log energy grid 1 eV–10 MeV, 1 s duration,
  200/400/800 cells and 250/500/1000 steps. Initial energy is distributed between
  adjacent cell centers preserving number and energy. All nine cases approach
  the Maxwellian mean 15 keV; electron/ion heat partition converges with grid.
  This is a fixed-bath mathematical/physics benchmark, not a burning plasma.

The benchmark deliberately retains thermalized test particles on the kinetic
grid to verify detailed balance. A production model must instead define and
validate thermal-pool transfer without double counting charge or energy.
The nonrelativistic assumption, fixed lnLambda and trace approximation must be
checked for each application; passing these tests cannot establish validity
for hotter backgrounds, relativistic products or significant fast pressure.

## Build and ABI verification

GNU Fortran and Intel ifx 2026.1.1 with default real and `-r8`, both Intel
variants using `-check all -traceback`, each pass 13/13 CTest cases. The
no-Fortran build passes 10/10. Logs are in `docs/validation/m3-*.txt`.
The installed no-Fortran static library was called successfully by the
C11 program `docs/validation/m3-c-caller.c` (compile with cc, link with c++).
`pb11::fusion_kinetics_fortran` exports explicit c_double/c_int types and
handles empty arrays and noncontiguous Fortran sections. No BALDUR stateful
kinetic coupling is claimed by these ABI checks.
