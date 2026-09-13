# Measured Maxwellian kinetic-to-fluid handoff

`fusion_handoff.h` provides an independent, caller-owned trial operation.
It follows the distribution-proximity approach investigated in
[THERMAL_TRANSFER_DESIGN.md](THERMAL_TRANSFER_DESIGN.md), while computing the
resolved distribution norm directly instead of estimating a particle entropy.
This operation forms a candidate fluid source; it does not advance a thermal
bath or assume that all particles in a thermal-scale kinetic component have
already equilibrated.

## Grid probabilities and eligibility

For target thermal energy T=kT in joules, the normalized energy Maxwellian is

```
m(E) = 2/sqrt(pi) * sqrt(E)/T^(3/2) * exp(-E/T).
```

`fusion_c_maxwellian_energy_grid` integrates this over every energy cell using
regularized Gamma(3/2). Below/above tails are returned explicitly; probabilities
are never renormalized to pretend a truncated grid covers the full Maxwellian.
The upper tail uses the survival function to avoid subtracting CDF values near
one. `represented_mean_energy_J=sum(q_i*Ecenter_i)` is the midpoint-represented
IN-GRID energy per full-distribution particle, not necessarily1.5T. It is a
useful reconstruction diagnostic, not a replacement for physical fluid energy.

For candidate cell-integrated populations Ni in m^-3, N=sumNi and
U=sumNi*Ecenter_i. `fusion_c_maxwellian_handoff_trial` measures

```
L1 = sum_i |Ni/N - qi| + q_below + q_above
mean_error = |U/(1.5*T*N)-1|.
```

Both must meet their explicit caller tolerances before the WHOLE candidate
is projected. This is COMPLETE L1, unlike the half-L1/total-variation curves
in the reference-study figures. The energy criterion is independent: a sparse
high-energy tail can have small probability mass but substantial energy.
The tests deliberately exercise that case. Empty candidates do not project.
These are bin-probability measurements, not bounds on an unknown continuous
shape inside each cell. Grid convergence remains necessary.

## Sources and atomic acceptance

On success with `projected=0`, the kinetic trial equals the old distribution
and fluid/correction sources are zero. Failure to meet a physical criterion
is a valid retained-kinetic result, not a numerical error.

On success with `projected=1`, the kinetic trial is zero and the returned
amounts are

```
Nfluid = Nold
Ufluid = 1.5*T*Nfluid
bath_correction = Uold-Ufluid.
```

The correction is signed energy TO the target bath. It is not additional
fusion Q or a second copy of Uold. The ledger checks actual returned-precision
N/E balances. A host must add Ufluid and the correction exactly once and must
validate that its thermal reservoir can accept a negative correction. The
whole transaction—kinetic removal, fluid source, bath energy and cumulative
counters—must commit or be discarded together. Returned quantities are trial
AMOUNTS, not rates; divide by the host's accepted timestep only at an adapter
which requires a rate. Repeated calls cannot mutate accepted state.

If the target temperature changes, evaluate a fresh trial from the accepted
candidate. The API alone does not prove that a projected Maxwellian remains
valid under later rapidly changing or unequal-temperature backgrounds. The
host fluid closure, bath evolution, spatial transport, charge/pressure and
any re-entry into kinetic treatment remain explicit coupling responsibilities.

## Validation

- Independent erf-based Gamma(3/2) CDF comparisons, exact-grid tails and
  normalization; non-renormalization of truncated grids.
- Resolved Maxwellian projection, rejected hotter/colder candidates, a sparse
  energetic tail with small L1 but excessive mean energy, empty candidates,
  positive/negative bath corrections, rollback and malformed/overflow inputs.
- C/Fortran interoperable4-double grid and13-double handoff ledgers; explicit
  `c_double`/`c_int`, checked extents before C_LOC.
- GNU and Intel default/r8 each31/31 CTest cases; no-Fortran20/20; installed
  C11 consumer and BALDUR Intel-r8 link. The Fortran test initially required
  bitwise zero correction for an ideal1.5T cell; the actual returned-precision
  budget retains a tiny rounding correction. The assertion now uses a
  four-machine-epsilon energy bound. No physical tolerance was loosened.
- Public API composition reruns the established4000-cell/2000-step10keV e/D/T
  alpha relaxation reference with epsilon.001. First projection is0.629s,
  four projection events occur, final fluid fraction is0.999999999979.
  Maximum half-L1 against full kinetics is1.982e-4; N/E residuals are below
  7.4e-15/1.8e-16. At the200 sampled times, public and original prototype
  particle/energy curves agree at stored precision; largest shape diagnostic
  difference is3.01e-17. This verifies API composition, not a new independent
  continuum benchmark. Previous grid/time/tolerance studies are retained.

The physical example has fixed equal-temperature trace baths; it is not an
EXL-50U case or full changing-background acceptance. Raw data, test logs and
retained test failure are in `docs/validation/thermal-handoff/`.

## Use and reproduction

C/C++ callers include `fusion_handoff.h`. Fortran callers use
`fusion_handoff_fortran` (prefer explicit ONLY imports when combining modules)
and link the exported CMake target `pb11::fusion_handoff_fortran`.
`examples/handoff_consumer.c` is an installed C11 usage example.

```
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j4
ctest --test-dir build --output-on-failure
c++ -O2 -std=c++17 -DFUSION_STUDY_PUBLIC_HANDOFF=1 -Iinclude \
    tools/study_fluid_handoff.cpp build/libpb11.a -o build/public_handoff
build/public_handoff 4000 2000 .001 > build/public-handoff.csv
```

Boost is a build dependency; there is no GSL or Python runtime dependency.
