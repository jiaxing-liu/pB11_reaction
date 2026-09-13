# Normalized three-alpha CM source grid

`fusion_c_alpha_spectrum_grid` integrates the declared coherent amplitude
kernel and returns **birth particles per reaction**, not a reaction rate or
an experimentally validated pB spectrum across incident energies. Its C and
Fortran interfaces are independent of BALDUR.

## Model and units

Mode 1, 2 or 3 selects the primary orbital angular momentum. Mode 13 mixes
l=1 and l=3 coherently after separately normalizing the two full phase-space
bases. If `a` and `b` are the cyclic spin amplitudes and
`N1=integral |a|^2 dmu`, `N3=integral |b|^2 dmu`, the mixed density is

```
| sqrt(k)*a/sqrt(N1) + exp(i*delta)*sqrt(1-k)*b/sqrt(N3) |^2 dmu.
```

It is subsequently normalized by its own integral, including interference.
`k` is therefore an input basis weight, not a final branching ratio. At k=1
and k=0 the shapes recover the pure l=1 and l=3 sources. The API reports raw
normalizations for numerical diagnostics; the mixed raw normalization uses
N1 as its overall amplitude reference. These quantities have no cross-section
interpretation. Spin averaging and overall orientation constants cancel in
the normalized scalar energy spectrum.

A is the total available CM kinetic energy in joules. The numerical API
accepts 0<A<=12 MeV, with a finite positive normalization required for success.
This numerical interval is not a validated nuclear model interval.
The measure is `dmu=sqrt(q*(A-q))*dq*dcos(theta)` in J^2. The nonrelativistic
Jacobi momenta and R-matrix conventions are those in ALPHA_AMPLITUDES.md.
Gauss-Legendre integration uses q=A*sin(u)^2, u in [0,pi/2], and cos(theta)
in [-1,1]. Node convergence and weight closure are checked. Each event
produces three energies whose sum is A; event probabilities sum to one.

The output energy grid has arithmetic cell centers. Adjacent-center linear
projection conserves particle number and energy. Particles outside the center
hull remain in explicit below/above number **and energy** ledgers, even if
inside the outer cell edges. They are not silently assigned to thermal ash,
escape or a boundary cell. Including this spill, number=3 and energy=A per
reaction. A caller may multiply by an explicitly selected event rate later.
The API performs no laboratory boost, deposition or state advance.

## Explicit cutoff

The old strict amplitude API continues to reject out-of-window Coulomb
arguments. `fusion_c_alpha_amplitudes_cutoff` accepts a caller-selected
1..10 keV cutoff and phase-space endpoints. It skips only a cyclic amplitude
whose primary or secondary pair energy lies below that cutoff, preserving
other permutations in the same event. It reports the number skipped. This
is a numerical approximation, not a measured nuclear cutoff. Other domain
errors still fail the call and clear the outputs.

`tools/check_alpha_cutoff.py` uses mpmath 1.3.0 at 50 digits and fixed-eta
Coulomb derivatives. It samples A=8.68,8.829,9.3,12 MeV, l=1,2,3, and endpoint
energies equal to cutoff, cutoff/2 and cutoff/10. The largest sampled radial
magnitude at a 1 keV cutoff is 8.71e-39 (secondary endpoint); at 10 keV it is
9.76e-13. All sampled magnitudes decrease toward the endpoint. These are
finite sample maxima, **not rigorous bounds** on omitted amplitudes or
integrated model error. The independent result is alpha_cutoff_reference.json.

## Numerical evidence and limits

The reproducible 84-case matrix in validation/alpha_grid_convergence.csv
checks A=8.829,9.3,12 MeV and modes 1,2,3,13, 32/64/128/256 quadrature points
per dimension, and 1/2/4/10 keV cutoffs. Histogram cells cover 0..8 MeV in
32 intervals; below/above spill is included in the L1 comparison.
The largest normalized histogram L1 change from 128 to 256 is 0.9245%; the
largest raw normalization change is 7.12e-10. The largest cutoff histogram
change is 9.91e-18. Particle/energy input residuals round to zero at the
reported double precision. No inference of zero physical or model error
follows from those residuals. New grids or parameter combinations require
their own resolution check.

Reproduce after a CMake build with BUILD_TESTING:

```
./build/validate_alpha_grid > /tmp/alpha_grid_convergence.csv
python3 tools/analyze_alpha_grid.py /tmp/alpha_grid_convergence.csv \
  --output /tmp/alpha_grid_convergence_summary.json
```

Tests additionally cover strict-API parity, retaining nonzero other
permutations at even-l endpoints, mixture endpoints, phase-dependent
interference, two-sided spill and error clearing. See the compiler validation
record for C/Fortran builds and installed C11 use.

## Nuclear model status

Kuhlwein Eq. (4) supplies a coherent mixture form, but the checked source
text does not fix the separate phase-space normalization convention for its
fitted percentages. The diagnostic k=.76 and delta=.67*2pi are **examples in
the explicit convention above**, not a claim that the author's fit has been
reproduced. See KUHLWEIN_MIXTURE_CONVENTION.md. The source-algebra discrepancy
in ANGULAR_CORRELATION_AUDIT.md remains open for physical validation.

This layer selects neither alpha0/alpha1 branching, nor an incident-energy
resonance mixture, nor final-state interaction corrections. It does not
provide the beam-axis angular correlations required for general boosted
laboratory spectra. Those choices and experimental comparisons are required
before it becomes the full pB production source. No complete BALDUR/EXL run
is validated by these numerical tests.

## Published-shape and resolution check

See PUBLISHED_SPECTRUM_CHECK.md and validation/published-spectrum/. For60keV
output bins,128x128 quadrature produces visible shape aliasing despite exact
particle/energy closure;256 to512 changes normalized bin shape by about0.4%.
Do not infer spectral convergence from moment closure. The published l3 and
mixture curves are not yet reproduced to the same accuracy as pure l1; the
paper's fitted coefficient is not certified in this API's unit-basis convention.
