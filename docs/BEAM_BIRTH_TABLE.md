# Fixed-projectile-energy full beam birth tables

D089 adds `fusion_beam_birth_table.h` and the optional ISO_C_BINDING module
`fusion_beam_birth_table_fortran`. The exported CMake Fortran target is
`pb11::fusion_beam_birth_table_fortran`. Existing APIs are unchanged.

Each immutable, caller-owned table fixes the channel, projectile slot, projectile
kinetic energy, source options and output energy edges. The sole interpolation
coordinate is positive target kT in joules. Construction copies inputs; evaluation
neither advances state nor accesses host data. Outputs use the existing complete
seven-product grid and birth coefficients, including both reactant energy debits
and number/energy below and above the grid. Multiply coefficients by distinct
beam and thermal densities outside this API; there is no identical-pool factor.

Adaptive geometric quarter/mid/three-quarter samples compare the interpolated
source with direct beam integration. The gates cover rate, both debits and the
maximum per-species number and energy L1 discrepancies, including spills. A
single convex log-temperature weight interpolates every quantity. These sampled
gates do not prove a uniform error bound; independent temperatures and coupled
convergence must also be checked. Out-of-domain evaluation fails, without
extrapolation. Cold target calculations remain available through the direct API.

Nodes retain every nonzero double grid value, including subnormals; exact zeros
alone are omitted. This is sparse storage of the full spectrum, not a moment
surrogate. The public output remains dense. Cumulative stored entries across all
direct construction samples are limited to 12,500,000 in addition to the explicit
knot/evaluation/depth budgets. This bounds spectral payload work; metadata, vector
capacity and temporary dense buffers are additional. Failure returns no partial
table. There is no renormalization, weak-source cutoff, energy interpolation or
hidden persistent cache.

Tests cover all five reaction channels, zero projectile energy, exact endpoints,
six independent temperatures per configuration, full spectra and spills,
particle/Q-energy closure, copied inputs, rejected extrapolation, atomic budget
failure, and an 800-cell endpoint containing subnormal spectrum entries. Fortran
tests cover shape validation, output clearing and lifetime as well as evaluation.
GNU and Intel Fortran builds are supported; the C++ compiler in both validation
builds was GNU. Caller must destroy each handle exactly once, after readers finish.

A representative DT T-on-D table (400 cells, fixed projectile energy
2.1399859254691007e-13 J, target 0.05–0.2 keV) required 17 knots and 95 direct
samples, retaining 2693 final nonzero values. Construction took 4.66 s; 1000
queries took 0.0484 s versus 0.489 s for 10 direct calls on the same machine.
Independent sampled rate differences were below 5e-7. These are a single-source
benchmark, not BALDUR performance. Tests use explicitly recorded controls, which
are not new production defaults. BALDUR/coupled-fast table dispatch and its
validation remain subsequent work; current host fast reactions still use the
direct source with the existing exact within-trial cache.

D090 update: explicit coupled dispatch is now available; see COUPLED_SOURCES.md.
The source-level benchmark above retains its original limited scope.
