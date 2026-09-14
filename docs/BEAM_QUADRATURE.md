# Fixed-coordinate beam quadrature (D083)

The actual DD-host secondary04 failure was captured in the beam source call:
DT channel3, incoming triton slot1, projectile1.1768500517434364e-14J,
target kT8.6808179066013538e-18J. Its independent beam reference returned4.
A standalone exact-input reproducer confirms this. The below-fit segment had
dimensionless integral3.8706120470433687e-317 and reported error
5.1876892813330887e-321; the relative-error gate rejected the model. The fit
segment itself had normal rate and energy identity agreement6.9e-17relative.
This is distinct from the D081 public SI moment and D082 inventory-spill cases.

Each existing physical integration subinterval is now mapped to[-1,1], with
its half-width Jacobian inside the integrand. Physical cuts, kernel, quadrature
rule, depth15, tolerance1e-11, nuclear model and downstream gates are unchanged.
The installed Boost1.90 implementation computes an error estimate for its
normalized rule while scaling the integral estimate. Passing narrow intervals
directly can drive needless subdivision and an inconsistent relative check.
Explicit fixed-coordinate integration keeps the Jacobian in both numerical
quantities. Recursive subinterval estimates remain subject to the unchanged
library checks; this is no claim of a rigorous quadrature error bound.

The exact host point now returns status0 and totalK2.868440203337139e-22m3/s.
A32-point T-on-D scan (projectile0..1000keV, target.02/.1/1/10keV) preserves all
success statuses; max relative change in total rate4.440892098500626e-16.
Measured total wall time falls from6.412812607s to.044196243s on this run.
This is a representative scalar benchmark, not full-host/pB acceleration proof.
The original scan CSV's projectile_kT_keV heading denotes projectile kinetic
energy, not a projectile temperature; comparison.json calls it E_keV.

Exact host-input and near-zero-projectile continuity regressions are added.
Existing subnormal regression is retained. GNU5/5 and Intel/r8 5/5 selected
rate-model/beam/beam-birth/coupled-fast/table suites and full host build pass.
The pinned host-point value is a numerical regression value; independent angular
and energy checks are provided by the existing tests, not by that pin alone.
A fresh actual DD secondary05 run is pending; no complete host window or energy
convergence follows solely from this scalar repair. See coordinator
BEAM_QUADRATURE.md and beam-quadrature/ for evidence and live handle.
