# Extended-precision beam intermediates (D084)

Actual DD diagnostic secondary06 completed with wrapper1 (solver0 but source4,
zone32, only10us of20us). Captured T-on-D DT input:
E=1.4508757491041221e-14 J, target kT=1.1122299932400609e-17 J.
Independent reproduction at external58585e7 returns4. Its fit segment has
K=5.4616646690171559e-22 m3/s and relative energy residual-1.07e-16;
below-fit dimensionless rate=1.3254004049682707e-314 and error=
6.5908357155222289e-321 fail the unchanged relative-rate gate. Both SI
rate/error round to zero. Normalizing the integration coordinate (D083) did
not prevent double intermediate quantization at this second actual point.

Keep Gaussian probabilities, reaction integrands, quadrature coordinates,
integrals and error estimates in long double until SI output conversion.
Use the same 61-point Gauss-Kronrod rule, depth15, tolerance1e-11, physical
cuts, nuclear cross-section API and final relative-rate1e-8 gate. Public C ABI
remains double. D081 public energy rounding bound is retained. No tail is
explicitly zeroed, no error gate is bypassed, and no nuclear model changes.
This validation uses Linux extended-range long double; it does not establish
equivalent behavior on platforms where long double equals double.

Exact point now returns0 with K=5.4616646690171568e-22 m3/s. The32-point
T-on-D scan retains all success statuses, max total-rate relative difference
1.1102230246251565e-15 from58585e7. One measured scalar scan changes .044196243s
to .048419155s; this is not a full-host performance measurement. CSV heading
projectile_kT_keV is inherited and means projectile kinetic energy.

The actual DD secondary07 run uses the rebuilt production host (no logging
overlay), requested20us,400cells,moving geometry,handoff0, DD thermal channels
and DT fast channel,modeltag84001. It must reach its full requested window
and pass independent budgets before actual-host acceptance can be claimed.

Validation: GNU5/5 and Intel/r8 5/5 selected rate-model/beam/beam-birth/
coupled-fast/coupled-table tests pass, including the exact input and independent
3x3 +/-0.1% E/kT neighborhood. Existing regressions retained. Full host Intel
build succeeds. Both C++ test configurations use GNU C++; Intel/r8 denotes the
Fortran integration build, not icx validation. No full actual-host acceptance
is inferred from these bounded tests.
