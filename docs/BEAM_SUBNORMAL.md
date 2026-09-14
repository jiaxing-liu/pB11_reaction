# Beam reference energy moments at IEEE underflow (D081)

The canonical D-on-T model at projectile483.993402keV and target10keV failed
with status4 in beam_segment's public energy-identity check. The above-fit
segment has a dimensionless quadrature rate1.3323724476802233e-289. Conversion
to SI produces subnormal energy moments; their independently rounded sum has
relative residual -.2 despite a converged normal-valued quadrature integral.
The old single-point reference returned no source; its failure is preserved in
the coordinator fast-cost/ archive. This also appeared at call236 of an evolved
fast-distribution controller fixture (the original whole fixture still fails a
separate charged-spill guard after this repair).

Retain the original relative condition1e-9*sum(abs energy moments), adding only
2*double denorm_min as an absolute representational bound for four independently
rounded public double moments (at most half denorm_min each). The coefficient
and moment calculations are unchanged; no tail is zeroed and no moment is
reconstructed from another. No nuclear data, quadrature tolerance or coupled
charged-spill guard changes. The existing put() contract already permits IEEE
underflow. This is scoped to beam_segment, not a generic tolerance relaxation.

The scalar corrected call returns status0 and rate4.4914501048856235e-22m3/s.
The regression asserts a positive subnormal above-fit rate, confirms the old
relative-only gate would fail, and checks the new measured rounding bound.
Selected rate-model, beam, beam-birth and coupled-fast suites pass4/4 on both
GNU and Intel/r8 build trees; full Intel BALDUR rebuild passes. Intel/r8 denotes
Fortran/host configuration; these C++ objects use the configured C++ compiler.

The evolved source fixture now reaches an independent charged alpha spill:
fastD4538.1174509712828keV yields above-grid coefficient
1.9733775364326608e-316m3/s on its25MeV kinetic grid. That rejection is retained.
Do not label this numeric repair complete dense-distribution integration or
relax the spill condition without an explicit physical/representation contract.
