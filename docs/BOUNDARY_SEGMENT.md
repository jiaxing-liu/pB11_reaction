# Cartesian segment geometry events

`fusion_c_boundary_segment` classifies a straight Cartesian path against the
D217 Fourier model LCFS. It does not evaluate a field or integrate a trajectory.
Projection (x,y,z)->(sqrt(x*x+y*y),z) is nonexpansive but is generally NOT a
straight R/Z line, so the implementation does not substitute an R/Z chord.

A subsegment whose length is smaller than a safe inside endpoint's clearance
cannot reach the boundary. Clearance subtracts the polygon/curve band and caller
allowance. An extra geometry allowance covers segment arithmetic. Only whole
subsegments proven inside can be skipped. Fixed-stack depth-first subdivision
visits the left half first; reaching an uncertain early leaf stops the search,
so later evidence cannot erase an unresolved earlier encounter.

Outputs distinguish CLEAR (whole path covered), GEOMETRIC_CANDIDATE (safe inside
and outside endpoints with the preceding prefix covered), and UNRESOLVED (no
safe outside witness). A candidate can contain multiple crossings. It locates
the earliest possible contact, not a unique intersection or a physical exit.
The retained enclosing candidate survives ambiguous refinement. Its width may
exceed the caller request: tolerance_met and termination_reason report this.
Tolerance refers only to segment parameter fractions, independently of boundary
approximation error. No silently widened tolerance is labeled as achieved.

The API bounds traversal by caller interval budget and fixed depth52. Initial
AMBIGUOUS returns UNRESOLVED; initial OUTSIDE is OUT_OF_RANGE. Errors clear output.
It inherits the classifier's guarded, non-interval error assumptions, and the
caller geometry allowance must additionally cover Cartesian interpolation and
length arithmetic. No allocation, hidden state, host state or I/O is introduced.

Analytic torus fixtures cover safely inside segments, exit, inside-to-inside
paths through the central hole, exact boundary endpoints, uncertainty and budget
exhaustion. This must later be combined with integration-error and timestep
refinement evidence before any discrete drift crossing becomes an orbit event.

Outside-endpoint refinement can shrink a retained candidate after left-first
traversal stops. Its proven lower endpoint is the earliest unresolved leaf's
safe inside start, with all preceding subsegments covered. Sampling inside or
ambiguous locations advances only a search hint, NEVER the proven lower endpoint;
only safe outside samples shrink the upper endpoint. Thus nonmonotone out/in/out
classification cannot erase an earlier possible crossing. Tolerance uses the
returned bracket width, not the shrinking heuristic search interval. All extra
queries consume the same bounded work budget. This is not a monotone root solve.

## Inside-prefix refinement correction (2026-10-18)

A small unresolved leaf does not ensure that a retained safe-outside ancestor
candidate meets the requested fractional tolerance. Previously the traversal
could stop at the leaf limit and retain a needlessly wide ancestor; shrinking
only the outside endpoint could not advance its unproved inside prefix.
Now it continues left-first prefix proof in that case. The algorithm retains
safe inside/outside witnesses, no guessed ambiguous classifications, depth
and traversal caps, and the existing boundary/roundoff allowances.

The analytic circle regression now meets0.001 fractional tolerance;1e-5
remains below the geometric uncertainty and keeps tolerance_met=0. Four
historical unresolved trajectories become stable model-boundary events at
unchanged physical tolerances. A stricter min_levels7 rerun keeps all four
classifications, with maximum coarse/fine event bracket spread3.826e-13s
below the existing1e-12s requirement. This does not establish physical wall
losses, the full birth measure, horizon sufficiency or FP/prompt ownership.
Raw inputs/results and the geometry-only isolated prototype are archived in
validation/prompt-prefix-refinement-20261018/. The final combined library
passes all120 C++/Fortran tests; formal installed consumers also pass.
