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
