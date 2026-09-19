# Finite static magnetic orbit driver (v1)

`fusion_c_prompt_orbit` composes the magnetic trial push with caller-owned
point/whole-segment domain callbacks. It is host-independent Cartesian SI code;
no BALDUR indexing, geometry arrays, allocation, mutable physics state or I/O.
See `include/fusion_prompt_orbit.h` for the C ABI. Installed with the library.

This is an empirical numerical model-boundary event/retention operator, not a
physical loss model. E=0, static magnetic field, no collisions. Initial proper
velocity is u=gamma*v. All refinement levels restart the same initial state and
horizon. The number of uniform steps doubles; rounded adjacent horizon fractions
define each actual dt. There are explicit finite step and refinement limits.

The caller must supply all tolerances before running and deterministic static
callbacks. Counters can change, but callback physics cannot depend on trial
history. The point callback uses the boundary classification enums. The segment
callback has the full D220 contract: CLEAR proves the complete segment including
endpoints, and CANDIDATE proves the inside prefix and an outside witness.
Endpoint point checks detect contradictions, but cannot prove an arbitrary
callback's whole-segment claim. Use the validated boundary implementation where
its geometric applicability conditions hold. Every new geometry is responsible
for those conditions and numeric allowances.

No field query occurs before the first half-drift is certified CLEAR and its
endpoints verified INSIDE. Candidates stop that level without querying an
outside field. The requested segment fraction tolerance allocates one quarter
of the time budget internally; actual geometric width is checked independently,
including its representation as absolute double times. The driver does not
reinterpret the callback's internal tolerance flag.

Return code zero permits three outcomes:

- UNRESOLVED (0): only reason/counters are usable. Never count this as loss or
  retention. Includes ambiguity, resources, time representation and refinement.
- RETAINED (1): two adjacent completed full-window paths satisfy caller final
  position and relative-u difference tolerances, after at least min_levels.
- EVENT (2): two adjacent localized candidates satisfy caller time/position
  difference tolerances. Returned bracket is the finer DISCRETE split path.

Time spread is center difference plus both geometric half-widths. Position
spread is the maximum over four endpoint-pair distances; overlapping candidate
intervals cannot falsely create zero estimated spread. Relative-u difference
uses the larger norm of the two final proper velocities, with both zero giving
zero. No fitted order, extrapolation, or true-error certificate is implied.

EVENT does not provide exit momentum: phase1 has no completed rotation and
phase2 rotated u also requires separate validation as physical crossing momentum.
Energy for this E=0 operator is the initial invariant. Do not use this restriction
to silently substitute a collision-inclusive or time-dependent loss model.
No wall/SOL/vacuum continuation or re-entry is implemented. A birth ensemble must
retain unresolved weight separately; outcome2 is not yet an ensemble loss rate.

Callbacks must not throw. Nonzero callback status propagates unchanged; malformed
classifications or contradictory endpoint claims return1004. Every nonzero API
status clears the result. Zero-initialized outcome is UNRESOLVED, not success.
Input state is immutable and input/output memory must not overlap.

Tests use independent analytic convex-sphere straight paths and a uniform-B
circular trajectory, with predefined absolute tolerances. Other checks cover
minimum levels, phase1 without field access, deterministic replay, ambiguous
initial state, callback failures/false CLEAR, resource limits and deliberately
unattainable tolerances. The uniform-B test checks true crossing time and position,
not just the same-method empirical spread. Independent underlying pusher,
boundary, field and inverse tests remain required. This first version makes no
formal adaptive/global error guarantee or complete confinement claim.
