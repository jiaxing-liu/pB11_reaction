# Static magnetic test-particle push, revision 1

`fusion_c_magnetic_push` evolves synchronized Cartesian x and u=gamma*v using
half a spatial drift, a frozen-field magnetic rotation, then half a drift.
The field callback is called at the computed spatial midpoint, not at the initial
position. `fusion_c_magnetic_midpoint` exposes that request for caller preflight
checks. No host geometry or global state enters either API.

The rotation uses theta=q*|B|*dt/(m*gamma) and
u_new=u_parallel+(u-u_parallel)*cos(theta)+(u cross b)*sin(theta).
This is the exact frozen magnetic rotation used in Boris-C; see Zenitani and
Umeda, equations6,11,12, https://arxiv.org/html/1809.04378v2 . The present API
composes this with symmetric position drifts at synchronized times. It is a
second-order spatial trajectory approximation, not an exact spatial helix or
an exact orbit in a nonuniform field. It assumes static B, zero electric field,
no collisions and no feedback. Relativistic proper velocity is used; this does
not change the separate classical collision model's domain.

Mass is positive kg, signed charge is C, dt is signed seconds, x is metres,
u is m/s and callback B is tesla. u is NOT ordinary velocity or momentum.
No gyro-angle clipping is performed; caller must resolve gyration and spatial
variation and demonstrate step refinement. Large-angle rotation algebra can
remain finite even when the spatial orbit is inaccurate. Zero B and zero charge
preserve u; dt=0 returns the input. q=0/dt=0 skip the field callback. A nonnull
callback remains required. Negative dt supports reversibility checks.

Callback failures propagate unchanged; nonfinite field values fail. All available
outputs clear on errors; inputs and outputs must not overlap. The callback must
reject out-of-domain queries before evaluating the field and must not throw.
No allocations, hidden history or output files are introduced by the library.
The second drift uses the actual rounded returned u to calculate gamma.

This is a trial primitive. It does NOT validate either drift segment, initial
state, endpoint or boundary event, and does not provide a trajectory error bound.
A successful call is not an accepted particle path and never a loss event.
Boundary handling must retain unresolved approach separately from a validated
model-LCFS exit bracket. No automatic drift acceptance is implied by midpoint B.

Tests cover zero field/charge/time, charge sign, relativistic and low-speed
uniform-field rotation, oblique parallel component, dt reversal, small/large
rotation angles, invalid inputs and callback failure. Uniform-field position
errors decrease fourfold for doubled step counts while |u| stays near roundoff.
A coordinating independent RK4 reference audit also checks a smooth nonuniform,
divergence-free synthetic B; it is not a device-confinement validation.
