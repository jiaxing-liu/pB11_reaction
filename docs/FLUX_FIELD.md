# Fourier flux-coordinate field evaluator (D210)

`fusion_c_flux_field_eval` is a host-independent, stateless C ABI implemented in
C++. It evaluates an up/down symmetric axisymmetric snapshot in `(xi,theta)`.
The header fixes array ordering, SI units, interpolation and error semantics.
Geometry coefficients use piecewise cubic Hermite with supplied native slopes;
all geometric derivatives come from that same representation. `F=RBphi` and
`C` are explicitly piecewise linear flux functions. Their interpolation is a
model-resolution choice requiring refinement, not reconstruction of unknown
sub-grid physics. No final field-component smoothing/projection is applied.

Signed `J=R*(Ztheta*Rxi-Rtheta*Zxi)` is retained. The caller owns coordinate and
physical-direction conventions. Native FFT array order is not part of the API.
No default loss coefficient, orbit evolution or dependency on BALDUR is added.

The current phase deliberately has no `(R,Z)` inverse lookup, magnetic-axis
regularization, global surface topology certificate or orbit integrator. Exact
xi=0 returns `FUSION_FIELD_AXIS_COORDINATE_SINGULAR` (1001); it is a coordinate
limitation, never a particle-loss event. Locally degenerate geometry returns a
numerical error. Outside supplied radial range returns OUT_OF_RANGE without
extrapolation. Grid-domain exit is not automatically LCFS exit: boundary identity
must be separately established. Revision1 fixes Hermite geometry / linear F,C /
analytic Fourier policy; do not silently change it. A returned success at one point cannot certify global nesting.
Full-domain handling is required before production particle-orbit use.

Validation includes a circular analytic solution down to xi=1e-8, a cubic radial
shape and derivative, periodicity, signed-J reflection, cylindrical divergence
using the analytic inverse circle, and error/output-clearing contracts. These
are in `tests/test_fusion_flux_field.cpp`. CMake/CTest and ASan/UBSan pass.
The circular divergence test is near roundoff; it is not a general shaped-field
finite-difference convergence study or an axis-limit test at xi=0.

Host-side D209 accepted-initial snapshot comparisons add960native points and
5760interior samples: maximum scaled field/geometry discrepancy1.04e-15,
minimum sampled signed J0.00815m^3, tangency numerator<=1.39e-17T m.
Only sampled local regularity is established, not a proof of global invertibility.
The host repository archives those case-dependent observations; no host fixture
or file I/O enters this library.

Build and run:

```sh
cmake -S . -B /tmp/pb11-field-build -DPB11_BUILD_FORTRAN=OFF -DBUILD_TESTING=ON
cmake --build /tmp/pb11-field-build --target test_fusion_flux_field
ctest --test-dir /tmp/pb11-field-build -R '^fusion_flux_field_tests$' --output-on-failure
```

The existing cache identity hashes all library sources/headers. Adding this
independent module therefore changes identities in newly configured builds.
Frozen running executables and their caches remain untouched; never relabel old
cache identity bytes to force reuse. No host production binary was rebuilt for
D210. The new module does not alter nuclear physics or old ABI entry points.
