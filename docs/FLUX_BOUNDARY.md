# Guarded model-LCFS classifier, revision 1

`fusion_c_flux_boundary` is an independent C ABI for the finite Fourier curve
R=R0+sum Rj cos(j theta), Z=sum Zj sin(j theta). It takes one surface, not a host
snapshot or inverse map. Lengths are metres. The header specifies sizes, units,
status codes and explicit caller error allowances. The model surface must be
chosen by the caller; this API has no material-wall information.

The implementation uses long double and an equally spaced polygon. Define
P=sum max(abs(Rj),abs(Zj)), M=sum j*j*max(abs(Rj),abs(Zj)), h=2*pi/N.
P bounds the distance from (R0,0); M bounds the second theta derivative norm.
The maximum curve-to-corresponding-chord deviation is bounded by M*h*h/8 in
exact arithmetic. W=(R-R0)*Ztheta-Z*Rtheta and abs(Wprime)<=P*M.
Only geometry passing all the following checks is supported:

- Positive R using R0-sum abs(Rj)>geometry allowance.
- min sampled W minus P*M*h/2 minus angular allowance is strictly positive.
- Center-to-polygon distance exceeds the curve/chord bound plus geometry allowance.
- Polygon winding about the center is exactly +1.

Under the declared error assumptions these conditions imply positive angular
progress, a homotopy avoiding the center and one winding. They are a sufficient
star-shaped geometry check, not a necessary condition for a valid plasma surface.
Clockwise, degenerate, multiply wound and unsupported shaped curves return
UNSUPPORTED_GEOMETRY, never OUTSIDE. The code does not use outward-rounded
interval arithmetic: this is a guarded numerical check, not a machine-verified
mathematical certificate.

The caller supplies positive aggregate geometry and angular error allowances.
The geometry allowance must cover evaluation, distance and bound arithmetic and
query uncertainty. The angular allowance must cover W and P*M*h/2 evaluation and
bound arithmetic. Input coefficient uncertainty, if present, must also be covered.
The library does not infer these from machine epsilon or validate their adequacy.
Tests use 1e-12m and1e-12m^2 for metre-scale fixtures on the tested compiler/libm;
these are explicit numerical test choices, not device error bars or a universal
accuracy guarantee. Extremely scaled inputs require new allowances/validation.

Query classification uses winding and distance from the polygon. Distance no
greater than M*h*h/8+geometry allowance gives AMBIGUOUS. Otherwise winding1 is
INSIDE and winding0 is OUTSIDE; incompatible winding gives unsupported geometry.
No inverse-convergence status is interpreted as location. No output is a wall
hit or a physical escape coefficient. Caller refinement may reduce the ambiguity
band; the API never arbitrarily resolves an uncertain boundary point.

Validation includes analytic circle/ellipse inside/on/outside queries, Fourier
fifth-harmonic curve samples, quadratic ambiguity-band refinement, multiply wound,
clockwise and degenerate rejection, invalid arguments and cleared error outputs.
The coordinating D217 report adds actual D209 LCFS samples, interior flux surfaces
and outward normal geometric probes without field extrapolation. CTest and
ASan/UBSan cover the implementation. No host production binary is changed.
