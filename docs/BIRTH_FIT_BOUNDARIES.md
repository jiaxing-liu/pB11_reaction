# Source quadrature at piecewise nuclear fit boundaries (D087)

The actual800-cell DD source trial fails at zone19 after10us of requested20us.
D086diagnostic22062 first located the rate/debit consistency gate. D08722832
then captured DTchannel3, triton projectile slot1, flattened800cell index1892,
E=2.1399859254691007e-13J, targetkT=1.4401072590598843e-17J, relativemax=
4.0054415850000003e-13J. Old accepted fast inventory=1.5071983349875445e-212m^-3.
Rate/debit discrepancies1.0422590944125288e-5/1.0422590943859333e-5 exceed the
unchanged1e-5 gates. Both diagnostics are terminal1, not complete windows.

The nuclear cross-section implementation switches rational fits at DT530keV
and DHe3900keV relative energy. Independent beam reference make_cuts already
splits these boundaries. Birth relative_knots omitted them, allowing a single
fixed-order quadrature interval to straddle the piecewise fit change. Add these
two physical integration splits when they lie below the requested upper limit.
This shared helper serves thermal/pair and beam birth. Nuclear formulas,
continuation, relative upper limit, node orders, tolerance and host gates are
unchanged. No tiny accepted inventory is used as a reason to ignore error.

A permanent test uses the exact captured DT input, plus D-on-T and D-on-He3
inputs at and +/-0.5% around their respective fit boundaries. All source rate
and reactant-debit discrepancies must be below1e-8, and complete outgoing
number/energy (including explicit spill) must close. The source-rate check is
independent of outgoing mapper grid; this scalar regression uses the existing
16cell mapping fixture and does not itself establish800cell host acceptance.
The old library fails the new test at the exact DT gate; the split version
passes. GNU6/6 and Intel/r8 6/6 birth/beam-birth/coupled-fast/thermal/table tests
pass, and the full host build succeeds. Intel/r8 uses GNU C++ with Intel Fortran.

Fresh matched production runs use this same code:400cells25970 and800cells39265,
20us requested,DDthermal+DTfast,modeltags87002/87003, movinggeometry,handoff0.
Outputs /tmp/dd-fast-knots400-d087 and /tmp/dd-fast-knots800-d087; cases with
same names. Neither is accepted until its complete output is checked.
Original110us400cell42952 and1600cellDT1s1244 continue with frozen executables;
last observed prefixes60us/5steps and.13316125s/21306steps are not final results.
A bounded Luna exactpoint/order scan is pending at /tmp/beam-gate-scan-d087;
use its eventual report as additional evidence, not as a prerequisite already
assumed passed. No full-goal/convergence/ash/loss-completion claim is made.
