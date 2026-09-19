# Fixed planar uniform-B ensemble (D230)

This benchmark isolates curved-orbit integration from birth quadrature. It uses
12 equally weighted proton10keV markers: initial radius.2/.5/.7m on +x, velocity
phase0/pi/2/pi/3pi/2 in the xy plane, and rho=.75m in uniform +z B. The sphere
has radius1m. Charge and canonical test mass are explicit in the fixture.
No isotropic ensemble or arbitrary-pitch reference is claimed.

With theta=Omega*t, x0=r, velocity phase phi, exact motion is
x=r+rho*(sin(phi)-sin(phi-theta)), y=rho*(cos(phi-theta)-cos(phi)).
The squared boundary residual is A+B*sin(theta)+C*cos(theta), with
A=r²+2*r*rho*sin(phi)+2*rho²-1, B=2*rho*r*cos(phi),
C=-2*r*rho*sin(phi)-2*rho². All derivative roots are atan2(B,C)+k*pi.
Enumerating them over the two tested horizons1 and2pi partitions the interval
into monotonic pieces. Left-to-right bisection locates the first crossing,
including a first exit followed by full-period return. This is an independent
analytic reference, not a check of endpoint signs over an entire orbit.

The actual prompt driver runs exactly two refinement levels, with intentionally
loose agreement tolerances to expose the discretization error at each selected
initial step. Initial Omega*dt=.2,.1,.05,.025,.0125. Actual step counts are
ceil(horizon/initial_dt) times the level refinement. Geometry budget1e-11/Omega;
max4096steps. Compare directly with true first-time/position/final-position, not
with the driver's adjacent-level empirical spread as a rigorous bound.

All levels classify the1rad case as4event/8retained and2pi as10event/2retained,
with zero unresolved. Weighted N/E closure passes. Maximum first-exit phase error
falls from6.5105e-4 to2.6007e-6 at1rad and1.5391e-3 to8.0819e-6 at2pi. Maximum
combined error decreases at each refinement; no fixed factor-four ratio is
required. Finest absolute position/phase gates are2e-5 in their unit scales.

A separate near-contact negative fixture has r=.25m,rho=.375-1e-8m,+y velocity,
full-period horizon. Its analytic maximum radius is1-2e-8m. Tight driver controls
allow retained or unresolved but prohibit a false event; actual outcome is
retained after3levels. This is not a proof for all grazing geometries.

The test validates curved planar numerical model-sphere exits and reentry
handling. It is not a physical wall-loss or EXL50U confinement result. Actual
accepted source/field epoch pairing, general3D application, beam/tail accounts,
collisions and production feedback remain separate requirements.
