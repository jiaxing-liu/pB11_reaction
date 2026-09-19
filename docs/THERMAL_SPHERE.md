# B=0 analytic composition benchmark (D229)

`test_fusion_thermal_sphere` composes the actual thermal marker builder, prompt
orbit driver and reducer in a unit sphere with B=0. The sphere callback proves
inside chords by convexity and brackets boundary contacts by bisection. An
independent positive quadratic root gives each marker's exact straight-ray exit.
Any orbit API failure invalidates the test. Discrete classifications must match,
exact discrete event weight must lie in the classification interval, and N/E
must reconstruct their inputs. Monochromatic number/energy fractions agree.

A continuous uniform-volume sphere and isotropic directions give an independent
integral target. For any fixed direction, points retained after distance L form
the intersection of two radius-a spheres separated by L. Its volume is
`pi*(4*a+L)*(2*a-L)^2/12`; division by `4*pi*a^3/3` and subtraction from one gives
`F=3*L/(4*a)-L^3/(16*a^3)` for0<=L<=2a. At L=a, F=11/16. A source at radius a/2
has exit condition direction cosine>=-1/4 and hence fraction5/8.

Spatial quadrature uses Gauss-Legendre in s=(r/a)^3 and cos(theta), uniform phi;
dV=V*ds*dcos(theta)/2*dphi/(2*pi). Volume normalization and first moments are
checked; isotropic second moments and mean r² are reported against3a²/5.
Direction rules use Gauss-Legendre cosine and uniform azimuth. No source or rule
is fitted to its resulting exit fraction. Proton K=10keV and T=a/v; driver
initial dt=T/2, levels2..4, max64steps, geometry time budget1e-9T, event tolerance
1e-7T/1e-7m, retained tolerance1e-8m/1e-10relative proper velocity.

Results from the current execution (all discrete classifications agree, zero
unresolved, maximum event-time error<4.48e-11T):

| Refinement | Fraction | Error relative to continuous target |
|---|---:|---:|
| Fixed radius, angular4 | .5000000000 | -.1250000000 |
| Fixed radius, angular8 | .6813418917 | .0563418917 |
| Fixed radius, angular16 | .5947253052 | -.0302746948 |
| Fixed radius, angular32 | .6430116038 | .0180116038 |
| Fixed radius, angular64 | .6205753705 | -.0044246295 |
| Spatial4, angular8 | .6899311156 | .0024311156 |
| Spatial6, angular8 | .6891467050 | .0016467050 |
| Spatial8, angular8 | .6869423691 | -.0005576309 |
| Spatial8, angular4 | .6874940264 | -.0000059736 |
| Spatial8, angular6 | .6876612619 | .0001612619 |

Spatial mean r² is.6005439071,.6001566695,.6000636121 at orders4/6/8. These
show convergence toward.6; angular indicator quadrature is not monotonic. The
particularly accurate coarse angular4 result is accidental and is not selected
as a production rule. Finest joint rule has524288markers and absolute fraction
error.000557631. This measures quadrature error; the classification interval does
NOT bound that continuous error. No claimed rigorous quadrature bound or assumed
Richardson order is attached to a discontinuous exit indicator.

This validates B=0 composition and separates discrete classification from
continuous sampling error. It does not test curved-orbit integration error,
physical wall loss, EXL50U loss, collisions, beam directions or unknown tails.
Uniform-B fixed-ensemble first-exit reference remains the next test.
