# Product kinematics and conservative birth mapping

This increment provides reusable kinematics and grid projection, not a
validated pB birth-spectrum model. C ABI: `fusion_products.h`; Fortran:
`fusion_products_fortran`. Energies are J, masses kg, momenta kg m/s,
particle rates m^-3 s^-1 and energy rates W/m^3. Both APIs are stateless.

## Sequential equal-mass kinematics

`fusion_c_three_equal_sequential_cm(m, A, q, cosine, out)` represents a
parent at rest with total energy 3mc^2+A, decaying to one product and an
intermediate of rest energy 2mc^2+q. That intermediate decays to the other
two equal-mass products. Inputs require m>0, A>=0, 0<=q<=A and |cosine|<=1.
The primary moves along +z; the supplied cosine is measured relative to the
intermediate motion (-z) in its own rest frame. An x-z decay plane is chosen.
All three CM kinetic energies and momenta are returned.

Writing r=mc^2, B=2r+q, W=3r+A and d=A-q, two-body conservation gives
T_primary=d(d+2B)/(2W), T_recoil=d(d+2r)/(2W),
p_primary*c=sqrt(T_primary*(T_primary+2r)), and
p_secondary_star*c=sqrt((q/2)*(q/2+2r)). A Lorentz boost by the intermediate
recoil produces the secondary momenta. Returned kinetic energies are
computed using pc squared divided by sqrt(r squared+pc squared)+r; this
avoids subtracting nearly equal energies near a stationary secondary.

This function chooses no intermediate-energy probability, angular density,
branching, nuclear Q, spin alignment, coherent amplitude or permutation
weight. It performs no CM-to-laboratory transformation. In particular,
A=Q+relative kinetic energy does not include the reactant CM kinetic energy;
a full laboratory source must account for that separately. Relativistic
kinematics do not make the existing classical collision model relativistic.

## Birth packets to the FP grid

`fusion_c_map_birth_packets` takes arbitrary nonnegative packet rates at
specified energies and projects them onto arithmetic centers of increasing
nonnegative energy edges. Adjacent-center linear interpolation preserves
both number and kinetic energy. Cell outputs are integrated birth rates,
not densities per unit energy; they match the existing FP source convention.

The interpolation domain is the hull of the cell centers, not the edges.
Packets below/above that hull remain in explicit below/above particle and
energy ledgers. They are not clipped, discarded, or automatically classified
as thermal ash or escape. A successful status includes computation of these
spill ledgers; every caller must account for them. Centers are evaluated in
extended precision, so a double-precision packet nominally at an endpoint
may lie just outside its exact arithmetic center. Grid design should avoid
relying on endpoint equality. Empty packet arrays are supported.

Finite nonnegative inputs are required. Errors clear available outputs.
Caller arrays must not overlap. C++ catches allocation exceptions at the ABI
boundary; no accepted state, RNG, accumulated counters or I/O is hidden.

## Validation and limits

`test_fusion_products.cpp` independently checks packet inventories, explicit
spill energy, empty/one-cell cases and invalid inputs. Forty-five sequential
configurations span intermediate energies from zero to all available energy,
including the near-stationary-secondary region, and collinear/general angles.
Checks reconstruct each mass shell, total energy/momentum and the intermediate
invariant mass from returned vectors; reversing the cosine exchanges the
secondary energies.

A composition test prescribes five equally weighted angular events at
q=0.3A with A=8.68 MeV, creating three alphas per event. This is a synthetic
source, not a fitted pB spectrum. It maps births to 240 logarithmic cells
from 0.01 to 10000 keV and advances 100 steps of 0.002 s in fixed 5 keV
electron/deuteron baths at 1e20 m^-3 and lnLambda=15. An escape frequency
0.3 s^-1 and first-cell removal frequency 1000 s^-1 are explicit bookkeeping
test policies, not validated device loss or ash criteria. Independently summed
final inventory, bath heat, escaped and transferred inventories close the
continuous birth budget. The observed relative number and energy residuals
are 4.7187e-17 and -8.83158e-18. See validation/m3-products-* for compiler and
ABI evidence. No BALDUR stateful coupling or EXL-50U prediction is established
by these tests.
