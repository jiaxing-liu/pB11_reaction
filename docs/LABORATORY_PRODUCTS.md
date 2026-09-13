# Laboratory product kinematics

`fusion_laboratory.h` defines massive particle four-vectors as mass (kg),
kinetic energy (J) and Cartesian momentum (kg m/s). No nuclear Q, species
identity, branching or angular probability is inferred. The API is stateless.

## Two-body center-of-momentum frame

For product rest energies r_a=m_a c^2, r_b=m_b c^2 and available kinetic
energy A, W=r_a+r_b+A. The first product has
T_a=A(A+2r_b)/(2W), momentum magnitude sqrt(T_a(T_a+2r_a))/c;
the second has exactly opposite momentum. Returned energies are reconstructed
from the returned rounded momenta using the mass-shell relation. The caller
supplies the direction of product zero; it must be unit within 1e-12 and is
normalized internally. Both masses must be positive and A nonnegative.

This supports the kinematics needed for DT, DD and DHe3 without selecting a
nuclear angular model. Neutrons must be recorded as escaping nuclear energy
by the eventual source layer, not automatically deposited into charged baths.

## Active Lorentz boost

`fusion_c_boost_particles` accepts already oriented input four-vectors and
an explicit boost velocity V with |V|<c. The sign convention is active:
an initially resting particle acquires velocity V. For total particle energy
E=mc^2+T, E'=gamma(E+V dot p),
p'=p+[(gamma-1)(V dot p)/V^2+gamma E/c^2] V.
The implementation evaluates (gamma-1)/V^2 as
gamma^2/((gamma+1)c^2), including V=0, and reconstructs T' from p' to avoid
negative cancellation near a stopped product. Input mass shells must agree
within 1e-10 relative to kinetic energy. No negative state is clipped.

The ledger independently transforms summed input four-momentum and compares
it with the sum of returned particle states. Energies exclude rest energy.
Residual acceptance uses the magnitude of individual transformed terms;
near a fully stopped state this remains meaningful when final kinetic energy
is near zero. Expected energy may carry a tiny signed rounding residual in
such a cancellation limit; it must not be reinterpreted as physical heating.
Arrays must not overlap and count>=1. Errors clear available outputs.
The caller owns orientation sampling, birth projection and all spill ledgers.

## Independent tests and model boundary

`test_fusion_laboratory.cpp` checks unequal-mass two-body energy/momentum and
mass shells at zero, 1 keV and 18 MeV available energy. Fifteen boosts span
zero to 0.8c along three axes. Tests independently check total invariant mass,
total lab energy for CM events and inverse recovery of individual momenta
and kinetic energies. A rotated three-alpha event with non-collinear boost,
classical small-speed limit, near-stationary deboost and invalid mass-shell/
velocity/direction cases are included. High-speed cases validate kinematics
only, not the nonrelativistic collision or Maxwellian reaction models.

The complete reaction source must adopt consistent nuclear masses and Q,
account for reaction-selected reactant energy, and specify a valid product
probability model. This increment does not supply those policies or activate
BALDUR stateful kinetics. See PRODUCT_KINEMATICS_AND_MAPPING.md for the
existing CM event and conservative birth-grid interfaces.
