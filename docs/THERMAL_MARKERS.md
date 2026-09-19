# Thermal mapped-bin tensor markers

`fusion_c_thermal_markers` builds deterministic Cartesian proper-velocity markers
for one accepted, extensive, charged thermal source bin. The caller supplies the
canonical nuclear mass, bin-center kinetic energy, positive physical-volume
nodes and positive solid-angle probability nodes. There is no BALDUR dependency.

Model1 means zero-drift one-particle isotropy, uniform birth per supplied physical
volume, and mapped-bin-center energy. These are explicit modeling assumptions;
a zone-total source does not determine a true poloidal localization. The API
cannot infer that its input was accepted or that the mass matches the species.
The adapter owns those identities and the corresponding immutable field epoch.

For each spatial/directional pair the particle weight is
`N * (volume_weight / source_volume) * direction_weight` and the energy weight is
that weight times K. N is already extensive; no extra volume factor is applied.
Direction is LAB Cartesian. Full solid-angle quadrature includes gyrophase.
Proper speed is `sqrt(K/m)*sqrt(2+K/(m*c*c))`, with c=299792458m/s. No vectors or
weights are silently renormalized. The scalar tolerance is caller-selected and
applies to volume/normalization, unit directions and first/second angular moments.
Moment consistency does not prove isotropic integral convergence.

Counts and product are limited to1e6. Output capacity equals the product.
N=0 returns the same fixed tensor shape with zero particle weights (valid
positions/velocities); callers can omit zero bins before invoking the builder.
All outputs clear before input validation for valid declared capacity. A local
candidate vector is copied out only after count/energy closure and finite double
representation checks. Positive underflow fails instead of discarding particles.
Memory allocation exceptions return exception status; no OOM injection is claimed.

Below/above tails retain their original unsampled N/E. Their mean energy is not
sufficient to build representative trajectories. Beam angular marginals are not
available and must not be made isotropic. Neutrons are noncharged/not-applicable,
not unresolved charged orbits. This API does not reconstruct pB three-alpha event
correlations. It neither traces orbits nor applies production source losses.

The unit test covers six-axis, rotated cube and nonuniform mixed angular rules,
unequal spatial volumes, proper-energy reconstruction at1e-30..1e-7J, deterministic
replay, reducer accounting, malformed measures/directions, positive underflow,
zero number/energy, nonfinite positions and excessive count products. CTest and
standalone ASan/UBSan pass; C11 inclusion compiles. These validate construction,
not physical confinement or the full quadrature/driver composition.

Next validation must separately refine spatial, angular and orbit discretization
against independent analytic boundaries, maintain unsampled accounts and bind
actual accepted source/field epochs before any production loss interpretation.
