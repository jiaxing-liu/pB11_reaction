# Geometry evolution of both kinetic components

`fusion_c_kinetic_geometry_trial` composes the existing independent radial
finite-volume operator and signed energy-work operator, in that order. It is
stateless: no nuclear source, thermal bath update, accepted-state commit,
BALDUR index, or implicit loss/boundary model is introduced. Its Fortran binding
uses explicit ISO_C_BINDING kinds.

## Layout and normalization

- S/T density arrays: C `[zone][species][cell]`, Fortran `(cells,6,zones)`.
- Radial conductance: C `[species][cell][face]`, Fortran `(zones+1,cells,6)`;
  nonnegative m3/s, shared by S/T at a given species and energy.
- Separate S/T boundary populations: C `[species][cell][inner_or_outer]`,
  Fortran `(2,cells,6)`, physical m^-3.
- Actual old/new volumes Va/Vb: m3. Time: s. Grid edges: J.
- Outward-positive advection: m3/s; moving coordinates supply minus their
  swept-volume rate. Compression C: s^-1; positive cools.

The caller supplies all coefficients and boundaries. Reflecting faces require
zero advection and conductance; vacuum boundary densities alone do not create a
reflecting boundary. Prescribed populations are not inferred from three scalar
moments. The C header states dimension/resource limits and the no-overlap
contract. All valid-size outputs clear on failure; results are published only
after all zones/species/components succeed.

The operator transposes populations into the radial kernel's component-major
layout, advances from Va to Vb, then applies energy drift at Vb. Shared face
amounts produce signed extensive spatial number and midpoint-energy exchange.
S/T contributions are combined per species. Work and numerical domain outflow
are multiplied by Vb before returning one `fusion_transport_ledger_v1` per zone.
The result can be passed directly to the volume-aware source-state stage.

If nuclear/collision evolution precedes this call at Va, retain that operator's
source-density ledger at Va. Do not reinterpret it as a Vb density increment or
as a midpoint-volume source. The host may normalize its source RHS at a third
explicit volume Vm. Geometry work does not become bath heat; numerical grid
outflow does not become physical escape or thermal ash.

## Numerical scope and tests

This is first-order radial-then-energy operator splitting, using the existing
implicit upwind kernels. It does not provide Strang splitting, a higher-order
transport scheme, a physical diffusivity/orbit model, or beyond-grid kinetic
re-entry. The resolution/diffusion limitations in RADIAL_TRANSPORT.md and
ENERGY_WORK.md still apply. Zero compression or zero population skips an exact
identity energy step only after common domain inputs are validated.

The C++ test compares all six species and both components with explicit
single-component calls to the underlying kernels. Independent physical-density
moments verify extensive N/U budgets and direct multizone v3 state acceptance.
An energy-unit covariance check repeats with 1e-16 J scaling. Identity,
uniform-density geometric conservation with explicit boundary populations,
invalid late inputs, and a genuine final-zone energy overflow after successful
radial transport are covered. No partial candidate escapes on that failure.
The Fortran test checks layout, a nonfirst species in S/T, identity/dilution,
and malformed-array clearing.

GNU/default and Intel-r8 each pass the two focused tests. Independent installed
C11 and Fortran consumers pass through the installed CMake target export.
Evidence is in `docs/validation/kinetic-geometry/`. This does not certify a
complete moving-geometry BALDUR nuclear run; the host controller/driver and
whole-window tests remain separate requirements.
