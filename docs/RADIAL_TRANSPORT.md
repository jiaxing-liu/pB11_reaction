# Independent conservative radial transport

This additive stateless numerical operator transports explicitly supplied
species/energy-bin components. It does not prescribe a tokamak fast-ion
transport coefficient, orbit-loss law, adiabatic work, or host geometry policy.
The caller supplies old/new cell volumes and face coefficients in SI units.
BALDUR-specific indexing, metric extraction and acceptance remain in the host.

## Discrete contract

For component p and outward face j, define

    L_j = max(a_j,0) + k_pj
    R_j = min(a_j,0) - k_pj
    F_j = L_j n_left(new) + R_j n_right(new)
    V_new,j n_new,j - V_old,j n_old,j = dt (F_j - F_(j+1)).

Here a is a signed volume rate in m3/s and k is a nonnegative diffusive
conductance in m3/s. Component arrays are component-major in C; Fortran
uses (zones,components), faces (zones+1,components), boundaries (2,components).
Boundary reservoirs supply the missing density at the two endpoint faces.
Zero a and k give a reflecting boundary. Returned face amounts dt*F are signed.
The inner ledger is positive inward; the outer ledger is positive outward.
Both signs may reverse. Initial and final inventories use their actual volumes.

Backward Euler with upwind advection gives a tridiagonal matrix with
nonpositive off-diagonals. Each column has strictly positive excess V_new,
with additional nonnegative boundary outflow terms, so this is a nonsingular
M-matrix for finite positive volumes. A long-double Thomas solve checks all
pivots and the final nonnegative representable result. Every component must
pass before any candidate output is published. Numerical overflow and invalid
inputs are reported, not converted to a successful zero solution.

The method is first order in time and uses upwind numerical diffusion. It is
not a convergence result for any production fast-ion trajectory. Extended
precision zone and total inventory checks use 1e-10 relative scale plus the
measured conversion error propagated through the discrete flux operator.
IEEE underflow may round trace values to zero; no population floor is applied.

## Geometry and energy boundary

For mesh swept-volume rate G, pure coordinate advection has a=-G. Uniform
physical density is preserved when V_new-V_old=dt*(G_outer-G_inner).
Volume change with zero face flux instead conserves each cell inventory and
changes its density. These are different tests. Energy-bin identities are
unchanged: weighting all components and their face amounts by their fixed bin
energy gives the corresponding transported energy. Adiabatic redistribution
between energy bins and signed work still require a separate operator and
ledger. Neither can be inferred from the density dilution operation.

## Verification and remaining integration

The native C++ tests cover moving-coordinate uniform-density invariance,
closed compression dilution, exact two-cell/two-component diffusion, time
refinement against an exponential exchange solution, signed inflow/outflow,
prescribed diffusive reservoirs, repeatability, failure atomicity, overflow,
and measured trace rounding. Fortran tests additionally check ABI layout and
rank/extent rejection. The C header and exported Fortran target are installed.

This operator is not yet called by the production BALDUR driver. Source-state
spatial exchange/work accounting, retry-safe geometry epochs, adiabatic energy
transport and complete coupled moving-grid windows remain required. See the
coordinator HOST_GEOMETRY_EPOCH.md for measured host epoch discontinuities.
