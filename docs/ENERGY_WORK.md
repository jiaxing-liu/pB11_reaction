# Signed energy-space work operator

`fusion_c_energy_work_trial` is a stateless conservative numerical kernel for
an explicitly supplied drift dE/dt=-C E. It preserves the existing scalar,
collision, radial and source-state ABIs. The Fortran binding exports the same
11-double/88-byte ledger and validates all extents before taking C addresses.

## Physics reference and derivation

Sokolov et al., JCP476 (2023)111923, Eq.(112)-(113), p.24, give the
angle-averaged phase-space measure and isotropic momentum compression term.
The original PDF equation was visually checked, including the p^3/3 variable.
This source studies energetic particles in the heliosphere, so it supports the
mathematical isotropic transport form, not tokamak orbit/loss validity.
[Primary paper](https://par.nsf.gov/servlets/purl/10493805),
[DOI](https://doi.org/10.1016/j.jcp.2023.111923).

With nonrelativistic E=p^2/(2m), the characteristic obeys

    dp/dt = -(div u) p/3
    dE/dt = -(2/3)(div u) E = -C E.

After converting from phase-space density to the integrated energy distribution
and separating physical-space advection/volume dilution, the remaining work
substep is

    partial_t g + partial_E(-C E g) = 0,
    N_j = integral_cell g dE.

This separation is the root's derivation, not a claim that the implementation
reproduces the paper's high-resolution Poisson-bracket algorithm. BALDUR's scalar
convrt compression coefficient provides the eventual caller-side reference;
its magnetic-coordinate divergence must be retained. C is frozen for one call.
The approximation is isotropic and nonrelativistic, with a prescribed scalar
compression rate; it does not evolve pitch angle, anisotropic pressure, magnetic
moment conservation or self-consistent MHD work exchange.

## Discrete method and ledger

Edges are increasing nonnegative joules. Populations are integrated cell
numbers per m3; energy is the arithmetic cell midpoint. At a face, b=-C*E_face
and F=b*N_upwind/width_upwind. Backward Euler gives a positive triangular solve
because b has one sign across the nonnegative energy grid:

    N_new,j - N_old,j = dt*(F_j-F_(j+1)).

There is no prescribed inflow at either energy boundary. Outgoing number is
reported separately at lower/upper edges, carrying that edge's energy. E=0
has zero drift flux. Grid outflow is a **numerical-domain inventory**, not a
physical orbit loss or thermal-ash source. A host must retain/account for it
or extend/reject its grid according to an explicit domain policy. No hidden
population removal, thermalization or volume rescaling is performed.

For face amounts A_j=dt*F_j and midpoint energies e_j, compute work independently
from the face sum

    W = (e_0-E_min) A_0 + sum_inner (e_j-e_(j-1)) A_j
        + (E_max-e_last) A_max.

This gives the discrete identities

    N_new + N_lower + N_upper = N_old,
    U_new + E_lower + E_upper = U_old + W.

Positive W heats particles; negative W cools them. This is external mechanical
work, not collisional heat to electrons/ions and not nuclear Q. Signed work and
spatial exchange still need representation in the source-state coupling layer.
The kernel does not yet change any BALDUR energy RHS or source-state context.

Candidate publication is atomic after every check. Finite trace values follow
IEEE rounding; measured state/flux conversion errors enter the extended-precision
balance checks. Invalid input, nonfinite results or ledger overflow clear all
valid-dimension outputs and return a status, without a density floor.

## Numerical accuracy and verified limits

The method is first order upwind and backward Euler. Moment conservation does
not prove spectrum accuracy. On a single [0,Emax] cell, cooling cannot move
population to a lower center and is unresolved; the explicit regression retains
this limitation rather than adding an undocumented ash sink.

For initial g=1 on[2,4] and zero elsewhere, frozen C=+/-1, t=.2 and grid[0,8],
the exact spectrum is g(E,t)=g(E/exp(-Ct),0)/exp(-Ct), with N=2 and
U=6exp(-Ct). Joint grid/time refinement from80 to640cells/steps gives:

| C | Relative energy error80 | Relative energy error640 | Spectrum L1 error640 / initial N |
| --- | ---: | ---: | ---: |
| -1 (heating) | 0.00327578 | 0.000408966 | 0.0645595 |
| +1 (cooling) | 0.00393550 | 0.000492435 | 0.0739570 |

The discontinuous spectrum remains visibly diffused at640cells. These are
convergence diagnostics, not a final high-accuracy kinetic grid choice. A
separate fixed-grid one-cell exchange checks time refinement against an exact
exponential solution. Tests also cover signed boundary work, zero drift,
nonuniform SI grids spanning11orders, energy-unit scale covariance, subnormal
rounding, independent cumulative budgets, input preservation, bad extents,
invalid coefficients and overflow.

GNU and Intel-r8 native/Fortran tests pass. Installed C11 and Fortran consumers
also pass. See validation/energy-work/. Production coupled trajectories,
operator-splitting convergence and source-state spatial/work accounting remain
required; no four-fuel window or full model acceptance is claimed here.
