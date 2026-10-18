# Immutable pB birth tables with an explicit source budget

`fusion_birth_table_budget.h` extends the independent library's immutable
common-temperature pB source tables. `fusion_c_pb_birth_table_create_budget`
takes the existing physical source options and interpolation controls plus an
explicit `fusion_thermal_birth_budget_policy_v1`. The caller owns the returned
`fusion_birth_table_v1` handle. Generic evaluate, info, destroy, pack/unpack and
content-digest functions operate on it. No table construction or evaluation
advances plasma, fuel, fast particles, cumulative diagnostics or accepted time.

The new Fortran module/installed target is
`fusion_birth_table_budget_fortran` / `pb11::fusion_birth_table_budget_fortran`.
It reuses existing source/control/coefficients types and generic table lifetime
and byte wrappers. New entry points expose budget construction, complete budget
information, exact request matching and each stored knot's source-budget report.
For Intel Fortran consumers select `LINKER_LANGUAGE Fortran` on the executable.

## Sampled errors and their denominators

Every expensive direct sample uses the explicit pB source-budget API. It keeps
the full relative-energy/CM reaction measure, with a conservative positive
mixture of angular orders. For pB the only product is alpha, so the expected
number is3*K and expected energy is Q*K+both reactant energy moments. The
sampled table/direct number and energy L1 distances include every grid cell
and all species' number/energy spills and use those expected source totals as
denominators. A zero denominator requires zero difference. No absolute source
floor or discarded low-rate region is introduced.

At each sampled quarter/mid/three-quarter temperature the composite estimate is

```
L1(interpolated budgeted table, finite full-order source)
 <= L1(interpolated budgeted table, direct budgeted source)
    + direct source's positive-mixture bound
    + separately declared source engineering roundoff allowance.
```

Endpoint source bounds must also meet the table's existing number/energy
controls. A table fails construction if either its direct-source or sampled
composite gates fail within the explicit work/depth limits. No partial handle
escapes. One convex log-temperature interpolation weight remains shared by
rate, both reactant debits, every grid value and every spill, preserving
positivity and linear particle/energy identities.

The base info's validated number/energy maxima now describe these sampled
composite estimates for budget tables. Additional info preserves validated
interpolation-only maxima, maximum effective source bounds (including the
engineering allowance for nonzero sources), maximum sampled allowances, the
complete source policy, and all existing controls and rate/debit diagnostics.
Each stored knot retains its complete parent counts, orders, ideal bounds,
allowances and source weights. Zero sources have zero effective source bounds.
Computed finite budget-branch error ratios and composite values are converted
upward to public double. Earlier arithmetic is not directed interval arithmetic.

The source allowance remains an engineering allowance, not a certified bound
on every floating-point operation. The table checks are sampled, not a
supremum error guarantee. Off-sample validation and source/energy-grid/time
convergence remain necessary. Neither angular/reference convergence nor nuclear
population/model uncertainty follows from a successful table constructor.

## Explicit identity and portable bytes

Budget-enabled tables use wire format2. The portable little-endian payload
adds the complete source-accuracy policy, sampled error metadata and per-knot
reports. Generic unpack checks exact byte length, current generated kernel
identity, checksum, numerical domains, particle/energy balance and internal
policy/report consistency, then publishes the handle. Import does not repeat
source quadrature or certify unsampled accuracy or data authenticity. Content
digests cover the policy, reports and complete canonical bytes.

Full-source constructors and format1 retain their existing behavior. Old public
request matching and the old coupled evaluator's internal model matching reject
all budget-enabled handles, even when the requested source budget is zero.
Call the explicit budget matcher to compare temperature limits, source options,
policy, controls and full grid exactly. Changing any policy field prevents a
match. Unknown versions and inconsistent zero-budget/low-order reports reject.
This prevents a table from acquiring an implicit accuracy policy through reuse.

Null/shape/error clearing and buffer-ownership rules follow the existing table
interfaces. New info, matches and knot-report outputs clear on failure;
knot indices are zero-based. Immutable evaluate/info calls share no mutable
library globals and permit concurrent callers.

This increment does not activate budget tables in the coupled source API or
BALDUR configuration. That requires explicit policy propagation, cache-owner
matching and restart metadata through the host's thin adapter. The generated
source identity changes honestly; old cache files must not be relabeled.
Current frozen host/library jobs keep their original paired binaries and data.
