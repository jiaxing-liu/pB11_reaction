# Thermal and fast-target reaction / kinetic trial

D078 extends the existing coupled operator through additive C/Fortran entry
points. It uses the validated beam-birth and shared-target depletion operators;
there is no new reaction spectrum, empirical heating fraction or accepted clock.
Existing thermal APIs and result layouts are unchanged.

## Operator sequence and ledger

1. Existing simultaneous thermal burn freezes its coefficients at old common Ti.
2. Mix remaining thermal ion energy, including explicit inert ions, to obtain
   the background Ti for fast-target coefficients.
3. Build reaction edges for every nonempty old canonical fast species/energy bin
   against each enabled remaining thermal target. Evaluate full beam sources,
   collecting K and the reaction-conditioned target energy moment from the SAME
   quadrature. Shared model/event controls come from thermal_options.birth;
   new fast options provide independent channel selectors and angular controls.
4. Solve all these fast/target edges simultaneously. For identical DD reactants
   use one fast/thermal orientation per nuclear branch, with no half factor;
   DT and other unlike reactants can have both projectile orientations.
5. Reconstruct S and T with the same per-bin nuclear survival denominator.
   Check their combined returned number against the network result. Keep per-edge
   event amounts to preserve weak consumption; never use old-minus-new as the
   only source diagnostic.
6. Reevaluate each nonzero event's full spectrum and multiply by events/K.
   Require repeated K/M to match the first evaluation. Sum thermal and fast
   nuclear births into S, then apply the EXISTING FP, escape, bath heat and
   measured handoff operations using surviving old S/T inventories.
7. Check per-species kinetic and total-local number/energy accounts including
   fast_consumed, thermal_consumed, full neutron energy and canonical Q. Only
   then publish all trial outputs. Inputs/accepted source state remain unchanged.

Fast kinetic energy removed by a reaction belongs to the product energy budget,
not direct fluid heat. Target energy likewise uses its selected M/K, not 3Ti/2.
New nuclear products and external injection enter after this step's fast burn;
within-step secondary reburning is not inserted. This is a first-order split,
not simultaneous implicit bath/reaction iteration. Fast-fast reactions are not
included. Time/order and distribution/grid convergence remain required.

The full existing source ledger receives both thermal and fast consumption.
External birth remains externally supplied birth; new nuclear products are not
relabeled external particles. The existing fluid-increment API therefore yields
target depletion and signed heat exactly once. Source-state stage/discard/commit
is reused; BALDUR must still commit fluid/geometry and kinetic candidates together.
Caller model tags must include fast selectors and quadrature/model options.

## Source storage and boundaries

Only O(cells + reaction edges) auxiliary storage is retained; no source spectrum
is stored for every input bin. Two deterministic direct evaluations per active
edge trade CPU for bounded memory. The thermal-table entry points accelerate
thermal birth only; fast source evaluation remains direct and can be expensive
for dense pB distributions. Prepared fast coefficients/tables and host integration
are still required for practical long runs; success is not a performance claim.

Existing source rate/debit accuracy gates apply to both source types. Charged
below/above-grid spill rejects the complete trial. Neutron spills are carried
explicitly in the neutron ledger. No particle is dropped, clipped, treated as
escaped, or thermalized solely because it misses the numerical grid.

The initial pB fixture incorrectly inherited the DT5MeV relative cutoff, beyond
the pB event domain, and returned OUT_OF_RANGE. With a supported2.5MeV cutoff,
its1eV first center still rejected below-grid alpha coefficient1.009662e-33m3/s;
a1e-6eV center still had8.559463e-80m3/s below-grid. The final257-cell fixture
extends the first center to1e-12eV and passes without changing any physical gate.
Only its first edge changed; the higher grid remains the declared hybrid grid.
This establishes coverage at the finite quadrature used, NOT a universal bound
on the continuous low-energy tail or convergence of a grid near zero energy.
All failed inputs and diagnostic logs are retained.

## Validation scope

Independent DT tests have nonzero old S and T fast deuterium, thermal tritium,
explicit neutron/alpha counts, nuclear Q plus both reactant debits, separate
external birth and collision heat, invalid-option clearing, deterministic retry,
and late charged-spill failure. Existing source state remains unchanged through
stage/discard, accepts one complete trial, and rejects repeated commit without
duplicating epoch/ledger. This validates that library owner, not full BALDUR
fluid/geometry rollback. Two frozen old identity/DT fixtures retain byte-identical
full results and population arrays against the pre-change library.

Root physical fixtures:
- Fast DT:1.9587534936013994e12 events/m3 in1e-4s; fastD,thermalT,neutrons equal.
- ThermalDT + fastDD(two branches) + fastDT(both projectile orientations):
  DTevents3.0679215747624816e13; total neutron number3.0689603982567863e13.
  Direct and exact-temperature-knot thermal-table output summaries match.
- Fast pB:3.7937899606692719e9 events/m3; alpha/event ratio3 and product energy
  including both reactant debits close at1.18e-15/1.23e-15 relative. No neutron.
  This uses one populated proton bin and30keV thermal target, no thermal burn.

The representative probes disable handoff; automatic ash transport is not
established by their alpha birth count. The full user PLAN, host activation,
neutral retry/restart, practical source evaluation and four-fuel windows remain
active work. Refer to docs/validation/coupled-fast for exact source/build evidence.

Final GNU4/4 and Intel/r8 4/4 selected suites pass, including a nonzero-fast DT
Fortran call, previous thermal/table regressions, and source-state rollback.
Full Intel BALDUR build passes. Standalone installed Fortran runs the same
nonzero-fast fixture through the exported target; strict C11 validates the
installed header, symbol and null-output contract. These do not activate the
new source owner in a BALDUR production input.


## Diagnosed entry points (D079)

The three existing fast trial names also have additive `_diagnosed` C entry
points with a required final `fusion_handoff_diagnostics_v1 *` argument.
The Fortran counterparts place diagnostics before status. They expose the same
candidate S-to-T transfer/thermal handoff observations as the existing thermal
operator; they do not add state, change physics, or alter acceptance ownership.
Null diagnostics is an error; diagnostics and trial outputs clear on failure,
including Fortran extent validation. Existing non-diagnosed interfaces remain.

The C++ test exercises all three diagnosed routes with nonzero fast DT,
compares complete output arrays/results to their non-diagnosed counterpart,
and checks candidate/transfer observations plus null/invalid clearing.
GNU and Intel/r8 selected C++/Fortran suites pass 2/2 each. BALDUR's D079
controller fixture additionally exercises diagnosed standard/effective-charge
Fortran dispatch and atomic zone rejection/retry. Host evidence is in
`docs/exl50u-program/controller-fast/`; this is not production driver activation
or a dense multi-step performance certificate.
