# Explicit pB source budgets in complete coupled trials

Include `fusion_coupled_sources_budget.h` and call
`fusion_c_coupled_sources_packets_budget_trial`. It retains all arguments and
ownership rules of the covered/packet trial and appends mandatory
`fusion_thermal_birth_budget_policy_v1 *pb_policy` and
`fusion_coupled_source_accuracy_v1 *accuracy` arguments. Inputs and outputs must
not overlap. No physical state is advanced: accept or discard fluid, kinetic,
packets and diagnostics together in the caller.

The explicit policy affects only the thermal pB source. Other thermal reactions
and all fast-target sources use their original full model. A null thermal-table
array selects direct thermal sources; otherwise an active pB source requires
an exactly matching budget table created with the D915 explicit constructor.
A full pB table, mismatched source/grid/policy, or out-of-domain thermal
temperature rejects. Old coupled entries continue rejecting budget tables,
even zero-budget tables. No policy is inferred from a handle or environment.
The policy is validated even when pB is disabled or lacks a reactant.

The separate accuracy output preserves existing options/result ABI:

| mode | Meaning | Populated nested member |
| --- | --- | --- |
| 0 | No pB thermal source evaluated | Neither |
| 1 | Direct budget source actually evaluated | `direct` |
| 2 | Immutable budget table used | `table` |

On successful whole trials `version=1`, including inactive mode0. Both estimates
are zero for inactive sources. Direct estimates are the ideal positive-mixture
number/energy bounds plus their declared engineering roundoff allowances,
rounded upward on conversion; zero-weight sources have zero estimates. Table
estimates are `table.base.max_validated_number_L1/energy_L1`: sampled composite
interpolation plus direct-source bounds/allowances. Table reports must not be
read as actual direct quadratures at the trial temperature. Only the populated
nested member is meaningful; the other stays zero.

These are numerical estimates against the same finite source reference,
including cells and spills, not nuclear/continuous-source uncertainties,
off-sample interpolation guarantees, or subsequent plasma response bounds.
The source rate/debit discrepancy gates remain separately enforced. Existing
charged-spill, numerical floor, beam coverage, FP, handoff, packet
representability, weak-source and conservation checks remain enforced.

All accuracy and physical outputs clear on failure for valid dimensions.
Provisional observations are held locally and published only after the entire
trial passes. No extra source evaluation, hidden cache, mutation, implicit
retry or process-global policy is introduced. Tables remain immutable and can
serve independently owned concurrent trials.

Fortran users link `pb11::fusion_coupled_sources_budget_fortran` and use
`fusion_coupled_sources_packets_budget_trial`. Its mapped-packet shape is
`(cells,7,5,2)`, including neutron; thermal/fast inventory arrays remain
`(cells,6)`. Malformed shapes clear outputs and reject before forming C pointers.
The module re-exports required common ABI types and policy/report types.

The engineering tests use tiny finite quadratures to verify interface, exact
zero-budget parity, matching, serialization reuse, failure atomicity and language
layout. They do not qualify physical source convergence or a BALDUR trajectory.
Explicit host configuration/owner selection, restart provenance and numerical
response convergence are separate integration steps.
