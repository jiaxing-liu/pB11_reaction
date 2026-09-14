# Explicit complete-source tables in coupled fast trials

`fusion_c_coupled_sources_trial` is additive to the existing coupled APIs. It
accepts optional five-channel thermal tables and a sparse array of borrowed
`fusion_beam_table_entry_v1` records. Each record identifies a channel, canonical
projectile slot, zero-based energy cell and immutable beam table. The Fortran
module/installed CMake target is `fusion_coupled_sources_fortran` /
`pb11::fusion_coupled_sources_fortran`.

Missing beam entries use the original full direct source. Present entries must
match every beam option, every grid edge, channel, slot and the cell-center energy
exactly. Duplicate records, null handles, disabled channels and the redundant DD
slot1 orientation reject, including when the population at that entry is zero.
A used table outside its temperature domain rejects; no fallback conceals a bad
entry. The caller owns all table handles and retains them throughout the call.
No constructor, mutable cross-step cache or physical state advance is hidden here.

One interpolated full spectrum supplies both target-network coefficients and
accepted births. The existing split remains: thermal burn, old-fast/remaining-
thermal burn, combined births, FP/escape and handoff. The old exact intra-trial
cache still reuses identical sources between coefficient and birth passes. All
source/Q/particle balances, spill representability checks and S/T shared survival
rules remain. Existing entry points retain their arithmetic when no beam tables
are supplied. Source construction/model choices remain external to the host.

The new usage output counts actual table and direct calls, including uncached
second passes, and reports used tables' sampled interpolation envelopes. These
are separate from result.max_source_* discrepancies, which retain their direct
quadrature meaning (sampled envelope for a table). Explicit table constructor
controls govern interpolation accuracy; no new default tolerance is introduced.
Sampled envelopes are not uniform error guarantees. All trial outputs, usage and
handoff diagnostics clear on any failure; nothing is published before success.

Validation includes three populated fast-D bins against thermal T. Empty entries
preserve direct outputs exactly; endpoint table/direct full S/T arrays, thermal
state and ledger agree exactly; a sparse two-table/one-direct mix also agrees.
Three independent temperatures check per-species full S/T shapes, grouped heat,
events, neutron energy and target debit. Thermal-plus-fast splitting is tested.
Every mapped model field, full output grid, duplicate/null/wrong energy/slot,
invalid selector and temperature domain failure is checked with clearing. C++
regressions pass in GNU Debug and GNU Release builds used with Intel Fortran.
Fortran tests exercise real table handles, layout, optional arrays and clearing;
GNU and Intel pass, as does an installed Intel consumer.

One representative 113-cell, three-populated-bin coupled benchmark took 1.31s to
construct three narrow-temperature tables. One hundred direct trials took2.89s;
one hundred table trials took0.108s, with identical endpoint physical states.
This includes local coupled FP work but is not a complete BALDUR speed benchmark.
Full-discharge table selection, construction budgets, temperature coverage and
convergence remain caller-level validation work. No assertion of ignition, ash
transport completion or nonthermal discharge convergence follows from this test.
