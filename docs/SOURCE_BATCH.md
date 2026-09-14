# Synchronous source-state batch acceptance

`fusion_c_source_state_commit_many(count, states, tickets)` is an additive C
API. It publishes a set of already validated per-zone trials together. This
is failure atomicity under exclusive caller ownership, not simultaneous
hardware writes, cross-thread synchronization or crash-safe persistence.

The caller supplies 1..1000000 distinct nonnull contexts and their current
staged tickets. All accepted times, trial times and accepted epochs must be
exactly equal. Grids and model tags may differ. Equal-time contexts at different
epochs reject: callers must not accidentally combine old and new zone states.

Validation includes duplicate detection using a separately allocated/sorted
pointer vector. All allocation and checks precede publication. Publication
uses only no-throw vector swaps and trivial field copies, shared with the
existing single-context commit. A failed batch changes neither accepted nor
pending state in any context. The caller may repair a bad member or ticket
and retry the group; a failed batch does not discard other valid candidates.

The operation transfers populations, cumulative base and inert ledgers,
format promotion, time and epoch. It does not advance physics, evaluate sources,
remap zones, enforce thermal fuel availability or commit host-owned thermal
and geometry arrays. The caller must finish all fallible host validation and
allocation before batch publication, then publish its own same-epoch state
without a further failing operation. This API alone is not full host atomicity.

Tests use different cell counts/model tags and mixed legacy/extended contexts.
They verify stale last tickets, duplicate/null members, absent staging,
mismatched target times and mismatched accepted epochs leave all accepted
snapshots unchanged. A corrected group then commits, repeat commit rejects,
and a good candidate remains usable after a failed group. Existing single-zone,
inert accounting, restart and physical coupled tests remain the regression
reference. See docs/validation/source-batch for recorded compiler checks.

## Validation evidence

GNU CTest passed 60/60 after the C++ change, before adding the Fortran batch
wrapper. After the final Fortran edits, its GNU binding test passed separately.
Intel ifx with default REAL=8 passed all four batch/state/inert/Fortran binding
tests on the final tree. The new Fortran test checks unequal array extents,
last-ticket failure, preserved nonzero inventories and exact cumulative ledgers
after successful group publication. Logs and final source hashes are archived
in `validation/source-batch/`.

The archived strict C11 consumer compiled with `gcc -std=c11 -Wall -Wextra
-Werror`, linked with g++ against the GNU static library, and returned zero.
It creates, stages and commits two contexts and checks the accepted time/epoch.
This was a source-header/static-library check, not an installed-package test.
The Intel build initially reported an unknown new target from an old generated
Makefile; reconfiguring the existing build directory resolved this before tests.
