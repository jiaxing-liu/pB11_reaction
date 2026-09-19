# Regular-coordinate inverse mapping, revision 1

`fusion_c_flux_inverse` in `fusion_flux_inverse.h` maps physical R/Z metres to
u=xi cos(theta), v=xi sin(theta), using the explicit regular-core snapshot.
It does not change the field-v1 or core APIs. Callers own arrays, model choice,
position tolerance, iteration budget and optional initial seed; no host state,
allocation or I/O is introduced. See the installed header for the full contract.

Newton works in u/v, with a native axis-linear or caller seed and up to 24
halved line-search trials. Seeds/trials are projected onto the exported disk.
The exact physical axis is handled without a finite exclusion region. Success
requires the supplied physical residual tolerance. Returned theta uses atan2;
compare u/v near the axis rather than an undefined axis angle. The iteration
budget limits residual/Jacobian evaluations; an update on the last iteration
needs a subsequent evaluation to be recognized as success.

NOT_CONVERGED (1002) means stalled/exhausted inversion, including possibly
external points. It never means physical escape. Invalid snapshots retain their
validation/numerical error codes. Every error clears the output. No reliable
inside/outside classifier, global uniqueness guarantee or orbit loss closure is
provided. Boundary classification remains a separate API requirement.

Validation: shifted-ellipse axis/core/outer-boundary roundtrips, two core choices,
explicit and automatic seeds, low-budget and outside nonconvergence, degenerate
axis rejection and output clearing; three field/core/inverse CTest targets and
ASan/UBSan. The coordinating BALDUR D215 report archives actual D209 snapshot
multi-seed, matching-surface and deterministic replay checks. These are finite
sample evidence, not continuum topology or physical confinement validation.
