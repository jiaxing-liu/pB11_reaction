# Separate numerical energy account

The source-state owner now accepts six signed numerical transfers to the ion
reservoir, one per fast species, in J/m3. These are separate from physical heat,
nuclear Q, inert-bath heat and transport/work. A negative transfer means extra
energy enters the kinetic grid. The per-species energy balance includes the
signed transfer; nuclear stoichiometry/Q checks remain physical and unchanged.

Use stage_numerical or stage_volume_numerical with explicit inert and numerical
arrays. The moving-volume call scales the numerical density by
source_volume/reference_volume, exactly as source/inert accounts. It rejects
nonfinite, overflowed or nonzero-to-zero scaled amounts. Accepted cumulative
accounts use reference-volume normalization. Discard/rejected stages cannot
promote the context or change accepted accounts. Commit/commit_many publishes
populations and all staged accounts together.

A successful commit through either new stage promotes the context even when
all numerical amounts are zero. Promoted contexts reject legacy stage/snapshot
calls; callers must explicitly preserve this account. New snapshots also read
legacy contexts, returning zero numerical transfers. Fixed-volume snapshot
failure leaves kinetic arrays untouched; volume snapshot failure clears them,
matching the respective existing API conventions.

Restart versions4(stationary) and5(moving) insert six signed finite binary64
words after the existing inert-heat words. They retain the existing checksum,
tag/time/epoch checks and full cumulative inventory validation. Promotion needs
an accepted step; numeric formats with epoch0 reject. Legacy contexts keep
versions1/2/3 byte layouts. Old readers reject the new version numbers.

The numerical state API does not calculate the correction or move host fluid
energy. Callers pass the coupled floor ledger's correction exactly once and
must commit/restore their ion/electron/geometry state at the same epoch. The direct fusion_source_state_numerical_fortran binding and
fusion_c_coupled_numerical_increment helper now expose these accounts. The
BALDUR stage adapter can explicitly accept them; controller/driver propagation
and full-host integration remain pending. This API is not full host rollback or physical convergence proof.

The numerical increment helper adds all six corrections to raw physical terms
in long double before converting the final ion increment to double. It preserves
zero-correction legacy outputs and rejects missing/nonfinite/overflowing inputs
with cleared output. It does not divide by time/volume or advance host state.
Never include the same numerical amount in physical heat as well.
