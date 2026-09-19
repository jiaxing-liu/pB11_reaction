# Trial-owned source spectral packets

`fusion_c_coupled_sources_packets_trial` is an additive form of the covered
coupled trial. Its last outputs are a caller-owned mapped array `[2][5][7][cells]`
and `fusion_birth_packets_v1` metadata. Source0 is thermal, source1 fast-target;
channel and species orders are the existing nuclear IDs, with neutron6 explicit.
All weights are STEP AMOUNTS per source m3, not rates or extensive particles.
Multiply by the same explicit source volume used by the production source ledger.

The thermal mapped amount is the actual `grid/rate * thermal_burn.events` used
by the trial. Fast mapped amounts are captured in the actual target-network
reaction-loss loop as `loss/rate * spectrum`, before source merging. The packet
contains the original below/above N/E by source/channel/species and separate
thermal/fast event amounts. It does not call a second source solver, recompute
an approximate shape, or normalize to mixed species totals. Neutron packets are
included for accounting; they are not charged orbit markers.

Original below-grid physical energy is distinct from floor-mapped center energy.
The optional numerical floor continues to put below-grid particle number in the
first kinetic cell and keeps its separate ion-energy correction. Neither that
mapped amount nor the correction replaces physical below-grid values here.
External birth, subsequent FP evolution, heat deposition, escape and handoff are
not part of these source packets. Above-grid charged spill still follows existing
production rejection rules. Angular/pitch information is NOT reconstructed.

The new path checks packet N/E sums against the existing physical nuclear-birth
and neutron ledger. Channel/species stoichiometry and shape have independent
regression tests. Fine-grained packet representation can fail even if a coarser
legacy aggregate is representable: positive unrepresentable packet values are
not clipped. Only callers explicitly selecting the new entry point encounter
these additional diagnostic gates. Existing entry points retain original source
arithmetic and no O(cells) capture allocations. Opt-in workspace is approximately
70*cells long doubles plus70*cells doubles, plus small fixed metadata.

## Trial is NOT accepted

Status0 means only that this trial completed. The stateless library has no host
step identity and cannot approve an outer transport attempt. The caller must
bind the output to its source ticket, source volume, field epoch, time interval,
and model identity. On retry/rejection, discard it; publish only after the same
outer attempt commits. Never re-publish an older internally consistent packet
under a new ticket. This output is not a checkpoint/physics state owner and
contains no callbacks or persistent diagnostic side effects.

Output buffers clear at entry (mapped extent only when cell count is valid).
Candidates remain local until every physics and packet gate passes. Any failure
clears physical and diagnostic outputs. Nonoverlapping arrays are required.
Input state is never advanced. A direct ISO_C_BINDING interface and reversed-dimension packet type are in
fusion_coupled_sources_floor_fortran; Intel default-real and -r8 layout/call/error
tests pass. Host ticket/lifecycle integration remains following work. The direct
Fortran binding has the same caller-owned pointer/extent contract as C.

Regression validation covers mixed DT thermal/fast source separation, exact
legacy physical-state parity, per-source reaction stoichiometry, direct thermal
spectrum comparison, number/energy reconstruction, deterministic retries and
output clearing. The pB floor fixture separately checks original below-grid
energy and packet clearing on a late floor rejection. Consult coordinator
STATUS.md for live/terminal test evidence; source implementation alone does not
prove the complete host lifecycle or physical ensemble validation.
