# Exact intra-trial beam spectrum reuse (D086)

The coupled fast operator evaluates each edge's beam birth source before the
shared-target network solve, then evaluates the identical input again while
assembling products. The second call does not use depleted target temperature:
both calls intentionally use frozen fastTi. Reuse the complete first result
and full7*n double spectrum within that same invocation only. The existing
source-rate/debit comparisons, spill bound, network solve and second-pass
coefficient identity checks remain. No interpolation or spectral truncation.

Cache at most32MiB of dense spectrum payload per local trial. Metadata/vector
allocation overhead is additional; this is not a32MiB total-process limit.
Once the payload cap is reached, continue the original recomputation path.
If cache construction throws bad_alloc, stop adding entries and recompute;
previous entries remain usable. Cache does not survive return/rejection, has no
shared mutable state, and cannot reuse a different zone/background/timestep.
Public ABI and caller model controls are unchanged. Existing zero-rate entries
need no cache because their accepted loss is zero.

A two-scenario benchmark uses three populated projectile bins, with thermalDT
on/off, and serializes ALL thermal/kinetic arrays plus the full trial result.
The harness statically checks its ABI has no padding before byte comparison.
Uncached, normal cache,8KiB cap and injected cache-allocation failure produce
the identical SHA25686caad2acaf7ffe13719884842731a069ac8aa63bae63960c2742a52279bc238.
Source-call counts:12/6/10/12; measured wall seconds
.118510118/.062428993/.100599323/.117388394. The reduced-cap/failure variants
are temporary test source overrides, not public model settings. This single
fixture measurement is not a full-host performance claim.

GNU3/3 and Intel/r8 3/3 coupled-fast/thermal/table tests pass, and the full host
Intel build passes. Existing accounting, invalid input, retry and source-state
acceptance checks remain in test_fusion_coupled_fast.cpp. Root's parity harness
and source overrides are archived for reproducibility.

A separate source benchmark at actual zone32 and adjacent energies reports
.033..039s per complete DT source call. The400cell full grid costs22400bytes;
the exact actual point has58nonzero bins (29He4+29neutron), zero spill fields.
Sparse696byte estimate excludes metadata and is only a potential optimization;
this implementation deliberately retains the full dense spectrum. These are
standalone normalized source coefficients, not accepted host populations.

Actual host cache validation is LIVE55916: /tmp/dd-fast-cache-d086, case
bald_ifort-stdf90/run/exl-baseline/dd-fast-cache-d086,20us400cells,moving,
DDthermal+DTfast,tag86001. Compare all sevenCSV files with accepted DDsecondary07
only when complete. Existing D085 runs42952/30586 and longDT1244 stay frozen;
they use their original linked executables and must not be restarted for this
optimization. Full goal remains active; convergence/ash/loss are incomplete.
