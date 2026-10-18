# D885: prepared transfer constants, reviewed isolated prototype

The internal helper computes fixed projectile/bath long-double logarithms once
per bath and keeps the original left-associated normalization,y2 formula,
exp(normalization-y2),conversion/tail policy and bath-wise accumulation. Public
scalar validation order and output clearing remain unchanged; coupled consumes
only validated prepared inputs. There is no public C/Fortran ABI change.

The frozen old source reference and74184 input comparisons have zero status/
double-bit differences.67993 valid inputs also compare the prepared evaluator,
including subnormal tails,zero outputs,numerical failures and invalid inputs.
Six affected existing tests pass. Six pairs of508-grid transfer timing show
2.2027x path speedup; this is not a host speedup claim.

The more relevant frozen complete coupled FP microbenchmark uses the actual
D874 source/508 edges and1000trials each. Six alternating pairs give median
6.41965s old versus5.63441s prepared,paired speedup1.13763x (~12.1% saved).
All50040 output bytes are identical, SHA256 matches the original D874 run.
NoLTO,mesh change,source cutoff,fast-tail deletion or physical tolerance change.
Bath values are representative fixed input, not captured runtime baths.

Root reviewed the formulas/call invariants and applied the proposed code plus
independently diagnosed D886 geometry fix to the external working tree for
D887 full GNU-C++/Intel-ifx regression. This changes future source identity;
current D880 native/CLI binaries remain frozen at54b and are not rebuilt.
Only completed full regression can qualify the next library revision. Full
host speedup and new-identity cache/trajectory acceptance remain unproven.
