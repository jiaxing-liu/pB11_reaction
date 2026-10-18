# D891: represent every thermal charged mapped bin

The former historical diagnostic selected one tensor item per bin and retained
that item's small original weight. It covered only about0.095% of mapped weight.
The new optional `--samples-per-bin 8 --seed 1729` uses the standalone library's
monoenergy categorical compression API on every validated bin tensor. Independent
caller-provided PCG64 uniforms select with replacement; each draw carries binN/m.
Duplicate draws are retained. The sum of N and N*K over every bin closes to
roundoff. This is a stochastic estimate of the full supplied **discrete tensor**,
not an assertion that7336 trajectories exactly integrate the continuous source.

Existing default subset selection and physical orbit settings are preserved.
Tails, beam-source and neutron accounts remain disjoint. The historical accepted
step1 fixture has917 eligible bins,7336 draws and no nonzero charged tails/beam;
it includes channels2,3,4 and is not a pB production result. Input packet and
field identity remain the previously verified accepted-bundle-d231 fixture.

The library API rejects disappeared positive CDF intervals and unrepresentable
weights/energies, uses staged publication and checks relative N/E closure within
64 double epsilons. Even a positive long-double interval can be unreachable by
a finite-precision uniform generator; this limit is explicitly documented. The
adapter records53-bit PCG64 uniforms and NumPy version. Samples do not preserve
angular moments exactly and never establish event-level nuclear correlations.

The analyzer groups compressed draws by zone/channel/species/cell, checks every
stratum and counts each bin only once. It reports separate event and
(event+unresolved) estimates, stratified empirical standard errors and a
simultaneous95% weighted-Hoeffding endpoint radius. Bounds assume independent
ideal uniforms; PCG64 approximates that model. They exclude finite RNG precision,
continuous quadrature/core/field/horizon/tail errors. Empirical zero standard
error does not rule out rare outcomes; the concentration radius remains positive.
These are model-LCFS classifications, not material-wall loss or escape rates.

Tests cover weighted CDF endpoints, zero/invalid/subnormal/overflow cases,
repeatability, N/E, exact known unequal-bin statistics, complete analyzer
accounting, unequal physical/reference volumes and legacy subset behavior.
Full library125/125 passes. ASAN/UBSAN passes with leak detection disabled;
initial LSAN ptrace/environment failure is retained, no leak-check claim made.

`run.py` launches the pinned actual7336-marker100ns diagnostic with unchanged
0.5ns initial step, level7 and1e6 per-level step budget,2GiB address cap and
48h wall budget. It refuses duplicate run-state and never automatically retries.
Source/binary/input hashes are captured in the live run-state.json. Session56447
was confirmed running; no completed orbit acceptance is claimed here. No active
D880 native binary/cache or source deck was changed.
