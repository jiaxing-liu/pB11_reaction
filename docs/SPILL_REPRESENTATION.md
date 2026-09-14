# Fast birth spill and representable trial inventory (D082)

The coupled fast API previously rejected any positive out-of-grid coefficient,
even when no possible reaction in the supplied inventory could produce a
representable out-of-grid particle or energy amount. An evolved-D fixture had
old fast population2.9601203997565862e-291m^-3 at4538.1174509712828keV,
K=1.7427633952465272e-22m^3/s and above-grid alpha coefficient
1.9733775364326608e-316m^3/s. Even complete consumption bounds the outside
number by3.3518233846020156e-585m^-3. The public spill energy coefficient rounds
to zero; the implemented bound also includes its upper representable neighbor.

Before depletion, for each active channel/bin and charged species, bound the
possible outside number and energy by accepted fast inventory divided by the
lower representable neighbor of K, multiplied by upper representable neighbors
of both below/above coefficients. Summing competing channels is conservative:
no channel can consume more than the whole starting bin. Sum in long double
before testing, so individually unrepresentable contributions cannot evade an
aggregate check. Guard both integrated inventories and step-equivalent rates
using max(1,1/dt). A factor2 margin covers positive arithmetic rounding; permit
only bounds whose double conversion is exactly zero. Nonfinite/overflowing
bounds and every representable spill still return OUT_OF_RANGE.

There is no tunable tail cutoff, renormalization, heat/loss relabeling or change
to the independent birth spectrum. Normal sampled spectra and all nuclear
rate/debit gates remain unchanged. Upper neighbors protect moments that rounded
to zero at the source ABI. The second deterministic spectrum evaluation uses
the already checked first-pass bound; it does not accumulate the bound twice.
This implements the existing SI-double representability contract, not a
physical bound on all continuous tails or an energy-grid convergence result.

The physical regression accepts the tiny-inventory case with nonzero DT event
and reactant debits; identical source/grid with old fastD=1e15m^-3 still rejects
and clears outputs. GNU and Intel/r8 coupled-fast/thermal/table suites pass3/3
each. The original complete two-step Fortran controller fixture now passes:
commit, evolve the accepted distribution, bad final-zone rejection, retry and
second commit. Full host build passes. The actual DD zone39 source4 failure is
separate and remains under diagnosis; do not claim full DD host acceptance.
