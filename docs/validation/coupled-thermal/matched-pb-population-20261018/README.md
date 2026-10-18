# D881 accepted-state pB population-model coverage diagnostic

This fixture samples all50 native zones at accepted steps0,1280,6400,12801
(0,1.9996875,9.9996875,20.00125ms) from the completed matched pB case D878.
The original compressed profile SHA256 is pinned in input-receipt.json.
These are accepted-state diagnostics, not exact pre-trial source-entry replay.
No sparse-snapshot integration is used to claim the total reaction yield.

`diagnose_pb_thermal_population.cpp` is the source of the independent public-API
consumer now provided in the library examples directory. It reads whitespace
separated label,kTa_J,kTb_J triples from stdin and returns all three effective
population rates plus the overlapping extrapolated diagnostic subset. Fractions
remain empty when total reactivity is zero; no denominator floor is introduced.
The four command-line options are continuation,pb_low,narrow_peak,continuum_scale.
The fixture uses1,0,0.051,1, matching the active host population parameters.

Rebuild against an ordinary current library build, then reproduce raw rates:

```sh
g++ -std=c++17 -O2 -I/path/to/pB11_reaction/include \
  /path/to/pB11_reaction/examples/diagnose_pb_thermal_population.cpp \
  /path/to/build/libpb11.a -o /tmp/pb-population-diagnostic
/tmp/pb-population-diagnostic 1 0 .051 1 < states.txt > populations-new.csv
MPLCONFIGDIR=/tmp/pb-population-mpl python3 analyze.py
```

`analyze.py` uses archived profiles.csv/populations.csv; compare regenerated
rates before replacing a fixture. It checks all200 labels, temperature mapping,
finite rates, population sums and subset bounds, then writes summary and figures.
The run took5.60s and does not mutate host cases or caches. The zero-rate,
noninteger-policy,negative-temperature and empty-input frontend checks passed.
Binary/static-library SHA256 and source identity are in tool-checks.json.

Native diagnostic Ti uses the historical1.6021e-16J/keV temperature conversion;
this does not change native energy diagnostics or assert exact trial-entry Ti.
Independent equal-temperature Maxwellian event rates are weighted by native
nH*nB*zone-volume at each snapshot. Source-weighted continuum-table extrapolated
fractions range99.5165–99.5597%; effective alpha0 peak fractions are about1.590–
1.592%. This is a coverage warning within the explicit endpoint-held Taskaev
proxy, not a99.5% numerical error or measured physical branching ratio. The
remainders are not declared pure alpha1. The extrapolated subset overlaps the
three populations and must not be added again.

At these cold states the absolute pB source is tiny. The diagnostic does not
replace source-spectrum/model sensitivity, physical-loss closure, actual source
entry validation or the fullM3 gate. No core library source/ABI/cache key or
running binary is changed. See docs/PB_POPULATION.md for the existing model
contract and data references.
