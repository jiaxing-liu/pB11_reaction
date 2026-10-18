# D869: pB local temporal refinement

Live workspace data: `../cache/benchmarks/pb-time-refinement-d869/` relative to
the parent of the BALDUR repository. Inspect the live manifest and the original
execution session before considering a restart. `launch-d869.json` is only a
launch snapshot. Each resolution has a 12-hour wall budget.

`completed-64.json` and `steps-64.csv` are terminal current-source evidence:
64 steps reach 1 ms in 916.69 s; maximum accumulated energy residual is
8.31e-16. All 28 numeric columns match the archived 64-step reference exactly.
Electron heat is 123768.168688 J/m³, network-ion heat 4113.227200 J/m³ and
carbon heat 24.243451 J/m³. This is a reproduction check, not the unresolved
128/256-step convergence result. Those resolutions run sequentially.

The128-step run has now terminated:1843.34s, maximum energy residual1.01e-15.
64→128 carbon-heat difference is0.24855J/m³ (1.0358%), network-ion heat0.9591%,
electron heat0.3033% and fast stored energy0.2596%. See completed-128.json and
steps-128.csv. The256-step resolution is now terminal in3736.68s. The128→256 differences
are fast energy0.13035%, electron heat0.15278%, network-ion heat0.48282%
and carbon heat0.52079% (0.124317J/m³). The three-resolution maximum
accumulated energy residual is1.005e-15. These decreasing differences support
approximately first-order temporal convergence at this point; they do not
close the fullM3 gate. Final manifest,256-step CSV, resources and PNG/PDF
figures are archived here.

Reproduce the completed three-resolution figures:

```sh
MPLCONFIGDIR=/tmp/pb-refinement-mpl python3 \
  docs/exl50u-program/checkpoint-pb-time-refinement-20261018/plot.py \
  /home/cloud/cloud/research/pB-baldur/cache/benchmarks/pb-time-refinement-d869/manifest.json \
  /tmp/pb-time-refinement-d869-figures
```

The plotter verifies the completed sequence, CSV hashes, finite histories,
increasing physical times, endpoints and energy closure before writing figures.
It produces PNG/PDF histories, separate electron/network-ion/carbon heat panels,
endpoint-difference bars and a numerical summary with absolute differences.
It rejects the actual running manifest without creating figures. A rendering
check using the existing archived 32/64-step CSVs reproduced their 2.04719%
carbon-heat endpoint difference; this does not replace the new finer runs.

Conditions are a local stress study: ne=1e22 m⁻³, initial Te=5 keV and Ti=100 keV,
carbon=0.001ne, 800 energy cells and nq=ncos=8, fixed 1 ms physical duration.
They are not EXL-50U operating conditions or a four-fuel material comparison.
Temporal differences do not bound grid, angular, spectrum or host errors, and
the M3 completion gate remains open. No time-zero samples are invented.
