# Standalone pB temporal refinement — terminal archive

All three resolutions (64/128/256 steps) are complete. This folder is the
independent-library copy of BALDUR checkpoint D869, paired with BALDUR
commit07ab785 on feature/exl50u-fusion-library. There is no live computation
to restart from this archive. Original paths in manifest/run.py preserve
execution provenance; run.py is the original host-side launcher, not a
relocatable library benchmark launcher.

To regenerate figures from this self-contained archive:

```sh
MPLCONFIGDIR=/tmp/pb-refinement-mpl python3 plot.py manifest.json /tmp/pb-refinement-new-figures
```

The output directory must not already exist. The plotter checks all CSV hashes,
completion, finite/increasing times,1ms endpoints and energy closure. Source
and binary provenance are pinned in manifest.json. Run conditions and the
bounded interpretation are documented in../../../COUPLED_THERMAL.md.

The128->256 carbon-heat difference is0.124317J/m³ (0.52079%); electron heat
and network-ion heat differences are0.15278% and0.48282%. The maximum
accumulated energy residual is1.005e-15. These are local temporal differences
at fixed energy/angular grids, not a wholeM3 or EXL-50U four-fuel validation.
