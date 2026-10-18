# Conditional cold pB thermal spectrum diagnostic

The standalone example diagnose_pb_thermal_spectrum.cpp uses only public C APIs.
Compile it against an ordinary installed/build library with C++17. It reads
cells, cells+1 SI energy edges, all17 birth options in declaration order, then
label,kT[J] records from stdin. Run with two new output filenames; optional
`base-only` limits the scan to the baseline remainder0/continuum-scale1.
The ordinary scan compares baseline, all-low, all-broad, continuum-scale0/2,
each at input nq/ncos and doubled orders capped at1024. It preserves all other
supplied relative/CM/model controls. CSV outputs contain scalar marginals,
all original spill amounts and rate/debit/energy diagnostics. No density factor,
modified total cross-section, full event correlations or host source change.

The actual four accepted-state maximum-rate zones (Ti0.099,0.294,0.415,0.478keV)
use the BALDUR508-cell grid. Forty evaluations completed with N/E closure below
4e-16 and exactly identical rates across variants/orders. BUT normalized cell
L1 changes8->16 are about1.26--1.31: conservation does not establish shape
accuracy, and these curves are explicitly UNCONVERGED. They cannot qualify a
physical model sensitivity or first-version complete spectrum. A separate
one-zone base-only32->64 calculation is in progress, not yet accepted.

The approximately99.5% continuum-table extrapolation proxy is model coverage,
not a99.5% spectrum error bar. Its population choices really enter the
thermal quadrature node by node; the four-state diagnostic is not exact trial
entry or full-trajectory source integration. Remainder alternatives are not
measured branch decompositions or confidence bounds.

Prototype CLI/known nonzero synthetic checks are archived. Final source adds
optional base-only control with unchanged default scan and passes strict build.
The library core/header is unchanged; kernel identity remainsa5b4a289..., and
running caches and matched native cases are untouched. The actual raw CSV,
provenance and analysis live in BALDUR checkpoint-cold-pb-spectrum-20261018.
