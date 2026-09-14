# Birth-table validation studies

This directory archives the completed numerical acceleration checks.
CSV and log artifacts are copied byte-for-byte from `/tmp`; their SHA-256
digests are in `summary.json` and `SHA256SUMS.txt`.

## Inputs and scope

- The birth-table microbenchmarks use four off-construction samples for DT
  and pB. The declared acceptance gate is 0.1% for each rate/debit/number/energy
  metric.
- The coupled comparisons retain all 64 accepted time points for each fuel and
  compare the table run with its direct-evaluation reference.
- `electron_heat_J_m3`, `network_ion_heat_J_m3`, and `inert_heat_J_m3` are
  reported independently; the last is the carbon-bath heat ledger.
- This is a numerical acceleration and parity study. It is not an EXL study
  or a fair performance ranking of fuels.

## Build provenance

- Birth-table study binary: `/tmp/study_birth_table`
  (SHA-256 `e9915054854e5d5b73acec8493c5d2c38ed84b3bc0f7961412cb749146fb7951`).
- Coupled driver binary: `/tmp/study_coupled_thermal_table`
  (SHA-256 `3a8a703cd307d12a8ef128c15c469053d4bb85d1c0beb980f4a6e8b6f0b6a1e2`).
- The microbenchmark binary predates the two sampled-direct diagnostic fields;
  the successful interpolation values are unchanged. No source hash is used
  as compile provenance.

## Results

- **DT**: off-construction max `0.00057909`, direct/lookup timing ratio `1660`–`1849`, table construction `0.46479 s`.
  Terminal independent heat relative errors: electron 0.00253415%, network_ion 0.00390554%, carbon 0.00395542%.
- **pB**: off-construction max `0.000137032`, direct/lookup timing ratio `7.253e+05`–`7.711e+05`, table construction `40.4608 s`.
  Terminal independent heat relative errors: electron 0.0306969%, network_ion 0.0293819%, carbon 0.0293487%.

The four-panel figure is `birth-table-validation.png` (also available as
`birth-table-validation.pdf`). Panel 2 uses a logarithmic timing axis and
annotates construction cost separately; panels 3 and 4 show parity errors
against the direct coupled runs.
