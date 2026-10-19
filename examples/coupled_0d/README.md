# Standalone coupled 0D C++ caller

Build against an **installed** library, independently of its source/build tree:

```sh
cmake -S examples/coupled_0d -B build-coupled -DCMAKE_PREFIX_PATH=/absolute/install/prefix
cmake --build build-coupled --target coupled_0d -j1
./build-coupled/coupled_0d
```

The only linked dependency is the exported public `pb11::pb11` target, found by
`find_package(pb11_reaction CONFIG REQUIRED)`. This project does not include
library source, internal or test headers, or BALDUR code. A C++17 compiler is
required. Any unsuccessful API call or failed numerical check exits nonzero.

Each fuel takes two 0.1 ms steps through `fusion_c_coupled_thermal_trial`:
thermal burn, charged birth, two-component implicit FP collisions, signed bath
heat and optional measured thermal-ash handoff. Neutrons remain a separate ledger,
never local heat. Carbon is an explicit inert bath, with a separate heat account.
The caller sets electron density from the total accepted positive charge of
thermal ions, both fast populations and explicit fully stripped carbon (mean
Z=6). Public nuclear-charge metadata supplies network-ion charges. The seed
alpha charge is included before initial Ue is formed at 5 keV. At each step
the caller refreshes ne from accepted inventories without resetting Ue; ne
is frozen inside that local trial. A relative 1e-12 total-charge check precedes
commit. This is a fully stripped tutorial assumption, not an ionization or
transport solver; mean Z is never inferred from a mixture's mean Z squared. External fast injection and escape rates are explicitly zero.

Each step first computes a candidate and discards it, checks accepted inventories
are untouched, repeats the identical trial, and compares all result fields. It
commits thermal numbers, both thermal energies and both kinetic arrays together
only after the API and ledger checks pass. This is caller-owned transaction
logic; the library has no hidden accepted time or cumulative state.

CSV columns report accepted SI number densities, electron/ion/fast energy
inventories, per-step reaction events, alpha handoff, signed thermal increments,
and residuals, including the frozen electron density and relative positive-
charge residual. `fusion_c_coupled_thermal_increment` obtains weak signed fluid
sources from the ledger without subtracting rounded final/old inventories.
Final Ue/Ui returned by the coupled trial are used for acceptance. Fast birth,
internal S-to-T transfer and fluid ash are different inventories: alpha handoff
may legitimately remain zero during this short hot-alpha demonstration.

Grid, constants, explicit source-model options and numerical gates reuse
`tests/test_fusion_coupled_thermal.cpp`; the pB 100 keV ion temperature, 2.5 MeV relative-energy cutoff and
1e-18 keV first grid edge follow
`tools/study_coupled_thermal.cpp`. Both cases seed 1e12 m^-3 hot alphas to expose
collision heat immediately. The pB relative cutoff obeys the public 12 MeV parent-domain limit (the
DT fixture uses 5 MeV). The 257-cell grid spans zero to 25 MeV; birth
relative/CM quadrature orders and angular orders are eight, with explicit 0.002
source discrepancy gates. Handoff uses the existing measured 0.001 shape/mean
energy gates. These numbers are tutorial settings, not established time/grid/
angular convergence or differential pB data accuracy. The source model uses the
explicit entrance-proxy remainder policy and bare nuclear rates. This example
is neither a device simulation nor whole-M3/M7 completion.

Checks preserve the existing coupled integration fixture's fast-number, fast-
energy, source-Q and total-energy tolerances. API success also enforces its
stored-inventory roundoff-aware residual contracts. There is no fixed empirical
heat fraction, discarded neutron heat, tolerance relaxation or failure-to-zero
fallback. This small example does not exercise independent checkpoint/restart;
existing M4 qualification is separate.
