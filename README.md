# pB11_reaction

A standalone C++17 library for proton–boron-11 fusion, with an explicit
thermal DT/DD/D-He3 network extension. The pB astrophysical factor, fusion
cross section and Maxwellian reactivity are described by:

> A. Tentori and F. Belloni, “Revisiting p-11B fusion cross section and
> reactivity, and their analytic approximations”, Nuclear Fusion 63 (2023)
> 086001. [doi:10.1088/1741-4326/acda4b](https://doi.org/10.1088/1741-4326/acda4b)

The public API uses `double` precision. Implemented modules include:

- [Thermal network, cross sections and particle/nuclear ledgers](docs/THERMAL_NETWORK.md).
- [Classical Coulomb coefficients and conservative energy kinetics](docs/KINETIC_PRIMITIVES.md).
- [Finite-temperature two-component transfer and controlled thermal-handoff studies](docs/THERMAL_TRANSFER_DESIGN.md); internal kinetic transfer is separate from fluid ash.
- [Thermal nuclear-window contribution study](docs/NUCLEAR_WINDOW_STUDY.md), separating continuation assumptions from dataset coverage and numerical convergence.
- [Coherent nuclear masses/Q and two-body birth](docs/NUCLEAR_DATA.md), with primary-source provenance and explicit compatibility conventions.
- [Finite thermal-fuel burn trial](docs/THERMAL_BURN.md), with simultaneous shared-fuel depletion and reaction-conditioned energy debits.
- [Measured Maxwellian handoff API](docs/THERMAL_HANDOFF.md), with separate distribution/energy criteria and signed bath-energy correction.
- [Beam/Maxwellian window rates and unequal-temperature reactant energy moments](docs/BEAM_AND_THERMAL_MOMENTS.md).
- [Conservative birth mapping and sequential product kinematics](docs/PRODUCT_KINEMATICS_AND_MAPPING.md).
- [Unequal-mass two-body kinematics and laboratory Lorentz transformation](docs/LABORATORY_PRODUCTS.md).
- [Nuclear Coulomb functions and coherent alpha-amplitude kernel](docs/ALPHA_AMPLITUDES.md), with an explicit finite numerical domain and unresolved experimental spectrum validation.

The simplified pB instantaneous-thermalization source API remains available.
These tested primitives do not yet constitute a complete pB fast-particle
burn model: validated birth probabilities, physical ash/loss closure and
coupled reaction-state evolution remain under development. Each module
states its data domain, state ownership and validation limits.

## Requirements

- A C++17 compiler
- CMake 3.16 or newer
- Boost 1.70 or newer (Boost.Math headers)
- Optional: a Fortran compiler with ISO_C_BINDING support (gfortran and Intel
  Fortran are supported by the CMake build)
- Optional: gnuplot for the comparison plots

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For an Intel oneAPI build, initialize the compiler and runtime environment in
the same shell before configuring or running tests, for example:
`source /opt/intel/oneapi/setvars.sh`.

When a Fortran compiler is available, the same build also creates the
`pb11::fortran` target and runs the cross-language binding test.  To build the
C++ library without probing for Fortran, use
`-DPB11_BUILD_FORTRAN=OFF`.  The generated module file is placed in
`build/fortran-mod/pb11_fortran.mod` and is installed under
`lib/fortran/modules` by `cmake --install`.

To build only the library and comparison driver, disable CTest:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build -j
```

## Public API

Include the public header:

```cpp
#include "pb11.hpp"
```

The four functions are declared in `include/pb11.hpp`:

| Function | Input | Output | Published range |
|---|---|---|---|
| `pb11_sfactor(E_MeV)` | Centre-of-mass energy in MeV | S-factor in MeV·barn | 0–9.76 MeV |
| `pb11_cross_section(E_MeV)` | Centre-of-mass energy in MeV | Cross section in barn | 0–9.76 MeV |
| `pb11_reactivity_integral(T_keV)` | Ion temperature `kT` in keV | Reactivity in m³/s | Reference calculations use 10–500 keV |
| `pb11_reactivity_fast(T_keV)` | Ion temperature `kT` in keV | Reactivity in m³/s | 10–500 keV |

`pb11_sfactor()` and `pb11_cross_section()` return NaN outside the published
energy range.  The cross section returns exactly zero at zero energy.
`pb11_reactivity_integral()` accepts positive finite temperatures, while the
analytic approximation returns NaN outside 10–500 keV.  At 70 keV the fast
function uses the high-temperature branch.

## C ABI and Fortran binding

Foreign-language callers should include `include/pb11_c.h` instead of relying
on C++ name mangling.  The C ABI functions use output pointers and return an
integer status code.  They catch C++ exceptions, reject non-finite inputs, and
set the output to finite zero on every error path.  A successful call returns
`PB11_STATUS_OK` and a finite result.  Status text is available through
`pb11_c_status_message()`; the ABI revision is returned by
`pb11_c_abi_version()`.

The explicit C entry points are:

```c
int pb11_c_reactivity_integral(double T_keV, double *result);
int pb11_c_reactivity_fast(double T_keV, double *result);
int pb11_c_reactivity(double T_keV, int method, double *result);
```

`PB11_REACTIVITY_INTEGRAL` selects the direct quadrature reference path and
accepts every positive finite temperature.  `PB11_REACTIVITY_FAST` selects the
published analytic approximation and accepts 10–500 keV.  `pb11_c_sfactor`
and `pb11_c_cross_section` accept 0–9.76 MeV.  Out-of-range, non-finite,
unknown-method, null-output, and numerical errors are distinguishable through
the constants in `pb11_c.h`.

The source module `fortran/pb11_fortran.f90` wraps these functions with
`ISO_C_BINDING` and is built as `pb11::fortran` when Fortran is available:

```fortran
program reaction_rate
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use pb11_fortran
  implicit none
  real(c_double) :: rate
  integer(c_int) :: status

  call pb11_reactivity(100.0_c_double, PB11_METHOD_INTEGRAL, &
                       rate, status)
  if (status /= PB11_STATUS_OK) error stop "pB11 reaction-rate failure"
  print *, rate                 ! m^3/s
end program reaction_rate
```

The module also exports `pb11_reactivity_integral`,
`pb11_reactivity_fast`, `pb11_sfactor`, and `pb11_cross_section`, each with a
final `status` argument.  For BALDUR, convert the returned SI reactivity once
at the source-term boundary (`1 m^3/s = 1e6 cm^3/s`) and check the status
before forming a density product.  The integral and fast methods are explicit
so a caller cannot silently change the physics model by changing a global
library setting.

The two reactivity paths are intentionally independent:

```text
pb11_sfactor(E) -> pb11_cross_section(E) -> pb11_reactivity_integral(T)

Tentori equations (7)-(15) -------------> pb11_reactivity_fast(T)
```

The integral function evaluates equation (6) up to 9.76 MeV using Boost's
adaptive 61-point Gauss–Kronrod quadrature.  It explicitly resolves the narrow
148 keV resonance and performs all mass, energy, barn, and SI conversions from
their physical definitions.

## Instantaneous thermal source (SI)

Include `pb11_source.h` to evaluate fuel consumption, thermal helium production
and a specified electron/ion split in one stateless call:

```c
#include "pb11_source.h"

pb11_instant_source_v1 source;
int status = pb11_c_instant_thermal_source(
    100.0 * 1.602176634e-16, /* kT in joules, not kelvin */
    1.0e20, 2.0e19,         /* proton and B11 densities, m^-3 */
    0.3, PB11_REACTIVITY_INTEGRAL, &source);
/* Use source only when status == PB11_STATUS_OK. */
```

`reaction_rate_m3_s` is R = np*nB*<sigma v>. The proton, boron and
helium-4 source fields are respectively -R, -R and 3R in m^-3 s^-1.
`fusion_power_W_m3` uses Q=8.68 MeV per reaction. Electron and ion power
fields sum to that released power. The electron fraction is prescribed by the
caller; this interface assumes immediate local thermalization and does not
calculate slowing, escape or reactant thermal-energy redistribution.

The Fortran module exports the same `bind(C)` derived type and
`pb11_instant_thermal_source(kt_j,np,nb,fe,method,source,status)` wrapper.
Inputs must be finite with kT>0, densities>=0 and 0<=fe<=1; the selected
reactivity method retains its temperature domain. All output fields are zero
on failure. Unrepresentable nonzero results return a numerical error.
For CGS hosts, multiply number densities by 1e6 before this call, multiply
returned particle sources by 1e-6 and powers by 10 to obtain cm^-3 s^-1
and erg cm^-3 s^-1. The original scalar interfaces remain compatible.

## Minimal example

```cpp
#include "pb11.hpp"

#include <iomanip>
#include <iostream>

int main() {
    const double energy_MeV = 0.600;
    const double temperature_keV = 100.0;

    std::cout << std::scientific << std::setprecision(8)
              << "S(E)       = " << pb11_sfactor(energy_MeV)
              << " MeV barn\n"
              << "sigma(E)   = " << pb11_cross_section(energy_MeV)
              << " barn\n"
              << "integral   = " << pb11_reactivity_integral(temperature_keV)
              << " m^3/s\n"
              << "analytic   = " << pb11_reactivity_fast(temperature_keV)
              << " m^3/s\n";
}
```

When this repository is included from another CMake project:

```cmake
add_subdirectory(path/to/pB11_reaction)
target_link_libraries(your_target PRIVATE pb11::pb11)
```

## Test and comparison drivers

`test_pb11` checks the 148 keV resonance, the 0.400 and 0.668 MeV piece
boundaries, zero/invalid inputs, representative Figure 1 magnitudes, the
10–500 keV reactivity accuracy, and the 70 keV LT/HT switch:

```bash
./build/test_pb11
```

`benchmark_pb11` prints the numerical integral and analytic approximation at
the key temperatures used for validation, followed by the maximum errors in
the low- and high-temperature regions:

```bash
./build/benchmark_pb11
```

It can also produce machine-readable CSV scans:

```bash
./build/benchmark_pb11 --energy-scan > energy_scan.csv
./build/benchmark_pb11 --reactivity-scan > reactivity_scan.csv
```

The energy scan contains `E_MeV,S_MeV_barn,sigma_barn`.  The reactivity scan
contains the temperature, both reactivities, signed relative error, and
absolute relative error.

## Plotting

The supplied gnuplot scripts can invoke the benchmark directly:

```bash
gnuplot -c examples/plot_energy_scan.gnuplot
gnuplot -c examples/plot_reactivity_comparison.gnuplot
```

They create `pb11_figure1.png` and `pb11_reactivity_comparison.png`.  To plot a
saved CSV and select the output filename, pass both as arguments:

```bash
gnuplot -c examples/plot_energy_scan.gnuplot energy_scan.csv figure1.png
gnuplot -c examples/plot_reactivity_comparison.gnuplot \
    reactivity_scan.csv reactivity.png
```

With the published, rounded table parameters, the analytic reactivity differs
from the direct integral by about 2% below 70 keV and less than 1% from 70 to
500 keV, consistent with the paper.

## License

This project is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE)
for the complete terms.

Explicit continuation and per-segment thermal/beam rate/energy APIs: [RATE_MODEL.md](docs/RATE_MODEL.md).

Published alpha-spectrum comparison and convergence limitations: [PUBLISHED_SPECTRUM_CHECK.md](docs/PUBLISHED_SPECTRUM_CHECK.md).

Accepted/trial kinetic contexts, explicit source/heat ledgers and portable
restart are documented in [SOURCE_STATE](docs/SOURCE_STATE.md). They wrap
physical operators and do not by themselves implement host coupling. The
[Taskaev2024 audit](docs/TASKAEV_2024_CHANNEL_EVIDENCE.md) extends the
channel-resolved pB evidence with target-energy and angular-model caveats.

[Explicit FSCI source policy](docs/FSCI_SOURCE.md) adds the documented
Refsgaard Model-II correction at16fm and source-shape refinement to1024
quadrature nodes. Conditional CM source sensitivity is validated separately
from incident-branch/laboratory composition and host feedback.
