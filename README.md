# pB11_reaction

A small, standalone C++17 library for the thermal proton–boron-11 fusion
reaction.  It implements the astrophysical factor, fusion cross section, and
Maxwellian reactivity described by:

> A. Tentori and F. Belloni, “Revisiting p-11B fusion cross section and
> reactivity, and their analytic approximations”, Nuclear Fusion 63 (2023)
> 086001. [doi:10.1088/1741-4326/acda4b](https://doi.org/10.1088/1741-4326/acda4b)

The library uses `double` precision throughout and is intentionally limited to
the S-factor, cross section, and thermal reactivity.

## Requirements

- A C++17 compiler
- CMake 3.16 or newer
- Boost 1.70 or newer (Boost.Math headers)
- Optional: gnuplot for the comparison plots

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

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

The two reactivity paths are intentionally independent:

```text
pb11_sfactor(E) -> pb11_cross_section(E) -> pb11_reactivity_integral(T)

Tentori equations (7)-(15) -------------> pb11_reactivity_fast(T)
```

The integral function evaluates equation (6) up to 9.76 MeV using Boost's
adaptive 61-point Gauss–Kronrod quadrature.  It explicitly resolves the narrow
148 keV resonance and performs all mass, energy, barn, and SI conversions from
their physical definitions.

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
