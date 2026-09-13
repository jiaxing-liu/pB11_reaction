# Thermal reaction network: scope and validation

The new C APIs are `fusion_network.h`, `fusion_rates.h` and
`fusion_cross_sections.h`. They are standalone and use SI units. The existing
pB C++/C/Fortran APIs remain compatible. `fusion_fortran` provides the new
Fortran interface separately from `pb11_fortran`.

Species IDs are p,D,T,He3,He4,B11. Channel IDs are pB->3alpha,
DD->T+p, DD->He3+n, DT->alpha+n, DHe3->alpha+p. Each DD branch has its own
reactivity. Its event rate is nD^2 <sigma v>/2; two deuterons are destroyed per
event. This implementation does not impose a 50:50 DD branch approximation.

All reactants in `fusion_c_thermal_rates` have one common Maxwellian kT.
The densities provided are thermal reactant densities, not fast-ion populations.
A mask of 30 enables both DD branches and the DT/DHe3 secondary channels.
Product birth is separate from reactant loss and is not yet thermalized ash.
The returned net source is their algebraic difference, useful for isothermal
instant-product network tests. No function here advances physical time.

## Nuclear models and domains

DT/DD/DHe3 use Bosch and Hale 1992, Nuclear Fusion 32 611,
[DOI](https://doi.org/10.1088/0029-5515/32/4/I07).
The coefficient audit and original-paper mirror are in [BOSCH_HALE_DATA.md](BOSCH_HALE_DATA.md).

| Channel | Fitted thermal kT (keV) | Cross-section CM energy (keV) |
|---|---:|---:|
| DD -> T+p | 0.2-100 | 0.5-5000 |
| DD -> He3+n | 0.2-100 | 0.5-4900 |
| DT | 0.2-100 | 0.5-4700 |
| DHe3 | 0.5-190 | 0.3-4800 |
| pB | Existing integral / FAST (FAST: 10-500) | 0-9760 |

New APIs receive energy in joules, not these tabulation units. Cross sections
return m^2; reactivities return m^3/s. Exact zero collision energy gives zero
cross section. Positive energies below a nuclear-data boundary fail explicitly.
No thermal fit is extrapolated beyond its range; even a selected zero-density
channel is checked. A one-ULP conversion-edge correction is documented.

DT changes from the low to high cross-section fit at 530 keV, following the
paper's p622 recommended intersection rather than silently extending the low
fit. DHe3 changes at 900 keV and retains the published small discontinuity.
Table VI lists only A1 in the numerator; unused higher A terms are zero, not
inherited from the low fit. This is consistent with p622's statement that not
all coefficients are needed and is checked against the DT intersection and
DHe3 positive jump. Fits are not refitted to make them artificially smooth.

## Nuclear energy ledger

`fusion_c_channel_q` returns nuclear release per event in J. The thermal
network also returns per-channel and total RQ in W/m^3. These values are not
charged-product energy, neutron energy or deposited electron/ion heat.
Reactant kinetic energy is not included in RQ.

For the four two-body channels, Q is calculated from the nuclear rest-energy
values in the [NIST CODATA 2022 listing](https://physics.nist.gov/cuu/Constants/Table/allascii.txt):
p=938.27208943, D=1875.61294500, T=2808.92113668,
He3=2808.39161112, He4=3727.3794118 and n=939.56542194 MeV.
Using nuclei avoids confusing atomic electron binding with nuclear release.
pB retains the explicit rounded 8.68 MeV value of the existing validated
instant-source baseline. The MeV-to-joule conversion is exactly 1.602176634e-13.
Constants are shared within the library (`src/fusion_constants.hpp`).

Full kinetic coupling must later account for kinetic energy of selected
reactants, product birth spectra, fast inventory, thermalization and escape.
It must not add RQ as thermal heat on top of separately deposited product energy.

## Verification and limits

- Table V published cross-section values at 10/50/100 keV test units and formula
  transcription. Table VIII thermal values at 1/5/10/20/50 keV test rate fits.
  These are publication regression checks, not new experimental validation.
- Independent Gauss-Kronrod cross-section integration versus the separately
  fitted thermal rates at 1/5/10/20/50/100 keV gives maximum relative differences
  of 1.9088% (DD Tp), 1.0930% (DD He3n), 0.7483% (DT), 3.1454% (DHe3).
  Both fits have published approximation errors; this is not integration error.
  Integration is over the documented cross-section domain, not undocumented
  extrapolation to zero/infinite energy. These selected temperatures suppress
  omitted tails; the comparison does not certify arbitrary-temperature moments.
- Channel tests check nucleons including neutron births, charge, independent
  birth/loss arrays, finite outputs, masks, density scaling, DD symmetry and errors.
- The closed-box test uses an explicitly isothermal, instant-product network
  with an RK4 test integrator. Single-channel depletion/yield agrees with analytic
  solutions; DD secondary burning converges under 80/160/320 steps. Nuclear
  release integrated from powers agrees with event counts times Q. This is not
  a self-heated zero-dimensional kinetic plasma or a BALDUR transport test.

Build normally with CMake and run CTest. The new tests are
`fusion_network_tests`, `fusion_rates_tests`, `fusion_cross_sections_tests`,
`fusion_closed_box_tests`, and `fusion_fortran_binding_tests` (when available).

Not yet implemented here: unequal-temperature or non-Maxwellian rates,
product spectra, collisions, finite-time ash, deposition or transport. TT,
He3-He3 and rare pB side channels are outside this initial network and must not
be inferred from the word network. The full EXL-50U program remains in progress.

## Independent caller example

`examples/thermal_network.cpp` calls only the public C ABI and is built as
`thermal_network`. Run `./build/thermal_network` after the normal CMake build.
It prints separate branch events and nuclear release for a 10 keV thermal DD
input, with no BALDUR dependency or global state.

Fortran species/channel IDs have the same zero-based numerical values as C.
The Fortran arrays use normal one-based subscripts, so index them with `ID+1`.
