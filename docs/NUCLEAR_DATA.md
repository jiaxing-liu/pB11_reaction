# Coherent nuclear masses and deterministic two-body birth

`fusion_nuclear_data.h` defines the fixed mass convention for the new composed
reaction model. It supplies bare nuclear masses, explicit reactant/product IDs,
Q derived from those same masses, and a deterministic two-body product call.
It does not change the compatibility APIs `fusion_c_channel_q`,
`fusion_c_thermal_rates`, or the old instantaneous pB source. Those retain their
previous table, including the expressly rounded pB8.68MeV. A caller must use
one convention consistently throughout its reaction/product energy ledger.

## Source selection and conversion

The root independently checked these primary sources on2026-09-14:

- [NIST CODATA2022 complete listing](https://physics.nist.gov/cuu/Constants/Table/allascii.txt):
  relative nuclear masses in u for p,D,T,He3,He4,n and the electron; the common
  atomic mass constant1.66053906892e-27kg and quoted standard uncertainties.
- [AME2020 mass table, IAEA mirror](https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt):
  B11 row72 is the NEUTRAL atomic mass11.009305166(13)u. Its mass excess is
  8667.708(12)keV. The micro-u table column is converted explicitly.
- [NIST ASD B ionization energies](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&format=0&ion_charge_out=on&ion_conf_out=on&level_out=on&order=0&seq_out=on&shells_out=on&sp_name_out=on&spectra=B&submit=Retrieve+Data&unc_out=on&units=1):
  the five successive energies sum to670.9838405eV. The underlying dataset is
  ASD5.12(2024), DOI10.18434/T4W30F. This is an element-table binding model,
  not an independently published isotope-specific B11 total binding energy.

The detailed source audit, table values and definition of atomic versus nuclear
Q are in [NUCLEAR_MASS_SOURCE_AUDIT.md](NUCLEAR_MASS_SOURCE_AUDIT.md). The new
set uses the more precise CODATA relative mass values, not separately rounded
MeV mass equivalents. Every kg mass uses the same atomic mass constant u.
With exact SI c=299792458m/s and e=1.602176634e-19C,

```
m_B11 = (M_B11_atom/u - 5*m_e/u)*u + E_binding/c^2
      = 11.006562986785523... u
Q     = (sum m_reactants - sum m_products)*c^2.
```

The binding term has a PLUS sign: the neutral atom's rest mass is smaller than
the separated bare nucleus plus five electron masses. Masses are rounded once
to canonical returned doubles; Q and rest energies are computed from those
same kg values in long-double arithmetic. Thus the kinematic inputs and nuclear
ledger cannot disagree because different Q/mass tables were selected. Table
rounding is far below the measurement uncertainty but does not justify hiding
inconsistent input conventions.

## Known-data uncertainty and model limitation

The `known_uncertainty_scale` fields propagate the published standard
uncertainties by a linear absolute-derivative sum. This is a conservative
standard-deviation scale when correlations are unavailable, NOT a probability
confidence interval or hard bound on the true physical value. For Q, the
shared u uncertainty multiplies the relative MASS DIFFERENCE; treating u as
an independent error for each full nuclear mass would lose that correlation.
Repeated reactants/products use their signed integer stoichiometric weights.

The B binding contribution uses the linear sum of the five reported ASD
ionization-energy uncertainties,0.0001326eV. That does NOT quantify the
isotope/model systematic of using an element table for B11. This limitation
is explicit, not a zero uncertainty claim. The AME B atomic-mass uncertainty
alone is about12eV; retaining an unquantified binding-model systematic is
adequate for declaring the chosen model's scope, not a metrology certification.

Independent45-digit decimal calculations before canonical kg rounding give:

| Channel | Q(MeV) | Known-data linear uncertainty scale(eV) |
|---|---:|---:|
| pB→3alpha | 8.68237827423 | 12.29331 |
| DD→T+p | 4.03266389874 | 0.13009 |
| DD→He3+n | 3.26885695083 | 0.47050 |
| DT→alpha+n | 17.58924790374 | 0.54298 |
| DHe3→alpha+p | 18.35305485165 | 0.15414 |

The table is derived from this chosen set, not five independent experimental
Q measurements. The old pB8.68MeV model differs by about2.378keV/event. Atomic
Q values also differ from bare nuclear Q by electronic binding; neither is
an extra heat source to add on top of the nuclear mass deficit.

Mass-data IDs0..5 coincide with the six thermal species. IDs6(neutron) and
7(electron) are data IDs, not new thermal species slots. Nuclear charge means
bare Z, not the coronal mean charge used by `natomc=2`. Bound-electron charge
and atomic radiation remain the host atomic model's responsibility.

## Product composition and validation boundary

`fusion_c_nuclear_two_body_cm` accepts a DD, DT or DHe3 channel, the TOTAL
reactant kinetic energy in its center-of-momentum frame, and a unit direction
for the first listed product. It computes Q from the canonical masses and
calls the existing exact relativistic two-body kinematics with Ecm+Q. It does
not choose a differential cross section or angular sampling probability.
pB is rejected by this two-product API: its separate three-alpha source is
still required. There is no equal-alpha substitute in this operation.

Root's24-case composition check reconstructs incoming and outgoing particles,
then applies the same Lorentz boost to both sides. Across all four two-body
channels, Ecm=0,.01,1MeV and zero/nonzero boosts, lab-frame kinetic energy gain
equals Q and total momentum is conserved. Maximum relative Q/momentum
residuals are5.84e-16/2.01e-17. This is a reaction-event kinematic check;
thermal/beam distribution averaging, source probabilities, slowing and host
coupling remain separate work. Raw CSV and decimal reference are retained in
`docs/validation/nuclear-data/`.

## Interface acceptance and reproduction

GNU/Intel default/Intel-r8 each pass36/36 CTest cases; no-Fortran passes23/23.
The installed C11 example and BALDUR Intel-r8 rebuild/link pass. The new
unit tests independently check decimal Q/uncertainty values, B electron-mass
and binding conversion, charge/baryon counts, the old pB Q compatibility value,
invalid data IDs and malformed two-body inputs. The Fortran wrappers use
explicit c_double/c_int, checked extents, and C-interoperable mass/channel
records; the product record is reused from fusion_laboratory_fortran.

```
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j4
ctest --test-dir build --output-on-failure
build/validate_nuclear_birth > nuclear-birth-four-momentum.csv
```

Fortran consumers link exported target `pb11::fusion_nuclear_data_fortran`
and use module `fusion_nuclear_data_fortran`; use ONLY imports when combining
modules with shared status/species identifiers. The product record can be
imported from `fusion_laboratory_fortran`. A standalone C11 caller is in
`examples/nuclear_data_consumer.c`. No runtime web/database dependency exists;
this is an explicitly frozen source set. The new operators are not yet active
in BALDUR's source assembly, and this checkpoint does not certify EXL cases.
