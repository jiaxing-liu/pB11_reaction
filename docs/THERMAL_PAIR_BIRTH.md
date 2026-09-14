# Unequal-temperature thermal product birth

The new `fusion_c_thermal_pair_birth_grid` extends the common-temperature source
with two independent, zero-drift Maxwellian reactant pools. It returns full
scalar-energy product spectra plus explicit spill and neutron coefficients.
The source remains stateless, uses canonical channel order/masses/Q and the
existing nuclear models, and does not advance fuel or host state. Both kT values
are positive joules. No density or identical-pair counting factor is included.

## Correlated Gaussian reduction

This is the Gaussian-conditioning derivation already documented and tested for
scalar reaction-conditioned energies in BEAM_AND_THERMAL_MOMENTS.md, now retained
inside the birth integration. Write u=va-vb and V=(ma*va+mb*vb)/M, M=ma+mb and
mu=ma*mb/M. The relative Maxwellian temperature and conditional CM distribution are

```
Tr = mu*(Ta/ma + Tb/mb)
a  = (Ta-Tb)/(M*(Ta/ma + Tb/mb))
V  = a*u + W
Var(W_i) = Ta*Tb/(ma*Tb+mb*Ta)
Tc = M*Var(W_i) = Ta*Tb/Tr.
```

W is independent of u. For each reacting relative energy Er, integrate the
residual CM energy M*|W|^2/2 and polar cosine between W and u. Then reconstruct
va=V+(mb/M)*u and vb=V-(ma/M)*u. Their actual momenta enter the existing parent
API under the classical reactant energy convention. The parent supplies total
energy, available CM product energy and the full vector boost. The source uses
the boost magnitude with its existing isotropic outgoing CM-event orientation.
This construction retains the correlation; using Tr alone with an independent
CM Maxwellian would reproduce scalar rate while generally giving wrong fuel
debits and product spectra.

The independent reference rate/debits use the existing scalar model, whose
conditional energy moments are analytically reduced. The new source instead
integrates both reconstructed reactant energies at each velocity node. Product
number and Q+both-reactant-energy closure include all grid spills and neutrons.

Equal temperatures dispatch through the original arithmetic path: no new angular
quadrature is applied. The public old function and struct layouts remain intact.
The correlation order is validated4..32 even when equal temperatures make it
inactive. The new Fortran wrapper uses explicit C kinds and the existing(cells,7)
array layout.

## Explicit model and numerical limits

The nuclear source and isotropic outgoing orientation are exactly the declared
common-T closures. This does not add a calibrated angular distribution, beam
pitch, anisotropic incoming populations, or relativistic nuclear amplitudes.
pB alpha1 still retains NR amplitude weights at Q+Er followed by the documented
on-shell kinematic remapping; alpha0 retains sequential narrow-width kinematics.

`relative_max_J` bounds Er. `cm_max_kT` now bounds the residual W energy in units
Tc, not total CM kinetic energy. Existing cm tail fields explicitly describe W;
they are not missing-product-energy bounds. Finite tails are not renormalized.
A conservative maximum of |a*u|+|W| supplies the pB parent-domain check.
Independent relative-energy, residual-speed, polar-angle, nuclear-event and
output-grid convergence remain caller obligations. An OK status proves finite
conservative evaluation, not universal numerical or experimental accuracy.

This API does not activate unequal-T in the current common-ion-temperature host
controller or table. Fast-target conditional birth, shared thermal-target burn
competition and exactly-once nuclear born/fast consumed ledgers remain separate
required implementation steps. An energy moment is never substituted for that
missing product distribution.

## Validation of this increment

- Four frozen pre-change O3 pB evaluations at0.002/0.006/0.1/2keV retain exact
  complete-output byte parity after refactoring. This compares the old executable,
  not merely two new functions using the same implementation.
- Existing GNU common-T birth/width/table regressions pass3/3.
- New C++ and Fortran tests pass GNU2/2 and Intel-r8 3/3 (including old binding).
  Tests cover all-channel equal-T parity, unequal DT/DD/DHe3 independent rates
  and both debits, explicit spills and conservation, DD temperature exchange,
  near-equal continuity and invalid-input clearing. Mapped number/energy tests
  use3e-10 relative tolerance; low-order rate/debit tests separately use2e-4/5e-4.
- Root source study uses128 uniform cells0..25MeV, relative16/residual12,
  nuclear-event8/8 and angles4/8/16/32. DT20/40keV and pB100/200keV have final
  rate/debit discrepancies6.96e-12/4.08e-11 and6.16e-9/3.35e-8 respectively.
  Angle16->32 number/energy L1 differences (including spills) are
  DT9.48044e-6/1.00153e-5 and pB8.07022e-7/6.10889e-7. This is representative
  polar convergence, not validation of all other quadrature/physical settings.
- The initial study used1000*keV to spell MeV, which rounded the one-keV cutoff
  below the existing exact lower bound; status3 was retained. Correct explicit
  MeV spelling resolves input failure without changing the gate.
- Initial build invocation omitted CMake regeneration and found no new targets;
  no-tests output was not accepted. The first regenerated Intel build selected
  the C++ linker and failed PIE relocation. Matching the established Fortran
  LINKER_LANGUAGE property fixes linkage. Final commands use no-tests=error.

Evidence is under docs/validation/thermal-pair-birth. Coupled multi-temperature,
monoenergetic fast-target births and full case acceptance are not claimed.
