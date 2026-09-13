# EXL-50U reusable fusion library development

Active work branch: `feature/exl50u-fusion-library`, based on `ede11c6b57061aef765c9b97343ddb6464fcb4ee`.

The approved program is coordinated in the BALDUR repository at `docs/exl50u-program/PLAN.md`. Its goal is a host-independent thermal/nonthermal reaction and fast-particle library plus four-fuel EXL-50U comparisons at electron chord-average density 0.6 nG. Instantaneous thermalization is only a regression approximation.

First increment: add an explicit SI instant-thermal source C ABI and Fortran binding while retaining all existing reactivity and cross-section entry points. The source API receives thermal energy kT in joules and densities in m^-3; output particle rates are m^-3 s^-1 and power is W m^-3. This moves already-used source arithmetic out of the BALDUR adapter; it does not claim a new product spectrum, slowing-down or orbit-loss model.

Further increments require documented literature/domain choices, source ownership, state rollback and particle/energy closure. No BALDUR common blocks, species slots, CGS conventions or implicit file output are part of the independent API.

M2 thermal kernel increment implemented: six charged species, five explicit
channels with separate DD branches and secondary DT/DHe3, Bosch-Hale cross
sections and common-temperature Maxwellian fits, particle birth/loss and
nuclear-Q ledgers, and C/Fortran interfaces. See THERMAL_NETWORK.md and
validation/m2-*.txt. GNU and Intel r4/r8 tests cover this stateless increment;
full finite-time product kinetics and BALDUR coupling are still outstanding.

M3 first kinetic increment: independent Maxwellian-bath Coulomb coefficients
and an isotropic energy-distribution implicit trial step with per-bath heat,
particle/energy ledgers and no hidden accepted state. See KINETIC_PRIMITIVES.md.
This is a primitive layer; product spectra, reaction/kinetic closure, physical
ash transfer and BALDUR integration remain under development.
