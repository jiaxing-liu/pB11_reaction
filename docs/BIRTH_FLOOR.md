# Candidate numerical treatment for below-grid births

This is a provisional numerical projection, not a new pB nuclear model and not
an installed BALDUR repair. Existing source and coupled APIs still reject
representable charged-product spill. Do not call the full pB case accepted.

For below-first-center birth amounts N>=0 and U>=0, retain N in the kinetic
species and place it at the existing first cell center E0. Define:

    U_grid = N E0
    delta_U_ion = U - U_grid
    U_grid + delta_U_ion = U

The physical nuclear birth ledger keeps U. The numerical kinetic source uses
U_grid. The caller must publish delta_U_ion as a separately identified signed
numerical energy transfer exactly once, together with the source; neither
nuclear Q nor Coulomb heating may silently absorb it. No particles become
thermal ash in this operation. Upper-grid overflow remains unsupported.

Inputs are amounts over the accepted candidate interval, not source rates or
conservative full-consumption bounds from D102. The library operates on all
six canonical charged species together so the ion energy budget is shared.
It checks both the unrounded and returned total correction and positive
remaining reservoir energy. It must not change or scale the particle counts.

Two caller-supplied acceptance limits constrain E0/kTi and the TOTAL numerical
energy borrowed relative to the ion reservoir energy. These bound numerical
energy displacement only. They are not error bounds on the spectrum, physical
transport, source quadrature or long-time evolution. No default limits are
selected by this utility. The caller certifies these moments come from the
below-center source; two moments alone cannot establish its support.

Choice of an ion numerical reservoir is explicit in the interface. It is not
a claim that physical cold-particle energy exchange acts only on ions. A
coupled promotion must record this numerical contribution separately, check
its magnitude, examine electron/ion energy sensitivity, and compare refined
grids/time steps without loosening the existing nuclear/FP/handoff gates.

The existing Coulomb kernel is based on the exact relaxation rates in the
2019 NRL Plasma Formulary, p31:
https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Plasma_Formulary_2019.pdf
Its low-energy limit includes energy diffusion/heating, so low energy alone
does not establish a Maxwellian distribution suitable for whole-population
handoff. The algebraic numerical projection above is OUR explicit numerical
choice; that reference does not validate it. Existing Maxwellian handoff
L1/mean-energy criteria remain unchanged.

Required before coupled promotion: C/C++/Fortran ABI and invalid-output tests;
source ledger inclusion of physical and mapped birth energy separately;
transaction rollback; finite reservoir handling; exact zero-spill parity;
actual pB retry/full-window tests; lower-edge/time/grid sensitivity; explicit
uncertainty limits. A standalone utility pass is not host/physics acceptance.

## Opt-in coupled source entry

`fusion_c_coupled_sources_floor_trial` aggregates actual thermal and fast-target
birth amounts after reactant debits and projects all six charged species jointly.
It returns the separate floor ledger alongside the physical source ledger. The
ion correction updates the trial reservoir before collision baths are assembled;
its temperature gate explicitly uses post-reactant-debit PRE-correction kTi.
Upper spill remains subject to the existing strict representability bound.
Old source entry points retain their behavior. A packet that becomes a zero
rate rejects; conversions are checked against measured source-scale IEEE rounding.

This is still an experimental numerical boundary, not validated physical closure.
The existing source-state stage/snapshot/restart APIs and fluid-increment helper
do NOT own this additional correction. Callers must not stage a nonzero-floor
result through those legacy APIs or reinterpret the correction as heat. Full
state-owner/host integration and convergence are pending. A caller using
this pure trial owns atomic publication of ALL returned accounts.

The coupled-floor C++ test uses a deliberately artificial local plasma to make
the correction representable:1keV first center,10keV ions, fastH1e20m^-3 and
thermalB1e19m^-3,dt1ms. It checks exact zero-spill parity, explicit-gate rejection
and identical retry, and independently reconstructs ion energy from the physical
ledger plus the separate correction. It is not an EXL actuator prescription.

The direct `fusion_coupled_sources_floor_fortran` ISO_C_BINDING module uses
explicit `c_ptr,value` pointers and interoperable types imported from existing
modules. Null optional arrays and failure clearing are tested with GNU/Intel.
