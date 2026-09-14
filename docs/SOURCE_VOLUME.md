# Volume-aware accepted kinetic state

The three `fusion_c_source_state_*_volume` entry points extend the existing
source-state transaction manager. They do not advance physics or infer mesh
motion. The same begin/discard/commit/commit_many/pack/unpack/destroy functions
operate on the new opaque contexts. Legacy fixed-volume contexts and v1/v2
restart layouts remain supported. Legacy stage and snapshot calls reject new
volume contexts, avoiding silent changes to their density interpretation.

## Explicit normalization

Let V0 be the immutable initial/reference volume, Va the accepted actual volume,
Vb the trial actual volume and Vs the explicit source normalization volume.
S and T remain **physical number densities** in each energy bin. Both are
kinetic populations; neither is thermal ash. Snapshot returns V0 and Va.

The existing 121-field cumulative source ledger and six signed inert-bath heat
fields are stored **per reference m3**: a supplied source-density increment is
multiplied by Vs/V0 before accumulation. An extensive source amount is recovered
by multiplying the cumulative field by V0, never by the current Va. The source
ledger retains stoichiometric/Q and per-species validation. Invalid original
source fields are checked before scaling; overflow or a nonzero amount rounded
to zero by this conversion returns an error. There is no population floor.

The new 42-double `fusion_transport_ledger_v1` is **extensive**, in particle
counts and joules:

| Array (six species each) | Sign and meaning |
|---|---|
| spatial_number, spatial_energy_J | signed net inward radial/coordinate exchange |
| work_J | signed mechanical work onto the kinetic populations |
| lower_number, lower_energy_J | nonnegative outgoing numerical lower-grid inventory |
| upper_number, upper_energy_J | nonnegative outgoing numerical upper-grid inventory |

Energy-grid boundary energy must equal the corresponding edge energy times
its outgoing number. These channels are distinct from physical escape,
external injection, nuclear birth, bath heat and fluid handoff. Zero net spatial
number with nonzero net spatial energy is legitimate when energy-dependent
counterflows occur. The context only records moments of grid outflow; it does
not evolve a kinetic population beyond its energy grid or provide an automatic
re-entry/thermalization closure. Significant outflow requires domain refinement
or a separately specified physical treatment.

## Checked balances and composition

For each species, define N and U as moments of S+T using midpoint energies.
Writing B for nuclear plus external births, F for fast reactant consumption,
L for physical escape, H for handoff, Qb for signed heat to all baths, X for
net inward spatial exchange, W for work and O for numerical domain outflow:

```
Vb Nb - Va Na = Vs (B_N - F_N - L_N - H_N) + X_N - O_N
Vb Ub - Va Ua = Vs (B_U - F_U - L_U - H_U - Qb) + X_U + W - O_U
```

Both step and cumulative identities are checked on the common V0 normalization,
using the existing relative 1e-10 absolute-term scale without a unit-sized
floor. All ledger accumulations must remain representable. Cumulative transport
is extensive; cumulative source and heat are per V0. No nuclear, heat or work
term is manufactured to force closure.

If source suboperators use different volumes, the caller must first express
their extensive source increments at one explicit common Vs. For radial then
energy-work splitting, radial face amounts are already extensive. The work
operator returns density amounts at its input population's volume; multiply
those work and domain ledgers by that volume before staging. A shared radial
face enters the adjacent zone ledgers with opposite signs. The test composes
the actual existing radial and work kernels and checks both local and global
budgets. It prescribes the work coefficients explicitly; it is not a claim that
arbitrary mesh motion determines a physical compression model.

## Acceptance and restart

A failed current-ticket replacement stage invalidates any earlier candidate.
No trial volume is accepted until commit. `commit_many` validates every ticket
before publishing any zone's vectors, volume or cumulative ledgers, with no
allocation during publication. Host fluid/geometry must be accepted at the same
epoch through the host's own transaction; this API cannot roll back host state.

Volume contexts serialize as version 3. After the v2 six inert-heat words, v3
adds reference volume, accepted volume, and 42 transport words in C declaration
order, before edges/S/T. Length is `8*(193+13*cells)` bytes. Fresh epoch-zero
volume contexts may be packed, with zero source/heat/transport and Va=V0.
Unpack checks checksum, dimensions, model tag, semantic fields and cumulative
balances. No pending trial enters a restart. v1/v2 layout is unchanged; readers
without v3 support reject it. Model-tag configuration and host checkpoint
coordination retain the existing SOURCE_STATE.md contract.

## Validation and limits

`test_fusion_source_volume` covers an exact two-bin sequence: volume dilution,
signed work, outward spatial flow, upper-grid outflow, differing source/actual
volumes and inert heat. It checks restart byte identity and identical continued
execution, checksum-correct semantic corruption rejection, batch atomicity,
illegal replacement stages and legacy API guards. A separate SI-scale synthetic
DT ledger checks normalized event/Q accounting; its test energy split is not a
DT product model. The actual radial+work composition checks internal-face
cancellation and independent extensive N/U balances.

The companion Fortran module is `fusion_source_volume_fortran`. Use the base
`fusion_source_state_fortran` lifecycle for its c_ptr handles. No BALDUR indices,
common blocks, file I/O, transport coefficients or physical-loss defaults are
introduced. This is a required accounting layer, not completion of moving-grid
BALDUR nuclear runs, global host energy accounting or full-host restart.

Final focused evidence: GNU/default and Intel Fortran `-r8 -check all` each
pass all five state/volume tests. Two consumers built against the installed
export pass (C11 and Fortran). See `validation/source-volume/`. The current
installation exports `pb11Targets.cmake`; it does not supply a `find_package`
Config file. The consumer uses that installed target export explicitly.
