# Accepted kinetic state, trial ledger and restart

`fusion_source_state.h` owns a per-zone accepted kinetic inventory. It is a
transaction manager around explicit physical operators, not a new burn or
transport solver. Products enter fast populations; this context does not turn
reaction yield into thermal ash. Both S and T components remain kinetic.

## Ownership and units

The common energy grid is in J. Each component is a six-species-major array
of integrated cell numbers in m^-3; cell energy is the arithmetic midpoint.
Species and five event-channel IDs are those of `fusion_network.h`. The
immutable `cells` query permits safe buffer allocation after unpack. Different
contexts have no shared mutable state. The same context must not be accessed
concurrently. A handle has one owner and must be destroyed exactly once.

The `fusion_source_state_fortran` module uses explicit ISO_C_BINDING kinds.
`c_int64_t` carries the bit patterns of C uint64 tags/tickets/epochs without
signed range restrictions. It queries the actual cell count before checking
snapshot/stage extents. A failed stage extent check invalidates any earlier
candidate for that ticket through the C API; it cannot leave stale results
committable. Destroy clears the Fortran handle; copying a handle does not
create another owner.

A step ledger contains AMOUNTS, never rates: numbers m^-3 and energies J/m^3.
It separates nuclear and external births, consumed thermal and fast fuel,
escape, fluid handoff, seven-bath heat, and neutron number/energy. Heat is
positive to the receiving bath and may be negative. Bath order is electron,
then p/D/T/He3/He4/B11; storage is fast-species-major `[6][7]`. Every other
entry is finite and nonnegative. Zero particles cannot carry nonzero energy.
External injection or host influx needs its own explicit number/energy entry;
physical escape or host outflux likewise needs both entries. Neither numerical
grid spill nor a failed operator call is automatically classified as escape.

Host thermal state, geometry, transport fluxes and coefficients remain host
owned. The caller checks fuel and bath energy availability, and must commit
or restore the host and library at the same accepted epoch/time. The library
cannot make a separately updated BALDUR common block transactional.

## Trial protocol

1. `begin(dt)` returns a fresh context-local ticket and supersedes old trials.
   Even an invalid replacement begin clears an earlier pending trial. `dt`
   must be positive, finite, and advance representable accepted time.
2. Advance physical operators from the **accepted** state. `stage` replaces
   the candidate populations and ledger for that ticket. Repeated nonlinear
   evaluations do not sum sources. Failed staging with the current ticket
   invalidates a previously valid candidate. A stale ticket cannot affect
   the current trial.
3. Accept with `commit` once. Only then advance time/epoch and cumulative
   diagnostics. Repeat or stale commits fail. Alternatively `discard` leaves
   the accepted state and cumulative diagnostics unchanged.

The manager has no implicit physical advance, file I/O or hidden global state.
It does not implement RESOLV's extrapolated trajectories; host I007 still
requires a compatible staging protocol or explicit rejection of unsupported
extrapolation in the new stateful mode.

## Checked identities

Nuclear stoichiometry and Q use `fusion_c_nuclear_channel`, whose mass set is
documented in [NUCLEAR_DATA.md](NUCLEAR_DATA.md). Legacy rounded-Q APIs must
not supply the new nuclear energy ledger.

- Nuclear charged births and neutrons equal the event stoichiometry.
- Thermal plus fast reactant consumption equals the same event stoichiometry.
- Nuclear charged-product energy plus neutron energy equals removed thermal
  and fast reactant kinetic energy plus the event-weighted new Q.
- For each fast species: final N = initial N + nuclear births + external
  births - fast consumption - escape - handoff.
- The analogous kinetic-energy identity also subtracts signed heat to all
  baths. S-to-T transfer is internal and cancels in these identities.

Every identity is checked to relative residual 1e-10, scaled to the absolute
terms in that identity without a unit-sized absolute floor. Internal sums use
long double; stored populations/ledgers and portable restart anchors use
binary64. Each prospective cumulative ledger must be representable and pass
the same balances against the original inventory before a trial is accepted.
This detects incompatible composition, but does not certify a physical model:
a self-consistent invented source could conserve and still be wrong.

## Restart schema 1

Packing is forbidden while a trial is pending. The byte stream is explicitly
little-endian: uint64 magic, version, cell count, caller model tag, accepted
epoch and ticket counter; binary64 accepted and initial time; initial N[6]
and U[6]; the 121 ledger doubles in header declaration order; edges[n+1],
S[6*n], T[6*n]; and uint64 FNV-1a checksum of all preceding bytes. No struct
padding, host pointer or unaccepted candidate is serialized. Length is
`8*(143+13*n)` bytes. Cells are limited to 1..1,000,000 before allocation.

Unpack creates a new context, checks exact byte length/version/checksum/tag,
finite values, ordering, epoch/time constraints, ledger stoichiometry, energy
and inventory balances, and admissible initial grid moments. It does not
modify an existing context. The checksum detects accidental corruption, not
adversarial tampering. `model_tag` is a caller-generated configuration
fingerprint; callers must include grid, mass, reaction, collision, loss and
other selected model settings. It is neither automatically inferred nor a
security credential. Host checkpoints must carry the same epoch and time.

## Validation scope

C++ tests cover trial replacement, rejection, stale tickets, independent
contexts, a real FP step, new-mass DT nuclear bookkeeping, and malformed or
incompatible restart. The accompanying continuous source-state study uses
80 logarithmic cells, 400 steps of 1 ms, externally prescribed electron and
proton baths, constant external alpha injection and uniform 0.3/s escape.
It deliberately rejects a trial every 17 steps, replaces a candidate every
5 steps, and restarts at 0.2 s. Final packed states equal an uninterrupted
reference byte for byte. The independently known backward-Euler total-number
solution is checked at every step.

The study reaches 0.4 s with maximum normalized number/energy ledger residuals
2.97e-16/5.81e-15. The source is mapped to the nearest existing cell center,
3.011005 MeV; it is an explicitly prescribed external source, not a claim
about the physical pB birth spectrum. Baths are prescribed and do not receive
self-consistent feedback. This tests composition/transactions, not EXL
performance, grid convergence of a new spectrum, or completed host coupling.
Raw source and CSV/log are in `validation/source-state/`.

Address/undefined sanitizer checks pass with leak detection disabled because
LeakSanitizer cannot operate under this execution environment's ptrace. The
initial environment failure and successful ASan/UBSan log are both retained;
no leak-sanitizer pass is claimed. Initial energy anchors that cannot be
represented consistently in the restart format are rejected at creation.

Final validation: GNU, Intel default and Intel `-r8 -check all -traceback`
each pass 40/40 tests; no-Fortran build passes 25/25. Final high-bit-tag/layout
binding tests pass in all three builds, plus GNU default-real-8 and
default-integer-8 together. Installed C11 and Intel-r8 consumers link and run;
BALDUR's existing Intel-r8 executable links successfully without new-mode
activation. The installed example uses the exported `pb11Targets.cmake`;
a `find_package` configuration wrapper remains later packaging work.


## Inert-ion extension and restart schema 2

The additive `stage_inert` and `snapshot_inert` interfaces preserve collision
heat to non-network ion baths without changing `fusion_source_ledger_v1`.
An extra six-element signed finite array contains the sum over all inert
baths for EACH FAST species; it excludes the seven baths in the base ledger.
Both step and cumulative kinetic-energy validation include these amounts.
No extra ion density, collision closure or thermal advance is hidden in the
state manager. Host common-ion energy must include the same heat once.

The original stage interface means zero inert heat for that step, even when
continuing an already extended state. It preserves earlier cumulative inert
heat. The extended interface requires its array, including for all-zero heat.
Only successfully committing an extended stage promotes an old context.
Rejected/replaced/discarded trials do not promote the accepted state. Repeated
stage calls replace, rather than add to, both base and extra heat accounts.
Once promoted, old `snapshot` rejects rather than return incomplete accounts;
`snapshot_inert` reads either format and yields zero extra heat for v1 state.
The C extended snapshot clears scalar/ledger/extra-heat outputs on failure and
leaves kinetic output arrays untouched. Fortran retains its stronger existing
policy of clearing all outputs. Invalid Fortran extents invalidate earlier
staging through the C null-stage path before forming any array addresses.

Unpromoted contexts serialize exactly as schema1. Schema2 retains all existing
header fields and inserts six binary64 inert-heat amounts after the121-double
ledger, before grid edges. Its length is `8*(149+13*n)` bytes, exactly48 more
than v1. The version field is2; the common checksum covers the entire payload.
New unpack accepts both; old readers reject2. Schema2 with epoch0 is invalid,
because only a commit can promote state. Nonfinite inert heat and coherent-
checksum payloads that violate cumulative energy balance are rejected.

The standalone real-DT/carbon regression composes8 coupled steps over1ms,
with changing thermal temperatures. A split run packs the kinetic state and
captures its host thermal tuple at epoch4, restores both at that epoch, and
continues. Complete final thermal values and packed kinetic/ledger bytes equal
the uninterrupted result exactly. At this numerical stress point the carbon
receives about448.984J/m^3. Dropping its actual heat from the staged ledger
rejects the candidate. This is a local ownership regression, not validation
of BALDUR's separately owned common blocks or restart files.

An old-library/new-library binary fixture with six nonzero fast species and
signed-bath energy transfer confirms v1 byte compatibility. Independent tests
cover signed extra heat, lifecycle replacement/discard, semantic corruption
with a recomputed checksum, zero extra heat, and continued restart evolution.
See `validation/inert-source-state` for commands and retained evidence.
