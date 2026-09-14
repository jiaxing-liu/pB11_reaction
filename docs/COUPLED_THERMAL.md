# Local thermal burn and fast-particle feedback trial

`fusion_c_coupled_thermal_trial` composes existing physical operators in the
independent library. It advances a trial of thermal fuel densities, electron
and common-ion energy, and the six fast-species S/T energy distributions.
The C and Fortran interfaces contain no BALDUR indexing, geometry or I/O.
The caller accepts or discards all outputs together; repeated evaluations
from the same accepted inputs do not accumulate events or heat.

This is a local first-order split operator. It is not yet a complete BALDUR
adapter or an implicit iteration of all bath and nuclear coefficients.
Thermal ions share one temperature, matching BALDUR's single Ti variable;
electrons have their own temperature. Unequal reactant-ion temperatures and
fast-fuel beam-target nuclear reactions remain separate required extensions.

## Physical and numerical sequence

1. Derive old Ti from total thermal-ion energy and the number of network plus
   explicit inert ions. Evaluate the complete common-temperature lab source
   for enabled channels with available fuel. Freeze its rate, selected-reactant
   debit means and energy distribution for this trial. `birth.pb_low` applies
   to pB only; the other channels use their own native low-energy model.
2. Solve simultaneous backward-Euler thermal fuel competition with those
   frozen coefficients. Both DD branches include their symmetry factor once.
   Products are born fast; T/He3 transferred to thermal fuel at the end of
   this step can burn only in a later step. Thermal fuel removal uses the
   reaction-conditioned energy debit, not a blanket 3Ti/2 per reaction.
3. Convert each channel event amount to its full product distribution using
   the SAME finite source rate. Charged grid spill rejects the entire trial;
   it is not clipped, called ash or deposited as heat. Neutrons retain their
   complete number/energy including source-grid spills, in a separate output.
4. Recompute the depleted common Ti, freeze the Coulomb baths, and apply the
   existing conservative two-component implicit FP trial to each fast species.
   External grid births enter S; explicit physical escape carries its own
   particle and energy account. The S-to-T transfer is internal kinetic
   bookkeeping, not fluid helium production. Boundary grid flux otherwise
   reflects, as in the existing FP operator.
5. Add all signed electron/ion collision heat to their actual thermal energy
   pools. The next accepted step derives its temperatures from the updated
   energies, so the background evolves in response to fusion/fast-particle
   transfer; it is not a prescribed bath history. Negative thermal energy
   rejects the trial; there is no temperature floor or hidden heat source.
6. Optionally apply measured whole-T-candidate handoff. For one candidate
   Nf,Uf and current common ion pool Ni,Ui, use the mixed target
   `Ti_mix=2*(Ui+Uf)/(3*(Ni+Nf))`. The shape and mean-energy gates apply at
   this temperature. If accepted, fluid energy plus signed bath correction
   equal Uf, and the resulting pool temperature is exactly the target. Species
   are processed in ID order 0..5, recomputing the pool after each acceptance.
   This declared order and first-order split require convergence/sensitivity
   checks; a rejected handoff leaves the complete kinetic candidate intact.

The old thermal-burn operator checks positive species energy while depleting
fuel. This can impose a stronger timestep restriction than the common-ion
pool alone requires; this implementation retains that conservative rejection
instead of redistributing an unavailable species debit silently.

## Inert ions, charge and source ownership

Carbon and other non-network ion baths are explicit inputs, with fixed density,
mass and mean squared charge. Their thermal energy belongs to the common Ui
pool. Their collision heat is returned separately, per fast species, in
`inert_ion_heat_J_m3[6]`. It is never mislabeled as hydrogen or electron heat.
Network bath heat retains the existing 6-by-7 layout (electron then ions0..5).
For charge mixtures, the collision input is `<Z²>`, not `<Z>²`.

The embedded source-state-v1 ledger alone does not include the extra inert-ion
heat. Pass it together with `inert_ion_heat_J_m3` to the additive source-state
`stage_inert` interface. Its version2 restart preserves both accounts; use
`snapshot_inert` to recover them. [SOURCE_STATE.md](SOURCE_STATE.md) defines
promotion, compatibility and same-epoch thermal ownership. The host still
must preserve its thermal and geometry state atomically with this context;
stateless operator success does not implement BALDUR's acceptance path.

The electron density is an explicit frozen trial input; only Ue is advanced
here. The host owns charge balance, ionization, ambipolar boundary electron
losses and their energy, neutral/beam injection accounting, radiation and
spatial transport. A fixed-ne standalone test with fully stripped thermal
fuel and no charged injection/escape is consistent because the nuclear network
conserves charge including fast products. The function does not invent an
ionization-energy model for partially stripped boron or other ions.

Fast particles use canonical bare nuclear charges/masses. Products are born
with the prior exact on-shell kinematics; the subsequent collision operator
remains classical, NR and dilute. Energy conservation of this combination is
not a claim of relativistic slowing-down accuracy. Pair Coulomb logarithms
are supplied explicitly and frozen for the step; no screening closure is
silently selected by this function.

## Accounting and failure semantics

All quantities returned in the ledger are amounts per step. For the closed
thermal/fast system with separate neutron and charged escape outputs,

`Delta(Ue+Ui+Ufast) = Q_events + U_external - U_neutron - U_escape`.

The implementation checks nuclear energy, each fast-species number/energy,
each complete thermal-plus-fast species balance, and the total local budget.
Physical collision heat, handoff energy/correction, neutron energy and
external sources enter exactly once. The extra inert-ion heat participates
in each fast-species energy balance and the ion pool update.

Relative ledger tests use actual terms. Complete state subtraction
also allows 16 machine epsilons times the actual stored thermal and kinetic
inventory magnitudes:
a source smaller than one ULP of a large bath cannot change that bath's stored
double. There is no arbitrary unit-sized absolute tolerance. Returned signed
residuals remain available for accumulation and case-level budgets.

Both numerical rate/debit discrepancy tolerances are explicit. Passing them
is necessary, but does not certify pB event-angle/source-proxy convergence.
The header freezes domains and array layouts. No input/output overlap is
permitted. Invalid calls clear output arrays for valid dimensions and never
modify accepted inputs. No state, file I/O or global cache is hidden here.

## Standalone evolving-bath evidence

`study_coupled_thermal` is built with the test-enabled CMake configuration.
Arguments are steps, cells, fuel (`dt/dd/dhe3/pb`), duration seconds,
handoff-enabled (0/1), and pB nq=ncos (default8). It performs actual local
trial/accept steps and accumulates neutron/Q/energy accounts independently.
The density is1e22m^-3 and a fully stripped carbon bath has nC=0.001ne;
this is a numerical coupling stress test, NOT an EXL device prediction or
a fair four-fuel performance ranking. DT/DD/DHe3 start Ti20keV and Te5keV;
pB uses Ti100keV. Coulomb logs are prescribed15 for the declared classical
model. No charged external source/escape is enabled in these closed tests.

The DT10ms study heats electrons from5keV to about20keV and ions from20keV
to about22.3keV. Complete accumulated energy residuals are around1e-16 in the
recorded runs. At800 cells,64→128 timesteps changes final Ue by0.187%, Ui by
0.0011%, fast energy by0.568%, and Q by0.0545%. At128 steps,800→1600 cells
changes Ue by0.0848%, Ui by0.0153%, fast energy by0.191%, Q by0.0121%.
These are measured numerical changes for this DT point, not universal bounds.
The strict handoff gates produce no projections in these continuous hot-source
runs; a thermal-scale T component is not automatically a thermal fluid.

A first pB run was correctly rejected at step2 because the lab source extended
below its first cell center. The preserved failure uses first edge1e-10keV;
extending it to1e-18keV, without clipping, allowed the4-step1ms study to finish.
That grid choice is not a universal guarantee of source containment. Dedicated
pB timestep, grid and angular resolution studies remain distinct from DT
convergence and from single-step accounting tests.

Full compiler/binding tests, independent invariant checks, installed consumers
and immutable study outputs accompany the published milestone. Atomic host thermal/geometry
state ownership, remaining nuclear-source closures and host acceptance are still
required for the original program. Verified common-T source interpolation is now
provided by the additive variant described below, for explicit validated ranges.


The archived pB1ms refinement (`validation/coupled-thermal/pb-refinement`)
compares8/16/32/64 steps at800 cells, with800→1600 cells and nq=ncos8→16
at8 steps. From32→64 steps, fast energy changes0.515%, cumulative electron
heat0.598%, network-ion heat1.891%, and carbon heat2.047% (0.49631J/m^3).
The corresponding total Ue change is only0.00610%; this much smaller number
must not be used to hide weak source sensitivity. Carbon's relative change
is still above2%, and all values are measured differences, not rigorous error
bounds. Independent angular/proxy source accuracy remains open.

The original continuous studies linked the pre-inventory-roundoff-guard
library. The successful evolution equations are identical after the guard
correction; an archived DT32-step rerun compares all common numerical columns.
Source hashes captured by the plotting script identify analysis-time source
snapshots, not necessarily the earlier executable's build inputs. The pB
refinement manifest preserves the launch-time binary/library hashes separately.


## Explicit table-backed variant

`fusion_c_coupled_thermal_table_trial` shares this physical operator with the
original direct entry point, but obtains complete birth coefficients from
immutable, validated tables. It checks channel, model options and exact grid;
there is no silent extrapolation or fallback. See [THERMAL_BIRTH_TABLE.md](THERMAL_BIRTH_TABLE.md)
for sampled error semantics, evolving direct/table comparisons and limitations.
The original entry point still evaluates the direct quadrature each trial.


## Explicit fluid increments for host source assembly

`fusion_c_coupled_thermal_increment` maps a validated complete source ledger
plus six inert-ion heat amounts to signed thermal particle and electron/ion
energy increments. The64-byte `fusion_thermal_increment_v1` has six number
amounts and two energy amounts; the existing result layout remains unchanged.
The Fortran wrapper uses explicit C kinds and validates inert-array extent.

For each thermal species dN=handoff-consumed. Electron dU sums the electron
heat column. Ion dU sums all network-ion and inert heat, adds handoff fluid
energy and subtracts consumed thermal energy. The mixed-pool handoff correction
is already in the corresponding heat column and is counted exactly once.
External birth and escape ledger fields describe fast particles and must not
be added directly to the thermal fluid. Units are per step m^-3 and J/m^3.

The helper uses explicit amounts, preserving weak sources that disappear in
final-minus-initial background subtraction. It checks finite fields and the
ledger's sign conventions, and rejects nonrepresentable output; it does not
validate full reaction stoichiometry/kinetic closure or advance any state.
Use a ledger from a successful validated operator. It does not supply host
time/volume conversion, particle transport, charge closure or acceptance.
Evidence: `validation/thermal-increment`.
