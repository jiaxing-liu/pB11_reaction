# Conservative thermal birth tables and coupled evaluation

`fusion_birth_table.h` adds an immutable, host-independent table for one
reaction channel's common-temperature birth coefficients. It accelerates the
existing direct quadrature; it does not replace the nuclear model, derive an
unequal-temperature/beam source, or advance plasma state. The direct API remains
available with its previous semantics.

## Construction and interpolation contract

Inputs explicitly specify the channel, temperature interval, energy edges,
complete source-model options, interpolation tolerances and resource limits.
The constructor copies configuration and grid. It evaluates direct endpoint
sources and geometric quarter/mid/three-quarter temperatures in each candidate
interval. Failed interpolation intervals bisect in log temperature; resource
or direct-source gate failure returns no usable context. The info query returns
configuration, accepted-knot and direct-evaluation counts, sampled interpolation
errors and actual sampled direct rate/debit discrepancies.

One convex weight interpolates the RAW rate coefficient K, both conditional
reactant-energy numerator moments Da/Db, every species/energy-cell coefficient,
and all below/above-grid number and energy moments. It does not interpolate
normalized spectra separately and does not independently interpolate logarithms
of different coefficients. Thus the linear identities are retained:

- Product number coefficient for each species is its multiplicity times K.
- Sum of product energies, including neutrons and spills, is Q*K + Da + Db.

The implementation checks these identities against canonical nuclear data for
direct nodes and every returned interpolation. All coefficients are nonnegative.
When the burn solver needs conditional mean reactant energies it divides the
interpolated numerator by the SAME interpolated K. Density products and the
identical-reactant symmetry factor belong to the burn/rate caller, exactly once.

## What the error checks measure

Rate and each debit numerator are checked relative to their direct references.
Number L1 includes grid-bin and spill-number differences. Energy L1 includes
arithmetic-center-weighted grid differences and spill-energy moment differences;
it does not reconstruct spectra outside the grid. Both L1 tests normalize each
product species separately and use the maximum over all seven species, so a
neutron-dominated energy budget cannot hide an alpha spectrum discrepancy.
A zero reference requires exactly zero difference. No unit-sized floor or
additional weak-source truncation is inserted.

These are finite sampled checks, not mathematical supremum bounds. Independent
non-construction temperatures and table refinement remain necessary; numerical
agreement with the direct evaluator does not establish its nuclear/source-model
accuracy. In particular pB angular quadrature and the named remainder-spectrum
proxies retain their own validation/sensitivity requirements. Direct finite-tail
and quadrature gates are distinct from interpolation gates. A truncated direct
source cannot be repaired by adding more temperature nodes.

Temperature queries outside the constructed interval reject, without clamping
or hidden extrapolation/direct fallback. Model-domain, negative/nonfinite input,
resource and conservation failures remain explicit. Exact knots return the
stored direct coefficients. The info/control/result layouts and all units and
limits are frozen in the header. Fortran uses ISO_C_BINDING and `(cells,7)` grid
arrays. A handle has one caller-owned lifetime; immutable evaluation may be
shared by independent zones, with no mutable globals or implicit host state.

## Use in the local coupled operator

The additive `fusion_c_coupled_thermal_table_trial` accepts five table handles.
Every enabled channel with available thermal reactants needs a table matching
its channel, complete source options and exact energy grid. Unused entries may
be null. The existing non-pB normalization of the pB low-energy selector is
retained. The table-backed and direct entry points share the SAME burn,
collision, handoff, escape and thermal-feedback implementation.

Table-mode source discrepancy outputs are maxima sampled while constructing the
table, not direct evaluations at the queried temperature. Inspect table info for
its separate interpolation gate settings and observed maxima. The caller owns
acceptable interpolation accuracy and the thermal interval. No table validity
statement promises an entire evolving trajectory will remain within it.
For example DT initially cools from20keV to about19.9981keV in the endpoint
regression because the reaction-conditioned fuel debit precedes delayed product
heating. A table starting exactly at20keV correctly rejects the next step;
explicitly constructing a19.5–23keV table covers the tested trajectory.

Tables are configuration, not accepted kinetic inventory. Source-state v2 stores
complete accepted ledgers including inert-ion heat; host checkpoints must pair
thermal/geometry state and preserve all source/table configuration in their
model fingerprint. BALDUR activation, source ownership, transport and restart
I/O still require a thin host adapter.

## Evidence and cost at the measured points

At400 energy cells, the DT19.5–23keV table with0.1% interpolation gates uses
5knots/23direct evaluations, built in about0.46s in the recorded microbenchmark.
Four independent temperatures pass; lookup is about11microseconds versus
20milliseconds direct. The pB100–101keV table uses2knots/5direct evaluations,
about40s construction; lookup is about11microseconds versus8s direct. These
are measured CPU costs in this environment and these ranges, not universal
speedups or construction-cost guarantees. Source nq=ncos8 is a declared model
resolution; interpolation accuracy does not certify that resolution.

Actual evolving-bath comparisons use800cells/64steps, DT10ms and pB1ms, with
otherwise identical direct/table inputs. DT table19.5–23keV changes final fast
energy by0.00293% and carbon heat by0.00396%. pB table99–101keV changes fast
energy by0.0311%, electron heat0.0307%, network-ion heat0.0294%, carbon heat
0.0293%, and Q0.0298%. Total local energy residuals remain around1e-16.
These are local high-density numerical stress tests, not EXL predictions or
fair four-fuel comparisons. The raw files and figures are in
`validation/birth-table/studies`; complete compiler and installed-consumer
checks accompany the milestone. Direct-mode trajectory regression, exact-knot
physical parity, mismatched model/grid rejection, and out-of-range failure are
covered separately from these finite-time sensitivity measurements.

DD and D–He3 companion comparisons (800 cells, 16 steps, 10 ms) additionally
exercise the other channels. Final fast-energy differences are0.01159% and
0.01744%; all three separately recorded bath heats differ by less than0.019%
through these histories. Both DD primary branches are exercised. New/previous
direct histories agree exactly in all shared columns. Raw data, exact commands,
absolute and relative errors are in `validation/birth-table/companion-fuels`.

Complete GNU, Intel default and Intel-r8 suites pass58/58 each; no-Fortran
passes36/36. Installed C11 and Intel-r8 Fortran consumers construct actual
non-null DT tables and call the coupled operator successfully. Host Intel
build/link passes without activating the new source mode. Full evidence is in
`validation/birth-table/build-checks`.
