# Explicit final-state Coulomb correction for conditional alpha1 sources

The new policy implements the phenomenological Model-II correction in
Refsgaard et al., *Physics Letters B*779 (2018)414–419, printed p415 Eq4.
See [the equation audit](PB_FINAL_STATE_COULOMB_EVIDENCE.md) for the primary
sources and the distinction from Laursen2016 Eq5. The later paper identifies
a calculation error behind the earlier reported10fm separation. It discusses
approximately15fm for short-lived excited-Be8 breakup and uses16fm in its
Model-II table. Here **16fm is an explicit model choice**, not a measured
universal distance or a new pB channel measurement.

## Formula and interfaces

For each coherent permutation, replace

`P_l1/rho1 -> (P_l1/rho1) C`,

`C = (rho1_tilde/P_l1_tilde) (P_2_tilde(E12)/rho12_tilde) (P_2_tilde(E13)/rho13_tilde)`.

The tilde functions use16fm; ordinary5.1/4.5fm radii, phases, level parameters
and CG algebra remain unchanged. The **amplitude** acquires sqrt(C), and the
three corrected amplitudes are then summed coherently. Logarithms combine
ordinary penetration and C before exponentiation; an underflowed ordinary
amplitude is not multiplied by a separately overflowing correction.

- `fusion_c_nuclear_coulomb_radius16`: independently generated16fm tables,
  same channel IDs and .001–12MeV numerical domain as the ordinary function.
- `fusion_c_alpha_amplitudes_fsci_cutoff`: policy0 NONE or policy1
  REFS2018_R16, with the explicit1–10keV numerical cutoff also applied to
  the two additional pair energies. Every suppressed permutation is counted;
  all three terms may be suppressed. Other errors are returned.
- `fusion_c_alpha_spectrum_model_grid`: normalized CM source with explicit
  policy. The old entry point permanently selects NONE. Mixtures normalize
  both bases **with the selected correction** before interference.
- `fusion_fsci_fortran` supplies the three ISO_C_BINDING wrappers, reusing
  the existing spectrum result layouts and checking extents before C_LOC.

Quadrature now permits4–1024 nodes per variable. The expanded bound was
exercised in the complete source study, and1025 is rejected. This is a
numerical control, not a guarantee of uniform pointwise shape accuracy.

## Independent numerical checks

The new-radius generator uses mpmath40-digit F/G functions, a Wronskian
check,17-point Chebyshev segments and seven independent validation abscissae
per segment. The tables contain34 segments for each alpha-Be8 channel and19
for each alpha-alpha channel. All stored thresholds pass. Thirty-five
independent reference points also test the compiled C evaluator, including
both energy endpoints. Six50-digit direct F/G references test the actual
unsymmetrized amplitude correction, including a factor106.7394 case.

The generator's cache fingerprint now includes charge product, reduced mass,
radius and partial waves. Alternate radii require separate output paths and
namespace. Spawned workers initialize the same radius explicitly. The first
attempt encountered this environment's forkserver socket restriction; the
generator now selects spawn directly. Cached regeneration reproduces the
versioned include byte for byte; no default-radius table was regenerated.

## Source-shape study

The study uses100 bins of60keV,0–6MeV, and retains separate below/above-center
number and energy. A is a declared conditional parent-CM energy,8.84MeV for
primary l2 and9.30MeV for the coherent l1/l3 example. It is not inferred from
a nominal excitation label or represented as a new experimental fit.

Distance is `sum(abs(delta N_bins), including spill categories)/(2*3)`.

| Conditional model | 256→512 nodes | 512→1024 nodes | NONE→FSCI at1024 |
|---|---:|---:|---:|
| l2, A8.84MeV | .00429/.00427 (NONE/FSCI) | .00199/.00200 | .0152854 |
| l1/l3, A9.30MeV, same unit-basis k=.76 | .00420 (FSCI) | .00157 | .00501409 |

Individual bin-shape refinement is still0.16–0.20% in this integrated metric;
this is not a pointwise confidence bound. The paired **policy difference**
changes much less on512→1024 refinement:3.85e-5 for l2 and1.32e-5 for the
mixture. Thus the displayed policy effects are distinguishable from the
observed quadrature changes. Total N=3 and total kinetic energy=A, including
spill, close to rounding. At1024 nodes, changing the cutoff1→10keV changes
the corrected shapes by less than1e-8 in this metric.

Holding unit-basis k fixed while changing the correction does not hold the
raw amplitude coefficients fixed. If N1,N3 are the uncorrected norms and
N1c,N3c the corrected norms, the latter comparison uses

`k_c = (k N1c/N1) / [k N1c/N1 + (1-k) N3c/N3]`.

Using the512-point norms gives k_c=.75863357611887294 for original k=.76.
At1024 nodes the fixed-raw-coefficient policy distance is.00482219. This
transformation does not resolve the previously documented experimental
mixture-convention/radial-model uncertainty. Raw tables, code, norms and both
comparisons are in `validation/fsci-source/`.

## Scope and next composition

This adds a documented alpha1-source model option. It does not yet select
incident-resonance weights or alpha0 branching, supply a laboratory source,
evolve thermal backgrounds or activate new BALDUR source ownership. The
[full source-family design](PB_SOURCE_FAMILY_DESIGN.md) records these remaining
steps and the analytic alpha0 marginal. Existing I009/I010 evidence is not
silently declared resolved by conservation or successful API tests.

## Final build and interface matrix

GNU, Intel default and Intel `-r8` each pass42/42 tests; the C++-only
build passes26/26. Fresh installed C11 and Intel-r8 consumers compile and
run successfully. BALDUR Intel-r8 builds; new source ownership is not yet
activated. The original Fortran boundary test expected513 to be rejected;
after expanding the documented limit to1024 it was updated to test1025.
Both initial failed matrix logs are retained alongside the final passes.
