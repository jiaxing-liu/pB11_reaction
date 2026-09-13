# Conditional pB CM birth composition

This implementation adds the missing alpha0 marginal and composes it with
existing alpha1 sources. It does not infer incident-resonance weights,
calibrate channel branches, or deliver a laboratory source. The full
[physical source-family design](PB_SOURCE_FAMILY_DESIGN.md) remains binding.

## Alpha0 formula

Use canonical alpha mass from `fusion_c_nuclear_mass(FUSION_HELIUM4)`, exact
SI c and the existing exact sequential kinematic API. For parent rest energy
`M c² = 3 m_alpha c² + A` and intermediate `B c² = 2 m_alpha c² + q`,
alpha0 is a delta at primary kinetic energy K0. The spin-zero intermediate
has isotropic secondary direction in its own rest frame, giving each of its
two alphas a uniform energy marginal on [Kminus,Kplus] in the parent CM.

The box endpoints come from the complete on-shell event with cosine1; its
energies are recovered from momenta without subtracting nearly equal rest
energies. In particular stationary-secondary and zero-width limits remain
nonnegative. The identity is `K0 + Kminus + Kplus = A`, with3 alphas/event.
The two box particles are correlated: this marginal does not constitute
three independent laboratory event samples.

Ground-state q=91.84(4)keV and natural width5.57(25)eV are the previously
audited Refsgaard Table2 inputs; see PB_SOURCE_FAMILY_DESIGN.md. The API
requires explicit q, recording the narrow-width approximation rather than
hiding another nuclear Q convention. Other q values are useful mathematical
or sensitivity inputs, not measured Be8 ground-state energies. The total
available energy for physical events must use canonical nuclear Q plus the
reactant CM kinetic energy. The API's0<A<=12MeV is a numerical domain.

The arithmetic-center grid projection integrates both linear basis functions
analytically over every box/center-interval overlap. A zero-width box is two
delta particles. Below/above-center-hull number and energy remain separate
ledgers: no clipping, implicit thermalization or physical escape is assumed.
A one-cell grid has only a point center, so a continuous box is entirely in
its below/above ledgers; this is intentional and conservative.

## Composition contract

`fusion_c_pb_alpha0_grid` supplies that per-event marginal.
`fusion_c_pb_cm_source_grid` takes total-event fractions f0 and flow, using
fbroad=1-f0-flow. The three components are alpha0, low-parent primary-l2
alpha1, and broad-region primary-l1/l3/mixture alpha1. This is an incoherent
mixture of declared event populations; interference inside each alpha1
amplitude is retained. It is not an incoherent replacement for interference
between overlapping entrance amplitudes of the same quantum numbers.

Both fractions must be nonnegative and sum<=1. All source controls are
validated even when a branch has zero weight; expensive alpha1 evaluation
is skipped for an inactive branch. The caller explicitly chooses FSCI0/1,
broad mode1/3/13, unit-basis coefficient k and phase. Alpha1 keeps its
published nonrelativistic R-matrix convention; alpha0 uses exact kinematics.
Global N/E closure does not certify equal angular or relativistic accuracy
between those models. Laboratory composition must address that distinction.

The result is14 C doubles: mapped/below/above N/U, two residuals, three
fractions and three alpha0 endpoint diagnostics. The Fortran module uses
explicit ISO kinds and checks actual array extents before calling C.
Inputs and outputs must not overlap. No hidden state or host indexing is
introduced, and old source interfaces remain unchanged.

## Independent reference and conditional illustration

A120-case mpmath80-digit reference uses the two-body invariant Kallen
formula and direct Lorentz energy endpoints, independently of the production
momentum-energy recovery. Separate antiderivatives integrate the grid hats.
Cases cover physical q, q=0/A, near-degenerate q, a stationary-secondary
neighborhood, one-cell/nonuniform grids and partial/complete spill.
Maximum absolute per-bin/spill count error6.4143e-15; maximum spill energy
error/A2.6600e-16; maximum endpoint error/A6.5052e-17. This is numerical
verification of the specified conditional model, not experimental fitting.

The illustration at A=8.84MeV uses the already refined1024-node l2 FSCI
source and alpha0 fractions0,5 percent,100 percent. Five percent is an
illustrative event fraction, not a calibration to the low-resonance measured
branch (whose ghost/selection convention remains relevant). The100-bin grid
has60keV spacing over0–6MeV. A plotted delta peak's height depends on grid
width; its area and energy are the physically meaningful quantities.

Remaining work: incident-energy/branch selection with source uncertainties;
reactant-selected laboratory motion/angular probabilities; composition with
burn, collisions and handoff; self-consistent baths and host atomic
acceptance; EXL four-fuel full-window validation. Numerical source closure
alone does not close I007/I009/I010 or complete M3–M7.

## Final build and installed interface checks

GNU, Intel default and Intel-r8 each pass44/44 tests. The C++-only build
passes27/27. Installed C11 and Intel-r8 consumers pass, and BALDUR Intel-r8
builds. The dedicated C++ test also checks8192-point midpoint packet mapping,
mixture corner parity, interior weighted sums and error-output clearing.

Initial matrix failures came from the test's keV-to-MeV multiplication
rounding its cutoff one ULP below the1keV endpoint. The test now uses the
same canonical MeV literal; no production domain was widened. The first
installed consumer caught an omitted new header in the install list; this
is fixed. These failures are retained beside the final passing logs. No
new source ownership or host time evolution is activated by these builds.
