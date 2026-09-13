# Thermal nuclear-window contribution study

This study quantifies candidate continuation contributions before selecting a
production closure. It is NOT an assertion that a small probability tail is
a cross-section-error bound. No production cross-section or rate API changes
in this increment. Root owns the formulas and numerical study; the source-only
window audit is in NUCLEAR_WINDOW_SOURCE_EVIDENCE.md.

## Explicit assumptions and independent integration

For each of the five channels, use the existing cross section only inside its
published fit interval. Outside it, study the following additional assumptions:

1. Constant endpoint S-factor with the same fixed Gamow exponent B_G:
   sigma(E)=sigma(Eb)*(Eb/E)*exp[B_G*(1/sqrt(Eb)-1/sqrt(E))].
   Here E and Eb are in keV, B_G in sqrt(keV). For BH use the positive low
   endpoint below the fit, and the upper endpoint above it. For pB the analytic
   piece already reaches0, so no extra low-end segment is added.
2. Above the upper endpoint, also calculate constant sigma(E)=sigma(Emax).
   This is an alternative sensitivity model, NOT an experimental upper bound.

The pB constant-S high tail matches the thermal-tail assumption discussed by
Tentori–Belloni2023 p5. This does not validate it for energetic beams. Applying
these continuations to BH is a separate modeling assumption: the BH paper
specifies finite fit domains and supplies no blanket permission to extrapolate
its rational S polynomial to zero. In particular, the DT polynomial's E->0
constant is very different from S at its0.5keV fit boundary. The study does
not select that unconstrained rational extrapolation.

For x=E/T, the independent energy-space quadrature is

```
K = sqrt(8*T_J/(pi*mu))*barn * integral sigma_barn(T*x)*x*exp(-x) dx
K*<Erel> = same_prefactor*T_J * integral sigma_barn(T*x)*x^2*exp(-x) dx.
```

Use the coherent new-model masses. Integrate fit, low continuation and high
continuation separately, with no probability renormalization. The numerical
integration stops at x=800 where the exponentially decreasing integrands are
below double working precision for these models; this is not a finite-support
plasma assumption. Explicit knots resolve Boltzmann scales, pB resonances and
all piecewise fit transitions. The fit-only result is independently compared
with the existing relative-speed-space thermal-pair integral.

The65 cases cover all five channels at common temperatures.01,.03,.05,.1,.2,
.5,1,3,10,30,100,190,500keV. Repeating with quadrature tolerance1e-10 then1e-12
changes every reported integral/mean by less than1.5e-14 relative. Independent
fit-window K and relative-energy moments agree within1.22e-14 and1.17e-14.
This validates numerical evaluation of these declared assumptions, not their
unmeasured cross-section uncertainty.

## Findings that constrain the next model decision

- With constant-S low continuation, at T>=.1keV on the sampled grid the
  DD low-end event-rate fraction is below4.9e-6, DT below7.5e-7 and DHe3
  smaller. At.03keV DD's low-end fraction is about12%, so a broad statement
  that BH windows are always complete would be false. Those very cold
  absolute rates are extremely small; full-case acceptance must distinguish
  relative reaction uncertainty from absolute particle/energy-budget effects.
- At500keV, the high-end rate contribution remains below5.2e-4 for the studied
  constant-S/constant-sigma alternatives across the four BH channels. Their
  relative-energy moment fractions are larger than event-rate fractions;
  rate-only checks cannot certify product energy moments. The constant-S
  high-tail relative-energy fraction reaches about0.246% for DD at500keV.
- For pB at500keV, the high constant-S tail contributes about1.78e-9 of the
  rate. This small high tail does NOT address pB low-energy model dependence.
- At pB temperatures1–3keV, essentially all calculated fit-window reactivity
  comes from E<140keV, below the lower SW dataset point used by TB2023. At
  10keV the fraction is about83.9%; at100keV it is about0.348%. This is a
  dataset-coverage diagnostic. Becker1987 reports measurements from22keV
  to1100keV; NS2000 uses Becker low-energy data. See
  PB_LOW_ENERGY_SOURCE_EVIDENCE.md for the independently checked source chain.
- The pB3.48–5.70MeV interpolation gap contributes about0.398% at500keV and
  becomes negligible at lower temperatures. Its role for nonthermal fuel
  remains separate from this common-Maxwellian study.

The leading next work is therefore a declared
continuation/acceptance policy with per-channel rate AND energy diagnostics.
Do not introduce an undocumented lower-temperature cut, set missing sigma to
zero, tune temperatures to hide a data gap, or reuse these thermal fractions
as bounds for arbitrary beam distributions. I008 remains open until the
composed production policy is implemented and its case-level consequences
are checked. This is not an EXL run or completed physical pB source.

## Reproduce

```
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target study_nuclear_windows -j4
build/study_nuclear_windows 1e-10 > nuclear-window-study.csv
build/study_nuclear_windows 1e-12 > nuclear-window-study-refined.csv
MPLCONFIGDIR=/tmp/pb-mpl python tools/plot_nuclear_windows.py \
  nuclear-window-study-refined.csv nuclear-windows
```

Raw CSVs and convergence JSON: docs/validation/nuclear-windows/. The plot's
left/middle panels show modeled continuation fractions, while its right panel
shows pB dataset coverage. Neither is a confidence band. The reactor results
must separately state the selected cross-section model and its limitations.

## Historical low-energy coverage and parameter sensitivity

Becker's reported experimental floor is22keV, not SW's140keV. The model's
rate fraction below22keV is82.97% at1keV temperature,0.8666% at3keV,
and1.4305e-6 at10keV. These are energy-coverage fractions, not uncertainty
bounds; they do not establish pointwise accuracy above22keV.

For E<=400keV, independently evaluate the NS low polynomial as transcribed
in TB2023 Table1 (C1=.240,C2=2.31e-4), keeping its inherited148keV narrow
term, and compare with TB's C1=.269,C2=2.54e-4. Above400keV keep TB
unchanged. This hybrid is explicitly a LOW-PIECE sensitivity, not the full
NS reactivity model. Also vary only TB C0 by+/-12MeV barn, motivated by
Becker's reported fitted S(0) uncertainty. The narrow-resonance tail means
C0 is not exactly the complete S1(0); this perturbation is not a propagated
experimental confidence band or a replacement for missing covariance.

| Common T (keV) | Below22keV / Kfit | NS low-piece change / Kfit | C0+12 change / Kfit |
|---:|---:|---:|---:|
|1|0.829696|-0.00270160|0.0590561|
|3|0.00866614|-0.00565202|0.0572409|
|10|1.43047e-6|-0.0106324|0.0448385|
|100|1.85283e-12|-0.0137014|0.0148922|
|500|2.45749e-14|-0.00219138|0.00214725|

C0-12 gives the opposite change because this contribution is linear. The
study checks that the mean of the plus/minus integrals recovers the central
TB low-piece integral. All added columns satisfy the same two-tolerance
convergence check (maximum relative change across all columns1.49e-14).
Numerical convergence does not eliminate the several-percent low-energy
parameter sensitivity. NS original full text and covariance were not obtained;
its parameters here are explicitly TB's transcription.
