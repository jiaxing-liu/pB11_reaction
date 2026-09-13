# Nuclear cross-section window evidence (source-only, 2026-09-14)

This note records published fit/data windows and the explicit limits on using
them. It does not choose a closure or extrapolation policy.

## Energy variable and windows

The external beam/thermal contract (`/home/cloud/research/pB11_reaction/docs/BEAM_AND_THERMAL_MOMENTS.md`)
defines the window in relative nuclear energy and says that the integrators
use only the fitted interval: they do not extend a fit, renormalize retained
probability, or treat missing data as zero cross section. A laboratory beam
energy is first transformed to relative CM energy. The host library implements
the same convention in `src/fusion_cross_sections.cpp:7,25-34` and
`src/pb11.cpp:13-17,55-58`.

| channel/model | published fit or retained data window | energy convention and evidence |
|---|---:|---|
| p–\(^{11}\)B, Tentori–Belloni (2023) | analytic piece: \(0\le E\le9.76\) MeV; source data SW: 0.14–3.48 MeV; BU: 5.70–9.76 MeV | \(E\) is CM. The paper gives incident-proton ranges \(E_p=0.15\)–3.80 MeV for SW (equivalent CM 0.14–3.48) and 6.00–18.0 MeV for BU, but retains the summed BU data only through CM \(E=9.76\) MeV (paper pp. 2–3, Eqs. (1)–(5)). |
| D(d,p)T, Bosch–Hale (1992) | 0.5–5000 keV | CM; Table IV pp. 620–621 and Table V p. 621 footnote. |
| D(d,n)\(^{3}\)He, Bosch–Hale | 0.5–4900 keV | CM; same tables. |
| T(d,n)\(^{4}\)He, Bosch–Hale | low 0.5–550 keV; high 550–4700 keV | CM; Table IV p. 621 and Table VI p. 622. The text says the two formulas meet at 530 keV, the best smooth transition; the C++ code consequently switches at 530 keV. This is distinct from the Table VI lower column value 550 keV. |
| \(^{3}\)He(d,p)\(^{4}\)He, Bosch–Hale | low 0.3–900 keV; high 900–4800 keV | CM; Table IV/VI. The text recommends switching at exactly 900 keV despite the small discontinuity. |

Bosch–Hale state on journal p. 612: “E denotes the energy available in the
CM frame,” with \(E_A=E(m_A+m_B)/m_B\) for a projectile A on stationary B.
Thus the listed Bosch–Hale domains are not laboratory projectile-energy
domains. Table VII is a separate Maxwellian reactivity fit: \(T_i=0.2\)–100
keV for DT/DD and 0.5–190 keV for D–\(^{3}\)He (journal p. 625); it is not a
cross-section-energy range.

## What the original papers say about extrapolation

Bosch–Hale Eq. (9) is a Padé fit to S-values from R-matrix cross sections. The
paper says a factorized form can be used mathematically at higher energy and
reports fitting calculated cross sections “up to 5 MeV” (journal p. 613), but
also states that Eq. (4) is strictly limited to incident S-wave dominance and
that the Gamow form is appropriate only well below the Coulomb barrier. The
same page reports that earlier fitted parameters gave “poor extrapolations of
the cross-sections to low energies.” Tables IV and VI explicitly label their
bottom rows as the validity ranges and maximum deviations (pp. 621–622).
There is therefore no source permission to evaluate a Bosch–Hale S(E) fit
below its tabulated lower endpoint (in particular, no published basis for
silently extending the DT 0.5-keV fit to E=0); any such continuation is a
separate sensitivity assumption.

Tentori–Belloni is more explicit about the pB boundary:

- The abstract says a reference cross section is proposed “up to 10 MeV”
  (journal p. 1), while the methods set \(E_3=9.76\) MeV and define the last
  piece as \(E_2<E\le E_3\) (p. 2). The conclusion repeats the rounded “up to
  10 MeV” wording but gives the actual retained BU interval 5.70–9.76 MeV
  (p. 8). The code's 0–9.76-MeV CM domain matches the explicit piece boundary;
  10 MeV must not be merged with 9.76 as if they were the same endpoint.
- The SW measurements begin at CM \(E\simeq0.14\) MeV, whereas the S1
  expression is written for \(E\le0.400\) MeV. The paper reports residuals
  for the fitted low piece but does not provide a low-energy uncertainty or
  an explicit validation statement below the first SW datum. The formula's
  mathematical E=0 evaluation is therefore an extrapolation gap, not a
  measured zero-energy domain. The 148-keV Breit–Wigner term is taken from
  Nevins–Swain because SW missed that narrow resonance (pp. 2–3).
- The paper states “no experimental data are available for 3.5<E<5.7 MeV”
  (p. 6) and says measurements connecting the two datasets are needed (p.
  8). Its S3 fit spans this gap, so values there are fit interpolation rather
  than direct measurements. It also reports a mean residual of about 30% for
  5.70–9.76 MeV and notes BU fluctuations (pp. 4, 8).
- For \(E>9.76\) MeV, the authors explicitly assume “the strongly decreasing
  behaviour” \(\sigma(E)\approx S(E_3)E^{-1}\exp[-\sqrt{E_G/E}]\) to estimate the
  thermal tail (p. 5). They conclude that exact knowledge above 3.5 MeV is
  not critical for their H–\(^{11}\)B ignition/burn temperatures, but in the
  same discussion stress that reliable cross sections well beyond 3 MeV are
  needed for nonthermal/beam fusion because yield depends on cross section and
  stopping power (p. 6). That thermal-tail estimate is not a validated
  high-energy beam cross section or a hard upper bound.

## Primary links and remaining gaps

- Bosch and Hale, *Nuclear Fusion* 32, 611–631 (1992), DOI
  [10.1088/0029-5515/32/4/I07](https://doi.org/10.1088/0029-5515/32/4/I07); audited
  local PDF: `/tmp/m2-physics-sources/bosch-hale-1992.pdf`.
- Tentori and Belloni, *Nuclear Fusion* 63, 086001 (2023), DOI
  [10.1088/1741-4326/acda4b](https://doi.org/10.1088/1741-4326/acda4b); accessible
  full-text mirror used for page-level checks:
  [ResearchGate article PDF](https://www.researchgate.net/publication/371625643_Revisiting_p-_11_B_fusion_cross_section_and_reactivity_and_their_analytic_approximations).

The source set supplies CM fit windows and qualitative tail statements, but no
general error bound for cross sections outside those windows. In particular,
the pB region below the SW floor (<0.14 MeV), interpolation gap
(3.48–5.70 MeV), and high-E (>9.76 MeV) require separate source accounting.
The companion PB_LOW_ENERGY_SOURCE_EVIDENCE.md establishes older Becker
coverage down to22keV; the SW floor is not the historical experimental floor.
A small thermal tail does not certify a beam/fast-particle rate.
