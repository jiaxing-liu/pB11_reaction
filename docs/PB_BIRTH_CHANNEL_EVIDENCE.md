# p+11B birth-channel evidence audit

Scope: primary-source evidence for the sequential `11B(p,3α)` labels
`α0 = 8Be(gs)` and `α1 = 8Be(first excited 2+)`, at the low 16.11-MeV
resonance and the broad 16.6-MeV resonance.  The `l` values below refer to
the decay step `12C* → α_primary + 8Be*`; they are not entrance-channel
`p+11B` partial waves.  No off-resonance branch fraction is inferred.

## Energy convention

`Ep` in the cited experiments is the kinetic energy of a proton in the lab
incident on a stationary `11B` target.  For this two-body entrance channel,

```
Ecm = m(11B)/(m(p)+m(11B)) Ep_lab ≈ (11/12) Ep_lab.
```

Laursen et al. state this convention explicitly in their Eq. (1):
`E = (3/2)E1 = (11/12)Ebeam + Q − E'`, with `Ebeam` the proton lab kinetic
energy and `Q = 8.682 MeV` for `11B(p,3α)` (printed p. 3).

The nominal low resonance is often called `Ecm ≈ 148 keV`; that corresponds
to `Ep_lab ≈ 161.5 keV` under the approximation above.  Munch & Fynbo fit
the thick-target scan with `Er = 162.6(5) keV lab` (printed p. 4, Eq. (4)),
or `Ecm ≈ 149.1(5) keV`.  Laursen's Table 2 uses 167–170 keV lab
(`Ecm ≈ 153–156 keV`); those are beam settings for a finite target and
should not be treated as a monoenergetic 148-keV-cm point.

Munch & Fynbo 2018 accelerated `H3+`; their stated 525-keV molecular beam
setting is `525/3 = 175 keV` per proton for the reaction-energy convention
(printed p. 3), and their resonance scan fit is quoted in proton lab energy.

For the broad state, Kuhlwein et al. title and text use `16.62 MeV`,
`2−, T=1`, and say it is 663 keV above the proton threshold (printed p. 1).
Their Fig. 1 labels the level `16.57` and `Ep = 675 keV` (printed p. 1),
while the experiment reports a `1361(2) keV H2+` beam (printed p. 3).  If
that molecular-beam energy is divided between the two protons, it corresponds
to about `680.5 keV` per proton and `Ecm ≈ 624 keV`; the figure label gives
`Ecm ≈ 619 keV`.  This internal source convention/rounding mismatch is why
the evidence below identifies the point only as the broad `~0.62 MeV cm`
resonance and retains the source's 675-keV-lab label.  It must not be read as
a 1361-keV single-proton beam.

## Channel and angular-momentum evidence

| parent state and source | channel | source-backed allowed primary `l` | secondary `l'` | branch evidence actually reported |
|---|---|---:|---:|---|
| 16.11 MeV `2+ , T=1`; Laursen 2016, printed pp. 2–3 | `α0`: `8Be(gs,0+)` | `l=2` | `l'=0` | Laursen corrected `8Be(gs)` branch: `5.4(1.1)%` from multiplicity-2 and `5.1(5)%` from multiplicity-3 (printed p. 10, Sec. 5.2). Their values exclude the reported `8Be(gs)` ghost contribution (about 20% in their Ref. 33). |
| same | `α1`: `8Be(2+)` | `l=0,2,4` by spin/parity | `l'=2` | No independent α1 percentage in Laursen. Their non-ground-state sample is called α1 only as a convenient label until the sequential interpretation is established (printed p. 7, Sec. 4.2). The measured shape supports primary `d` wave (`l=2`) dominance; allowed `s,d,g` remain the formal set. |
| 16.11 MeV `2+ , T=1`; Munch & Fynbo 2018 | same two channels | agrees with the above | agrees with the above | At the resonance, they report `σ(p,α0)=2.03(14) mb`; triple-coincidence `σ(p,α1)=38(5) mb` using a Laursen-model acceptance; recommended `σ(p,α1)=39(3) mb`, total `σ(p,α)=41(3) mb`, and reported ratio `σ(p,α1)/σ(p,α0)=19.9(14)` (printed pp. 5–7, Eqs. 12–16 and 24–25). These are channel cross sections at the low resonance, not an off-resonance branching curve. |
| broad 16.62/16.57 MeV `2− , T=1`; Kuhlwein 2021 | `α0`: `8Be(gs,0+)` | **forbidden** for an isolated `2−` parent | n/a | Fig. 1 explicitly says the `8Be` ground-state contribution is forbidden by angular momentum and parity (printed p. 1). No α0/α1 branch ratio is reported. |
| same | `α1`: `8Be(2+)` | `l=1,3` | `l'=2` | Kuhlwein's exclusive Dalitz analysis fits a coherent `l=1 + l=3` mixture. Their Table 2 gives `85(1)%/15(1)%` for the constrained phase-zero fit and `76(5)%/25(5)%` with free relative phase `δ=67(1)%·2π` (printed p. 5). These percentages are fitted partial-wave weights, with statistical errors only; they are **not** α0/α1 branching fractions. |

The formal `l` sets follow directly from the source angular-momentum and
parity assignments: α and `8Be(gs)` have `0+`; `8Be(2+)` has `2+`; and the
`2+ → α+α` secondary decay requires `l'=2`.  For a `2+` parent, the gs
channel therefore has `l=2`, while coupling `8Be(2+)` to total `J=2` permits
even `l=0,2,4`.  For a `2−` parent, the gs channel would require the same
`J` coupling but has the wrong parity, whereas the `2+` intermediate permits
odd `l=1,3`.  Kuhlwein's Fig. 1 supplies the direct experimental-paper
statement of the gs prohibition; Laursen's Sec. 2 supplies the low-state
allowed sets verbatim.

## What is measured versus what remains unknown

### Low 16.11-MeV resonance (`Ecm ≈ 148–149 keV`)

* Laursen's complete-kinematics run used 167–170 keV lab protons, detector
  energy resolution 40 keV FWHM, and 10-keV single-α spectrum bins.  The
  `8Be(gs)` sample is identified event by event from an αα relative-energy
  peak at 92 keV; non-gs events are only provisionally labelled α1.  The
  α1 multiplicity-3 efficiency is only 0.2–0.8% and is obtained from their
  selected sequential model (printed pp. 5–7 and p. 10).
* Laursen's robust direct branch result is therefore the small gs/α0 branch,
  `5.1(5)%` (multiplicity-3), with the stated ghost omission.  Treating
  `100% − 5.1%` as an independently measured α1 branch would overstate what
  the paper establishes; it is only a conditional complement after adopting
  the two-channel bookkeeping.
* Munch & Fynbo's all-three-particle measurement removes the old single-
  detector normalization ambiguity.  It provides the strongest low-resonance
  channel numbers above, but explicitly says the α1 acceptance and resulting
  `38(5) mb` cross section are model dependent.  Their `39(3) mb` value is a
  recommendation combining their result with an earlier measurement, not a
  new energy-dependent α1 data series.

### Broad 2− resonance (`Ecm ≈ 0.62 MeV`, source label `Ep ≈ 675 keV lab`)

* Kuhlwein collected about `6×10^5` exclusive 3α events, required each α to
  exceed 250 keV, removed the `8Be(gs)` locus, and fitted the remaining
  Dalitz distribution.  The fitted coherent `l=1/3` mixture is evidence about
  the primary orbital waves in the α1 channel; its 76/25 or 85/15 fit weights
  cannot be substituted for an α1 decay branch.
* Kuhlwein says roughly 50% of the 2− state's decay goes to 3α and the rest
  returns to `p+11B` (printed p. 1).  This is a total decay-mode statement,
  not a measured `(p,α1)/(p,α0)` cross-section ratio; the isolated 2−
  selection rule makes the 3α branch α1 in this sequential assignment.
* No primary source in the supplied set reports a broad-resonance α0/α1
  branch fraction or angle-integrated α1 cross section.  The source reports
  no machine-readable energy-binned α1 yield from which one could construct
  such a fraction.

## Cross-section coverage and the off-resonance boundary

Munch et al. 2020 measured **α0 only**, angle integrated, for `Ep_lab =
0.5–3.5 MeV`; their Table 1 includes `σ(p,α0)=4.3(5) mb` at 683 keV lab
(printed p. 4).  The paper shows a sharp α0 peak and a broad α1/secondary
distribution in the CM energy spectrum (Fig. 3, printed p. 2), but did not
extract an α1 cross section.  It states that the α0 data alone cannot decide
whether the below-2-MeV strength is subthreshold `2+/3−` ghosts or a broad
`1−` contribution plus a broad `0+` tail, and calls for α1 and Dalitz data
(printed p. 6).  Consequently, the 683-keV-lab α0 number is not an α0 branch
of the isolated broad `2−` state; α0 is forbidden for that state.

Sikora & Weller 2016 is an evaluation, not a primary branch measurement.  It
uses a sequential-spectrum ansatz with primary `l=3` at `Ep_lab=0.675 MeV`,
divides total measured α yield by three, and treats α0 separately only above
2 MeV (printed pp. 1–3).  It is useful as a model/evaluation citation but
cannot close the missing broad-resonance α1 branch or supply off-resonance
branch fractions.

## Primary source index

1. K. L. Laursen et al., “Complete kinematical study of the 3α breakup of
   the 16.11 MeV state in 12C,” *Eur. Phys. J. A* **52**, 271 (2016).
   [arXiv:1604.01244](https://arxiv.org/abs/1604.01244),
   [doi:10.1140/epja/i2016-16271-2](https://doi.org/10.1140/epja/i2016-16271-2).
   Relevant: printed pp. 2–3 (allowed `l,l'` and energy convention), p. 6
   Table 2 (beam energies), pp. 7–8 (channel identification and spectra),
   p. 10 Sec. 5.2 (branching), p. 11 conclusion.
2. M. Munch and H. O. U. Fynbo, “The partial widths of the 16.1 MeV 2+
   resonance in 12C,” *Eur. Phys. J. A* **54**, 138 (2018).
   [arXiv:1805.10924](https://arxiv.org/abs/1805.10924),
   [doi:10.1140/epja/i2018-12577-3](https://doi.org/10.1140/epja/i2018-12577-3).
   Relevant: p. 1 Table 1; pp. 4–6 Eqs. (6), (10)–(16); p. 7 conclusion.
3. M. Kuhlwein et al., “Exclusive decay study of the 16.62 MeV (2−, T=1)
   resonance in 12C,” *Phys. Lett. B* **825**, 136857 (2022).
   [arXiv:2109.07886](https://arxiv.org/abs/2109.07886),
   [doi:10.1016/j.physletb.2021.136857](https://doi.org/10.1016/j.physletb.2021.136857).
   Relevant: p. 1 Fig. 1 and 2− selection rule; p. 3 experiment/data cut;
   p. 5 Table 2 and conclusion.
4. M. Munch et al., “Resolving the 11B(p, α0) Cross Section Discrepancies
   between 0.5 and 3.5 MeV,” *Eur. Phys. J. A* **56**, 14 (2020).
   [arXiv:1908.04064](https://arxiv.org/abs/1908.04064),
   [doi:10.1140/epja/s10050-019-00016-8](https://doi.org/10.1140/epja/s10050-019-00016-8).
   Relevant: pp. 1–2 (state assignments and α0/α1 spectral labels), p. 4
   Table 1, p. 6 conclusion/interpretation limit.
5. M. H. Sikora and H. R. Weller, “A New Evaluation of the 11B(p,α)αα
   Reaction Rates,” *J. Fusion Energy* **35**, 538–543 (2016),
   [doi:10.1007/s10894-016-0069-y](https://doi.org/10.1007/s10894-016-0069-y).
   Evaluation/model only; not used as a direct branch-fraction source.

### Evidence boundary for production use

The closed primary-source facts are: low resonance parent `2+`, low α0 gs
branch about 5.1% with a documented ghost omission, low α1/α0 channel cross
sections only with model-dependent α1 acceptance, broad parent `2−`, broad
α0 forbidden, and broad α1 primary waves `l=1,3` with a coherent mixture
favoured by exclusive data.  A channel-resolved branch fraction versus
incident energy outside the low-resonance measurement and the broad-point
partial-wave fit is absent from these primary sources.

## Follow-up source beyond the bounded source set above

Root located Taskaev etal2024, DOI10.1016/j.nimb.2024.165490, with an
institution-hosted full text at
https://bnct.inp.nsk.su/publics/2024/NIMB_2024_3a.pdf . The completed
[TASKAEV_2024_CHANNEL_EVIDENCE](TASKAEV_2024_CHANNEL_EVIDENCE.md) audit and
132-row transcription cover both channels. These are target-depth averaged
proton energies; angle integration uses a two-angle angular ansatz. They
extend channel-resolved coverage but are not unfolded pointwise or
resonance-separated partial-wave data. The absence statements above apply
to the older enumerated source set, not all publications.
