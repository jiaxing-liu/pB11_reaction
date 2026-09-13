# Taskaev et al. (2024) p+11B channel-data audit

Source: S. Taskaev et al., “Measurement of the `11B(p,α0)8Be` and the
`11B(p,α1)8Be*` reactions cross-sections at the proton energies up to
2.2 MeV,” *Nuclear Instruments and Methods in Physics Research B* **555**
(2024) 165490.  The institution-hosted [full PDF](https://bnct.inp.nsk.su/publics/2024/NIMB_2024_3a.pdf)
and [DOI record](https://doi.org/10.1016/j.nimb.2024.165490) are the primary
source.  The transcription of Tables 3–8 is versioned as
[taskaev-2024-cross-sections.csv](validation/birth-channel/taskaev-2024-cross-sections.csv).
The downloaded PDF and text extraction are temporary evidence caches; the
institutional link above is the reproducible source.

This is a channel-cross-section measurement over a broad energy scan.  It
does not assign the scan to individual `12C` resonances or report new
parent-state `Jπ` or primary `l` values.  It therefore extends the
channel-resolved data coverage, but it does not replace the resonance-specific
branch and partial-wave evidence in Laursen/Kuhlwein.

## Energy and reaction convention

The accelerator proton beam is stated to cover 0.1–2.2 MeV, in 100-keV
steps.  The paper's Tables 2–8 use `E` as the **average kinetic energy of a
proton interacting with boron nuclei** and `ΔE` as the standard deviation of
that energy (Table 2 caption, printed p. 7).  Thus the table energies are
proton lab-frame energies at the reaction depth, averaged over the target
energy loss; they are not the accelerator setting at the target surface.  The
authors state that the proton loses about 28 keV at 2.2 MeV and 150 keV at
0.15 MeV in the boron/nitrogen layer, and explicitly warn that Table 2 gives
the average energy in the layer rather than the incident surface energy
(printed p. 9).  The measured setting differs from the set value by at most
2 keV and the beam stability is 1–2 keV, but this does not remove the target
energy spread.

For comparison with a stationary-target entrance convention, one may use
`Ecm ≈ (11/12) Ep_lab` for a proton on `11B`; no such conversion has been
applied to the tables here because the published `E` is the in-target
average.  In particular, the low `E=75 ± 75 keV` row is a broad target-energy
average, not a resolved 148-keV-cm resonance point.

## What the paper calls α0 and α1

The reaction equations and energy peaks are given on printed p. 7:

* `α0` means the primary α in `11B(p,α0)8Be(gs)`.  It appears as a sharp
  peak near 5.5 MeV; the two secondary α particles from `8Be(gs)` cluster
  near 1.4 MeV.
* `α1` means the primary α in `11B(p,α1)8Be*(2+)`.  It appears as a broad
  peak near 3.7 MeV; its two secondary α particles cluster near 2.4 MeV.
* The paper writes the sequential energy releases as
  `11B+p → α0 + 8Be + 8.587 MeV`, `8Be → α01 + α02 + 0.094 MeV`,
  `11B+p → α1 + 8Be* + 5.557 MeV`, and
  `8Be* → α11 + α12 + 3.124 MeV`.  It also quotes 8.681 MeV for the
  total `11B(p,3α)` release.

The α1 table is a reconstructed **total α-particle yield for the α1
channel**, converted to a reaction cross section with `N=3` in the yield
formula.  Up to and including 1.4 MeV proton energy, they count the α
particles above the line-B energy (about 2.4 MeV), multiply by `3/2` because
two of the three α particles are above that line, and use that to estimate
the total α1-channel yield.  Above 1.4 MeV, backscattered protons enter this
region, so they count from the α1 maximum to the α0 peak and multiply by
`2.95 ± 0.05`, calibrated from spectra at energies ≤1.4 MeV (printed p. 8).
This is a channel-yield reconstruction, not an event-by-event tag of every
primary α1.

The paper states that the observed spectrum is predominantly sequential and
that direct decay, if present, is much smaller.  It does not provide a direct
decay cross-section measurement.

## Extracted integrated channel cross sections

Tables 5 and 8 give angle-integrated channel cross sections in **mb**.  `E`
and `ΔE` retain the in-target proton convention above.  The `Δσ` columns are
the paper's reported statistical variance/uncertainty values; the paper also
gives larger combined accuracy budgets below.

| E (keV, in-target mean) | ΔE (keV) | α0 σ (mb) | α0 Δσ (mb) | α1 σ (mb) | α1 Δσ (mb) |
|---:|---:|---:|---:|---:|---:|
| 75 | 75 | 0.017 | 0.008 | 0.56 | 0.23 |
| 134 | 66 | 0.27 | 0.05 | 6.4 | 0.8 |
| 247 | 53 | 0.64 | 0.10 | 40 | 5 |
| 355 | 45 | 1.40 | 0.24 | 148 | 19 |
| 461 | 39 | 2.02 | 0.34 | 357 | 42 |
| 565 | 35 | 2.53 | 0.42 | 598 | 88 |
| 668 | 32 | 3.31 | 0.55 | 668 | 89 |
| 771 | 29 | 3.45 | 0.57 | 386 | 58 |
| 873 | 27 | 3.04 | 0.50 | 234 | 34 |
| 975 | 25 | 2.90 | 0.49 | 171 | 22 |
| 1077 | 24 | 2.61 | 0.44 | 152 | 20 |
| 1178 | 22 | 2.27 | 0.39 | 156 | 19 |
| 1279 | 21 | 2.07 | 0.37 | 161 | 18 |
| 1380 | 20 | 1.66 | 0.31 | 160 | 17 |
| 1481 | 19 | 1.27 | 0.24 | 146 | 22 |
| 1582 | 18 | 0.97 | 0.18 | 129 | 21 |
| 1683 | 17 | 1.35 | 0.24 | 117 | 22 |
| 1783 | 17 | 1.96 | 0.34 | 110 | 18 |
| 1884 | 16 | 2.43 | 0.40 | 102 | 17 |
| 1985 | 16 | 8.08 | 1.30 | 102 | 20 |
| 2085 | 15 | 9.00 | 1.52 | 84 | 20 |
| 2186 | 14 | 4.09 | 0.68 | 84 | 19 |

Examples directly relevant to the broad-resonance energy region are the
`E=565 ± 35 keV` row (`α0=2.53 ± 0.42 mb`, `α1=598 ± 88 mb`) and the
`E=668 ± 32 keV` row (`α0=3.31 ± 0.55 mb`, `α1=668 ± 89 mb`).  These are
measured channel cross sections averaged over the target energy distribution;
they are not a decomposition into the isolated 16.62-MeV `2−` resonance and
other amplitudes.

The CSV also contains the differential tables, with exact units and angles:

* Table 3: α0, CM differential cross section at 135°, `mb/sr` (printed p. 8).
* Table 4: α0, CM differential cross section at 168°, `mb/sr` (printed p. 8).
* Table 5: α0, angle-integrated cross section, `mb` (printed p. 10).
* Table 6: α1, CM differential cross section at 135°, `mb/sr` (printed p. 10).
* Table 7: α1, CM differential cross section at 168°, `mb/sr` (printed p. 11).
* Table 8: α1, angle-integrated cross section, `mb` (printed p. 12).

## Angular-distribution assumptions and errors

The spectrometers were at 135° and 168° in the laboratory.  The authors
transform the differential yield to the beam-target CM system.  If emission
is isotropic in CM, they use `σ = 4π G dσ/dΩ`.  They instead find visibly
anisotropic radiation in both channels and parameterize it with the even
form

```
dσ/dΩ(θ) = [dσ/dΩ(90°)] [1 + A(E) cos² θ],
σ = 4π [1 + A(E)/3] dσ/dΩ(90°).
```

`A(E)` is determined from only the two measured angles (α0: printed pp. 9–10,
Figs. 8–10; α1: printed pp. 10–12, Figs. 11–13).  Thus the angle-integrated
values depend on the stated two-angle angular ansatz; they are not an
isotropy-only extrapolation.

The paper's quoted accuracy budgets are:

* α0: boron thickness 10%, proton fluence 1%, solid angle 5%, statistical
  uncertainty 2–30%, combined 12–33% (printed p. 9).
* α1: boron thickness 10%, proton fluence 1%, solid angle 5%, counting the
  high-energy α1 events 3%, statistical uncertainty 0.1–2%, combined 12%
  (printed p. 10).

They assume α detection efficiency `k=1` (100%) after calibration with
standard α sources (printed p. 2).  The reported table uncertainties should
therefore be read together with the channel-specific systematic budgets,
target energy spread, and α1 yield-reconstruction rule.

## Direct-decay statement and its limits

Taskaev et al. do not report a direct-decay cross section.  Their estimate is
given on printed pp. 11–12:

* At `E≈0.6 MeV`, detection angle 168°, they count 117 α particles in the
  6–9 MeV region and assume all are direct-decay events.
* Depending on assumed low-energy continuation, they estimate 351 or 468
  direct-decay α particles, compared with 1,049,811 sequential-decay events
  (`3,936` α0 and `1,045,875` α1).
* They therefore quote a direct probability 2,000–3,000 times lower than the
  sequential probability.

This is explicitly an estimate based on assigning the 6–9 MeV tail to direct
decay and assuming shapes below 6 MeV.  It is not a measured direct-channel
branch fraction and should not be used as one.

## Evidence boundary for resonance use

The paper supplies a useful tabulated α0/α1 channel data set from mean
in-target proton energies 75–2186 keV, with differential data at two angles
and angle-integrated values in mb.  It establishes the authors' sequential
spectral labels and reports a small, assumption-dependent direct-decay
estimate.  It does **not** provide:

* a surface incident-energy-to-`Ecm` table for each row;
* a resonance-separated α0/α1 branch fraction;
* a new assignment of the scan to the low `2+` or broad `2−` `12C` parent;
* primary α orbital-angular-momentum (`l`) amplitudes; or
* a machine-readable direct-decay cross section.

Therefore these tables can document measured channel cross sections across
the scan, including the broad-resonance energy neighborhood, while branch
fractions and `l` choices remain source-limited to the resonance-specific
evidence already audited.
