# Published three-alpha spectrum comparison

This check changes the next physical-source decision. It does not certify a
complete source, refit experimental data, or activate a host model.

## Reference and extraction

Kuhlwein et al., arXiv:2109.07886v1, printed p4 Fig3, contains three plotted
model curves and the same measured single-alpha spectrum. Obtain the source
from https://arxiv.org/abs/2109.07886. PDF SHA256:
`1f521b3332780aab60647addebd19e88b31f1bc3ae00936db2d006ae5c9d3b22`.

Root inspected the rendered page and extracted the original vector line
vertices via pdftocairo SVG, not by pixel tracing. There are100 points per
curve, at approximately30..5970keV. The script's axes and transform constants
are specific to that PDF. Black is selected experimental data; orange is the
published model. The extracted repeated data means agree within0.012keV,
providing an independent check of panel-axis calibration. These curves do not
supply binwise uncertainties or original exclusive events.

The experiment selects three detected alphas above250keV and removes ground-
state Be8 events. Thus the black curves cannot be used as uncorrected source
probabilities or as an exclusive-event likelihood. This check primarily
compares the library to the paper's orange model curves. The paper does not
fully specify the numeric radial parameters or separate basis normalization
needed to claim exact reproduction of its fitted k/phase convention.

## Declared library calculation

Use the existing Laursen-based radial parameters in ALPHA_AMPLITUDES and
ALPHA_SPECTRUM: A=9.3MeV, cutoff1keV,100 bins of60keV over0..6MeV, modes1,
3 and13; mixture k=.76, phase=.67*2pi in the library's UNIT-BASIS convention.
A=9.3MeV is a stated comparison setting near the experimental incident energy,
not a replacement for the new mass-derived Q or an assertion that the nominal
16.62MeV level fixes every event's energy. The observed three-particle mean is
about3.103MeV; the library's full event mean is exactly A/3 before display
truncation. Each original library output retains all grid spill ledgers.

For shape comparison ONLY, normalize each displayed0..6MeV curve to unit area.
No density or reaction-yield normalization is changed in the library. The
shape distance is half the integrated absolute difference of these normalized
curves; it is not chi-square, parameter uncertainty, or a physical error bound.

| Mode | n256 to n512 shape change | n512 vs published model | n512 vs selected data |
|---|---:|---:|---:|
|l1|0.00390|0.00806|0.13166|
|l3|0.00468|0.04606|0.09016|
|mixture unit k=.76|0.00421|0.06564|0.05473|

The refinement column uses discrete cell populations normalized by their sum;
reference comparisons use trapezoidal areas over the plotted sample centers.
At n128, grid-projection oscillations give misleading shape discrepancies of
6–8% even though particle/energy totals are accurate. Refining to256/512 is
necessary for this60keV binning. Moment closure alone is not shape convergence.
At512, mapped particles plus recorded spill equal3; spill is about1.5e-4 to
5.0e-4 particles per reaction across the three modes and is never dropped by
the API. The n512 l3/reference discrepancy remains much larger than the last
quadrature refinement change.

## Test of a normalization hypothesis

If the paper's k were a coefficient of the library's RAW amplitudes with equal
fixed reduced-width scale, then the corresponding unit-basis coefficient is

```
k_unit = k_raw*N1 / (k_raw*N1 + (1-k_raw)*N3).
```

This follows by rewriting the amplitude before interference, not by assigning
post-interference branching fractions. For the current radial model,
N1~4.35297e-26J^2, N3~3.04303e-26J^2; k_raw=.76 gives
k_unit~.8191625759. A direct512 calculation at that value reduces the distance
to the published mixture model from0.06564 to0.05003; it does not eliminate it.
The selected-data distance becomes0.04434. These are hypothesis checks, not
permission to relabel the library parameter as the published fit coefficient.

## Decision for subsequent source composition

Do not freeze .76 as a measured unit-basis mixing fraction. The numerical
amplitudes remain independently tested, but exact experimental-spectrum
acceptance is unproven. Pure l1's close model-curve agreement does not certify
pure l3, the mixed model, or a different parent resonance. Before promoting
a physical source, account explicitly for the radial-model choice, normalized
mixture convention, incident-energy dependence, allowed ground-state branch,
and any final-state interaction approximation. I010 now includes a measured
shape discrepancy, not merely an absent textual normalization statement.

The source-code study, extracted curves, raw128/256/512 spectra, comparisons
and plotting scripts are in docs/validation/published-spectrum/. Full M3/M4
composition and EXL case acceptance remain required.
