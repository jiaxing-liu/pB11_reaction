# Reproduce the published-spectrum shape check

See ../../PUBLISHED_SPECTRUM_CHECK.md for the source PDF identity, model
assumptions, distances and limitations. No original PDF or detector-correction
model is distributed here. The extraction constants apply to that exact PDF.

From repository root, with a built library and plotting Python environment:

```
pdftocairo -svg -f 4 -l 4 /path/to/kuhlwein-2021-preprint.pdf /tmp/kuhlwein-p4.svg
python3 docs/validation/published-spectrum/extract_kuhlwein_spectra.py \
  /tmp/kuhlwein-p4.svg /tmp/kuhlwein-fig3-curves.csv
c++ -std=c++17 -O2 -Iinclude docs/validation/published-spectrum/study_published_spectrum.cpp \
  build/libpb11.a -o /tmp/study_published_spectrum
/tmp/study_published_spectrum 128 > /tmp/published-spectrum-n128.csv
/tmp/study_published_spectrum 256 > /tmp/published-spectrum-n256.csv
/tmp/study_published_spectrum 512 > /tmp/published-spectrum-n512.csv
/tmp/study_published_spectrum 512 0.8191625758768345 > /tmp/published-spectrum-raw-coefficient.csv
MPLCONFIGDIR=/tmp/pb-mpl python docs/validation/published-spectrum/compare_published_spectrum.py \
  /tmp/published-spectrum-n512.csv /tmp/published-spectrum-final \
  /tmp/kuhlwein-fig3-curves.csv /tmp/published-spectrum-raw-coefficient.csv
```

`study_published_spectrum.cpp` outputs per-event grid populations; stderr
records spill and normalization. Display comparison normalizes the retained
plot interval only and never changes these original populations. The plot
is a bounded reference comparison, not an EXL result or a new empirical fit.
An initial transcription of the raw-coefficient hypothesis used0.819157327;
that exploratory result was superseded by the formula-derived0.8191625758768345
reported here. No production code or tolerance was changed in this study.
