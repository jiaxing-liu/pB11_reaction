# Effective pB population weighting reproduction

See ../../PB_POPULATION.md for physical assumptions and diagnostic meanings.
The20-row numeric table is derived from the already audited birth-channel CSV.
The75/134keV target-averaged rows are deliberately excluded, not discarded
from the evidence repository.

Build/run study_pb_population.cpp against include/ and the built libpb11.a;
stdout is pb-population-study.csv and stderr is pb-population-study.log.
The CSV preserves12 physical input pairs x5 outputs.

Run reference_pb_population.py with mpmath1.3.0 in this directory. It reads
the study CSV, independently integrates40 component cases and writes
pb-population-reference.json. The unscaled log preserves the first
reference's tiny-integral relative-accuracy failure; the scaled log is final.
Plot with plot_pb_population.py (numpy/matplotlib), which writes PNG/PDF.

The four *-test.log files are final full-suite passes. consumer/ is an
installed C11/Fortran CMake project: set PB11_PREFIX, ifx and -r8, build and
run both executables. Run logs were checked for terminal success.
