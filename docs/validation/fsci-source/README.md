# FSCI validation reproduction

See ../../FSCI_SOURCE.md for equations, domains and physical limitations.
Run the CMake test suite with PB11_BUILD_FORTRAN=ON under GNU and Intel;
Intel default and -r8 builds use -check all -traceback. The C++-only build
uses PB11_BUILD_FORTRAN=OFF. The three *-pass.log files and nof-test-final
log are final results. Earlier *-test-final.log files preserve the outdated
513-node rejection assertion failure; the final rejection boundary is1025.

Regenerate tables using tools/generate_coulomb_tables.py with mpmath1.3.0:

```sh
python3 tools/generate_coulomb_tables.py --jobs 5 --common-radius-fm 16 --table-namespace fusion_fsci_table_data --cache-dir /tmp/pb11-fsci-r16-cache --output-inc /tmp/fsci-regenerated.inc --output-json /tmp/fsci-regenerated.json
```

The successful generation-spawn log and standard-CLI cached regeneration
log are retained; the first generation log preserves the forkserver socket
failure that motivated explicit spawn. The regenerated include matches the
production table byte for byte. Direct-reference Python/CSV use50 digits.

Compile study_fsci_spectrum.cpp against the library for the source scans;
its command-line arguments and output columns are defined in that source.
The n128 through n1024 CSVs preserve numerical refinement. The cut10 files
use10keV instead of1keV cutoff. The mix and fixed-raw variants distinguish
normalization conventions; see fsci-mixture-convention.json.
analyze_fsci_spectrum.py records paths and computes the reported metric;
change its local input/output paths when reproducing elsewhere.

consumer/ contains the installed-target C11 and Intel-r8 consumer sources;
configure with PB11_PREFIX pointing to the installation. Build and run both
consumer and fortran_consumer. Final log exits were checked individually.
