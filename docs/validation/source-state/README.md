# Reproduce the source-state validation

Build the library with CMake and run CTest. `PB11_BUILD_FORTRAN=ON` enables
the binding; explicit ISO kinds also work with changed default kinds. Use a
separate build directory for each Fortran compiler/default-kind combination.
The versioned logs include the first failed helper test, the environmental
LeakSanitizer failure, and the successful final matrices.

Run the continuous-source study against the current static library:

```sh
g++ -O2 -std=c++17 -I../../../include study_source_state.cpp /path/to/build/libpb11.a -o /tmp/study_source_state
/tmp/study_source_state > m3-source-state-study.csv 2> m3-source-state-study.log
MPLCONFIGDIR=/tmp/pb-mpl python3 plot_source_state.py
```

The plotting script needs numpy/matplotlib. It writes PNG/PDF next to the
CSV. Copies are in the paired BALDUR documentation. This is prescribed-bath
continuous **external injection**, not a pB spectrum or self-consistent host
calculation. The final rebuilt study reproduces the versioned CSV exactly.

For installed consumers, install the library and matching Fortran module,
then configure `consumer/` with `-DPB11_PREFIX=/path/to/prefix` and the same
Fortran compiler. Intel validation also supplied `-r8 -check all -traceback`.
Both executables must return status zero. They do not rely on source-tree
include directories or module files.
