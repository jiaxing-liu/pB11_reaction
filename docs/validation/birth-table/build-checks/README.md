# Build, binding and installed-consumer evidence

Complete suites pass: GNU58/58, Intel default58/58, Intel-r858/58, and
no-Fortran36/36. Each subdirectory retains configure/build/CTest and full
LastTest logs. Intel configurations use `-check all -traceback`; the r8 build
also uses `-r8`. Separate TMPDIRs avoid compiler scratch-file collisions.
The first focused three-test run passed before expanding to these suites.
Only expected ifx #10182 debug/optimization warnings were found in these
build logs; the host build also retains its usual debug-option remark.

The C11 installed consumer creates a real table and invokes a DT coupled trial
with carbon. The Intel-r8 installed Fortran consumer independently creates a
table through the Fortran wrapper, passes a five-handle array to the coupled
wrapper and checks actual fuel depletion, events and carbon heat. Both pass.
Their source and CMakeLists are retained. They include exported installed
targets under `${PB11_INSTALL}/lib/cmake/pb11_reaction/pb11Targets.cmake`,
with the Fortran linker selected for the Fortran caller.

Host `bash tools/build_intel.sh` passes with the library changes. This proves
build/link compatibility; the new table/kinetic source owner is NOT activated
inside BALDUR by this milestone. Actual local table/direct evolving runs are
retained in the sibling `studies` directory, with their own binary provenance.

The source hashes in manifest.json were captured after these build/test runs,
with no subsequent numerical source edits. They are separate
from earlier microbenchmark binary hashes; the latter predate two added info
diagnostic fields. The full direct API regression passes after the coupled
source-provider refactor. Tests also preserve explicit rejection of a100keV
truncated DT integration window and of real early cooling below a table that
starts exactly at the initial20keV.
