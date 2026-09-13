# Coupled-thermal build checks

Captured 2026-09-14 (CST) from the four requested build trees. The complete
CTest suites passed: GNU 53/53, Intel ifx-r4 53/53, Intel ifx-r8 53/53, and
no-Fortran 32/32. Each tree's `Testing/Temporary/LastTest.log` is retained
under its matching subdirectory.

The installed consumer records an actual C11 DT birth and thermal-feedback
check, plus the Intel Fortran consumer check. `installed-consumer/` contains
only its `CMakeLists.txt`, `main.c`, and `main.f90`; no consumer build products
are archived. The host record is a link-only build of the pB11 library into
the host executable; host state activation or runtime execution was not
performed.

Compiler diagnostic scan:

- `ifx-r4-build.log`: 23 instances of `ifx: warning #10182` (debug checks
  disable optimization).
- `ifx-r8-build.log`: 23 instances of the same expected debug-build warning.
- `consumer-build.log`: 1 instance of the same warning.
- `host-build.log`: 1 `ifx: remark #10440` about debug options without an
  optimization level.
- No other warning or error diagnostics were found in the archived configure,
  GNU, no-Fortran, install, or consumer-test logs.

The archived `representative-*.log` files are the final passing representative
test. No separate pre-fix coupled-thermal failure log was present in `/tmp`
at capture time. The known pre-fix status-4 report came from the particle
ledger roundoff check omitting old/new seeded fast inventories from its ULP
budget; it was not a physical negative thermal-energy result.
