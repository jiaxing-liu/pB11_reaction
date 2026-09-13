# Target-burn validation commands

The target-burn C++ and Fortran bindings were configured from
`/home/cloud/research/pB11_reaction` with `BUILD_TESTING=ON`.

```text
cmake -S /home/cloud/research/pB11_reaction -B /tmp/m3-target-burn-gnu -DPB11_BUILD_FORTRAN=ON -DBUILD_TESTING=ON
cmake --build /tmp/m3-target-burn-gnu -j2
ctest --test-dir /tmp/m3-target-burn-gnu --output-on-failure
```

GNU CTest: 25/25 passed; output is in `m3-target-burn-gnu-ctest.log`.

For Intel ifx r4 and r8, the same configure/build/CTest sequence was run
after `source /opt/intel/oneapi/setvars.sh`.  The build directories were
`/tmp/m3-ifx-r4` and `/tmp/m3-ifx-r8`; both suites passed 25/25.  Their CTest
outputs are in `m3-target-burn-ifx-r4-ctest.log` and
`m3-target-burn-ifx-r8-ctest.log`.

For the no-Fortran configuration:

```text
cmake -S /home/cloud/research/pB11_reaction -B /tmp/m3-no-fortran -DPB11_BUILD_FORTRAN=OFF -DBUILD_TESTING=ON
cmake --build /tmp/m3-no-fortran -j2
ctest --test-dir /tmp/m3-no-fortran --output-on-failure
```

The no-Fortran suite passed 16/16; output is in
`m3-target-burn-no-fortran-ctest.log`.

The installed C11 smoke caller was compiled against
`/tmp/m3-target-burn-installed/include` and linked with its installed
`lib/libpb11.a`.  The source and exact commands/output are in
`m3-target-burn-c11.c` and `m3-target-burn-c11.log`.

Inspected Intel cache settings: Release, ifx 2026.1, Fortran flags
`-check all -traceback` (default REAL) and `-r8 -check all -traceback` (r8).
For a fresh Intel build, pass these explicitly with
`-DCMAKE_Fortran_COMPILER=/opt/intel/oneapi/compiler/2026.1/bin/ifx` and
`-DCMAKE_BUILD_TYPE=Release`. GNU uses empty extra Fortran flags/build type;
no-Fortran uses Release. The companion BALDUR Intel-r8 executable rebuilt
and linked this library; its log is a link check, not source activation.
