# Thermal source increments: compiler and accounting evidence

The complete GNU suite passes59/59 (66.88s). Intel-r8 focused coupled, table,
handoff and Fortran binding tests pass4/4. The no-Fortran increment test passes1/1.
A strict C11 consumer compiles and calls the exported C function successfully;
its exact source is included. This consumer links the built static library,
not an installed package. Prior full matrix coverage is preserved in the birth-
table milestone; this increment does not claim a new full Intel/no-Fortran matrix.

Independent hand-calculated signed ledgers test particle removal/handoff, all
heat baths, explicit fast external/escape exclusion, weak1e-30 signals against
background1e20, nonfinite/negative fields, null inputs and overflow rejection.
Actual DT/carbon direct/table trajectories reconstruct thermal unknowns using
the returned increments within stored-inventory rounding. An actual mixed-pool
alpha handoff checks number transfer and signed heat correction counted once.
The Fortran wrapper checks64-byte layout, signed mapping and incorrect extent
before C_LOC. Numerical sources were frozen before these completed suites.

Reproduce with the existing configured builds:

```sh
cmake --build /tmp/m3-target-burn-gnu -j2
ctest --test-dir /tmp/m3-target-burn-gnu --output-on-failure
cmake --build /tmp/m3-ifx-r8 --target test_fusion_thermal_increment test_fusion_coupled_table test_fusion_coupled_thermal test_fusion_coupled_thermal_fortran -j2
ctest --test-dir /tmp/m3-ifx-r8 -R 'fusion_(thermal_increment_tests|coupled_table_tests|coupled_thermal_tests|coupled_thermal_fortran_binding_tests)' --output-on-failure
```

Intel requires sourcing oneAPI and its existing bounds-check/r8 configuration.
No new physical approximation or host state advancement is introduced here.
