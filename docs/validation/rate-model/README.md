# Explicit rate-model numerical validation

`study_rate_model_velocity.cpp` integrates directly over target speed and its
angle to the beam, independently of the production relative-speed/Langevin
kernel. It calls the declared cross-section evaluator (so this is an integral
check, not an independent nuclear-data fit). Cases exercise pB and DHe3 beam
rates, DD below its low fit endpoint, and DT across its upper endpoint.
The artificial10MeV DT case is a numerical test, not an EXL operating point
or a claim of accurate classical kinetics at that energy.

The CSV compares rate, reaction-conditioned target-energy moment and CM-energy
moment. Max relative difference7.25e-12. Target Maxwell speed integration
extends10 thermal speeds; tails at that point are negligible for these models.

Reproduce from repository root after building the library:

```
c++ -std=c++17 -O2 -Iinclude docs/validation/rate-model/study_rate_model_velocity.cpp \
  build/libpb11.a -o /tmp/study_rate_model_velocity
/tmp/study_rate_model_velocity > /tmp/rate-model-velocity.csv
```

The C++ regression test separately compares thermal segments to the independent
energy-space results frozen in ../nuclear-windows/. The existing-tests log
records compatibility immediately after the shared-kernel refactor; final
compiler-matrix logs in this directory supersede it for acceptance.

The expanded new-API common-T grid has65 cases at.01..500keV. Independent
reference rate differences: fit1.22e-14, low1.88e-14, high1.37e-13; conditional
relative-energy means differ by less than1.36e-15. Compile study_model_grid.cpp
as above, then run compare_rate_model_grid.py with its CSV and the versioned
../nuclear-windows/nuclear-window-study-refined.csv as arguments. A1e-300
absolute comparison denominator protects underflowed rate tails; every
reported nonzero worst case is above that floor.

Final matrix: GNU38/38, Intel default38/38, Intel-r8 with runtime checks38/38,
noFortran24/24. Logs are included. Installed C11 and Intel-r8 Fortran consumers
both returned exit0. Reproduce consumer build with CMake argument
-DPB11_PREFIX=<installation prefix>; source oneAPI and select ifx/-r8 when
using Intel-compiled modules. It imports the installed pb11Targets.cmake,
not build-tree include/module paths. BALDUR Intel-r8 link is recorded by the
paired host commit; no host activation or reactor simulation is claimed.
