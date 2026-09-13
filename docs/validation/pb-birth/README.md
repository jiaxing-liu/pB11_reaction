# Conditional pB CM birth validation

See ../../PB_BIRTH.md for physics, source conventions and outstanding work.

Independent high-precision test (Python with mpmath1.3.0):

```sh
c++ -std=c++17 -O2 -I ../../../include alpha0_reference_driver.cpp /path/to/build/libpb11.a -o alpha0_reference_driver
python3 alpha0_reference.py
```

Run from docs/validation/pb-birth. The script locates
its driver alongside itself and writes alpha0-reference.json. It uses120
fixed-seed cases and80-digit arithmetic with the canonical rounded mass.

Conditional source illustration:

```sh
c++ -std=c++17 -O2 -I ../../../include study_pb_birth.cpp /path/to/build/libpb11.a -o study_pb_birth
./study_pb_birth > pb-birth-study.csv 2> pb-birth-study.log
python3 plot_pb_birth.py
```

Requires numpy/matplotlib for the plot. The alpha1 source is evaluated at1024
nodes per dimension. Plotting a grid projection of a delta peak does not
imply that the physical peak has the displayed width. The5percent branch is
illustrative, not calibrated. No laboratory motion or host feedback occurs.

consumer/ is a separate C11/Fortran CMake project; set PB11_PREFIX to the
installation, CMAKE_Fortran_COMPILER=ifx and CMAKE_Fortran_FLAGS=-r8. Build
and run both consumer and fortran_consumer. Full test matrices and host build
logs are copied here after completion; retained failure logs are identified
separately if any arise.

Final results are *-pass.log and consumer-run.log. Initial *-test.log files
preserve the test-unit rounding failure; consumer-build.log preserves the
missing-installed-header failure. All are explained in ../../PB_BIRTH.md.
