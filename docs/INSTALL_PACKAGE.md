# Installing the independent library

The CMake project installs the C/C++ library, public headers, optional Fortran
modules, exported targets, and a versioned package configuration. A consumer
can use the installed package without knowing the source tree:

```sh
cmake -S . -B build -DPB11_BUILD_FORTRAN=ON
cmake --build build -j1
cmake --install build --prefix "$PWD/prefix"
```

The installed package is discovered with:

```cmake
find_package(pb11_reaction 1.0 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE pb11::pb11)
```

The package target keeps the stable C ABI (`pb11_c.h`) and the public C++ API
available through the installed include directory. Fortran targets are
exported when a Fortran compiler was enabled during the producing build. The
installation smoke test in `tests/test_install_package.py` builds a clean
consumer and executes a finite-rate call; it does not claim BALDUR or native
full-window acceptance.

## Installed Fortran consumer

Use a compatible Fortran compiler for both the library's modules and the
consumer. Fortran `.mod` files are compiler-specific; keep separate build and
installation prefixes for GNU Fortran and Intel ifx.

For an Intel build, activate oneAPI and select the compiler explicitly:

```sh
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-ifx -DPB11_BUILD_FORTRAN=ON \
  -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release
cmake --build build-ifx -j1
cmake --install build-ifx --prefix "$PWD/prefix-ifx"
```

A standalone Fortran project discovers and links the installed target:

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_consumer LANGUAGES Fortran)
find_package(pb11_reaction 1.0 CONFIG REQUIRED)
add_executable(my_consumer main.f90)
target_link_libraries(my_consumer PRIVATE pb11::fortran)
```

Configure it with `-DCMAKE_Fortran_COMPILER=ifx` and
`-DCMAKE_PREFIX_PATH=/absolute/path/to/prefix-ifx`. The exported target carries
the module include directory and static C++ runtime dependency.

The concrete consumer source and CMake file in `docs/validation/` accompany
`ifx-validation-20261018.json`. They were compiled from the installed package
and executed with ifx 2026.1.1, requiring a successful status and a finite,
strictly positive pB rate at 100 keV. The same build passed all 120 configured
C++/Fortran CTest tests. This is independent-library/interface evidence, not
BALDUR full-window or device-physics acceptance.
