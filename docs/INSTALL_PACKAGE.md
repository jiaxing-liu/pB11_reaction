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
