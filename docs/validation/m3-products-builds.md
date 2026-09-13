# Product increment build matrix

GNU15.2 Fortran: 17/17; Intel ifx2026.1 default and -r8: 17/17 each;
C++ without Fortran: 12/12. Intel tests use -check all -traceback.
The product C++ test includes the independent 100-step inventory check.

Configure source with cmake -S . -B BUILD -DCMAKE_BUILD_TYPE=Release,
then cmake --build BUILD -j2 and ctest --test-dir BUILD --output-on-failure.
For Intel, source /opt/intel/oneapi/setvars.sh --force first and set
-DCMAKE_Fortran_COMPILER=ifx plus the flags above. For no-Fortran, use
-DPB11_BUILD_FORTRAN=OFF. All module sources are in the repository.

Fresh installed-library smoke: configure with -DBUILD_TESTING=OFF and
-DCMAKE_INSTALL_PREFIX=PREFIX; build and cmake --install BUILD. Compile the
versioned m3-products-smoke.c with cc -std=c11 -Wall -Wextra -Werror
-I PREFIX/include; link its object with c++ PREFIX/lib/libpb11.a -lm.
Run the resulting executable. The caller checks both public product APIs
using only installed headers/library. No BALDUR dependencies are required.
