# Coherent spectrum numerical-layer validation

GNU/Intel default/Intel-r8:21/21; no-Fortran:14/14. Intel uses
-check all -traceback and -r8 for the r8 matrix. Commands follow the earlier
m3-products-builds.md. Nuclear table generation uses pinned mpmath1.3.0:
python tools/generate_coulomb_tables.py --jobs 5
The checked-in table is produced with CACHE_VERSION3 after moving precision
initialization before all numerical constants. All28 initial tasks completed;
final channel segment counts are34/34/34/19/19. Off-node phase errors are below
9e-8, logP errors below3.4e-15, scaled shift errors below1.7e-14. These are
sampled errors, not rigorous interval bounds. JSON records all reference values.

Reproduce the independent35-point,45-digit direct-derivative fixture with
python tools/generate_nuclear_reference.py. The compiled test checks these
values and27events for cyclic relabeling, secondary exchange and rotational
spin-norm invariance, plus the independently derived unsymmetrized polynomial.
This establishes the stated algebra/numerical implementation; it does not
resolve the source-level alternate angle formula or establish a complete
pB incident-energy-dependent normalized spectrum.

Install using cmake --install BUILD --prefix PREFIX. Compile m3-spectrum-smoke.c
with cc -std=c11 -Wall -Wextra -Werror -I PREFIX/include; link the object using
c++ PREFIX/lib/libpb11.a -lm and run. Root also rebuilt the BALDUR Intel-r8
executable after adding src/*.inc to its external-library prerequisites;
this confirms build integration, not activation of new stateful fusion physics.
