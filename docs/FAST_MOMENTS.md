# Independent fast-ion moments

`fusion_c_fast_moments` and Fortran `fusion_fast_moments` integrate both S and T
cell populations. Results have six number densities, six midpoint kinetic
energy densities, sum(Z n), sum(Z^2 n) and isotropic NR pressure=2/3 sum(U).
The 120-byte C record uses SI throughout. Charges come from canonical nuclear
data; no BALDUR indices, geometry, thermal pool, default electron density or
rest energy enter this function. S/T repartition at fixed sum preserves all
these moments. T remains kinetic until actual handoff.

This is the moment of the existing arithmetic-center finite-volume state,
not a subcell spectrum reconstruction. It deliberately uses the same NR scalar
pressure convention as the energy FP closure and legacy host fast pressure.
It does not claim tensor/relativistic pressure or infer momentum anisotropy.
A caller using a different kinetic closure must supply appropriate moments.
All values are finite/nonnegative; output overflow or nonzero underflow
rejects and clears the complete record. Zero populations return exact zero.

C++ tests hand-count all six charges and nonuniform-grid energies; check S/T
sum invariance, negative populations, nonincreasing edges, overflow, weak
energy underflow and zero state. Fortran tests verify the 120-byte ABI, alpha
charge/pressure and pre-C_LOC extent checking with cleared results.

GNU full matrix passed 62/62. After correcting the newly added Intel test
link settings to POSITION_INDEPENDENT_CODE ON and LINKER_LANGUAGE Fortran,
final GNU and Intel-r8 focused tests each passed 2/2. Installed header/library
strict C11 consumer compiled, linked and ran successfully. Logs and final
source hashes are in validation/fast-moments.

Two initial failures are retained: the first overflow test incorrectly used
a single max-density proton whose returned moments remained representable;
it was changed to five-charge boron to force charge overflow. Intel first
failed at PIE linking, before running the binding; the test target now uses
the same configuration as existing Fortran tests. Neither required changing
the moment formulas. This increment does not implement host epoch selection,
spatial transport, or any new physical simulation.
