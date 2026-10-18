Isolated transfer preparation prototype; all writes confined to this directory.
Public scalar validation ordering and output clearing retained unchanged.
Only fixed projectile/bath log normalization is prepared once; energy ratio
retains (energy/temperature)*(bath_mass/mass), long double arithmetic,
left-associated normalization and exp(normalization-y2). No fast math/LTO.
Coupled uses double(centers[j]) and original bath-wise accumulation order.
The internal prepare/evaluate helpers require validated inputs; scalar validates
as before, while coupled has checked edges, logs, species, temperatures and baths.
No public ABI, physical threshold, mesh, Coulomb diffusion, or tail policy change.

reference.sha256 identifies untouched transfer_reference_frozen.cpp.
transfer_reference_link.cpp is identical apart from renamed public symbol.
prepared-transfer.patch contains all proposed production changes (3 files).
manifest.json identifies proposed sources and fixed 509-edge input.
compare_bench.cpp reproduces comparison and six paired transfer-only timings.
summary.json and compare_bench.log record results. Benchmark includes preparation
once per projectile/bath/repetition. 6 projectiles x 7 ion baths, 508 actual centers,
100 repetitions per mode/pair; alternating mode order, one warm-up each.
Baths are representative cold DD/helium/carbon values, not captured runtime baths.
Transfer-path timings are not a full-coupled or host speedup claim.

Independent build: cmake -S source -B build -DCMAKE_BUILD_TYPE=Release
 -DPB11_BUILD_FORTRAN=OFF -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF
Build selected targets with cmake --build build -j2.
Affected-test output is in affected-tests.log.
