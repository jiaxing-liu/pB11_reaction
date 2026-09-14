# Inert-heat transactional state validation

The four complete suites pass: GNU55/55, Intel default55/55, Intel-r855/55,
and no-Fortran34/34. Build trees are /tmp/m3-target-burn-gnu,
/tmp/m3-ifx-r4, /tmp/m3-ifx-r8 and /tmp/m3-no-fortran. Configure, incremental
build and full CTest logs are retained for each; the Intel configurations use
`-check all -traceback`, with `-r8` in the latter. Reproduce using CMake and
CTest with the same options; separate TMPDIRs permit parallel compilers.

New C++ tests validate signed inert heat, old-snapshot rejection, promotion
only upon commit, replacing/discarded/failed candidates, byte-level v1/v2
schema and semantic corruption with a recomputed checksum, and continued
restart evolution. The actual coupled DT/carbon test advances8steps/1ms;
at epoch4 it saves a thermal tuple alongside packed kinetic state, restores
both and continues. Final thermal state and serialized kinetic accounts match
the uninterrupted calculation exactly. Carbon heat is about448.984J/m^3;
missing this account makes staging fail. This is not BALDUR restart acceptance.

`source_state_legacy_bytes.cpp` was linked once to the pre-extension static
library from milestone d1a1e6a, and once to the new source-state object ahead
of that archive. Both generated the retained v1 bytes with SHA256
`0ecb57eb1da899b1c16bba4b5e990c09f317d1d253940c061eabe258d7e9498d`.
The fixture has six seeded species and a nonzero signed-bath heat ledger,
and exercises create/begin/stage/commit/pack. Its source and both outputs
are retained so compatibility is not asserted only from a formula for length.

The install prefix is /tmp/m3-inert-state-installed, from the Intel-r8 tree.
The C11 consumer invokes both new APIs and v2 roundtrip; the Intel Fortran
consumer runs the full binding lifecycle suite through installed modules.
Both pass. The consumer source includes the exported `pb11Targets.cmake`
under `${PB11_INSTALL}/lib/cmake/pb11_reaction`, matching current packaging;
there is no `pb11Config.cmake`. Its first incorrect `find_package(pb11)`
configuration and first CXX-linker-selected Fortran link failure are preserved.
Selecting the existing Fortran linker convention fixed the latter. These
were consumer setup errors, not state-API or physics test failures.

The separate manual GNU Fortran default/real8/integer8/combined-kind tests
also passed in the binding agent's work; matrix logs here are the primary
retained build evidence. No new physics approximation is introduced by this
state-storage extension. Host thermal/geometry/transport ownership is still
required before full BALDUR activation.
