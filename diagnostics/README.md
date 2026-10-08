# Optional diagnostic foundations

Enable `PB11_BUILD_DIAGNOSTICS=ON` to build `pb11::diagnostics`. This separate
target currently implements streaming SHA256 and an explicitly owned table
registry with a C interface. It has no reaction, collision or burn formulas.
Its sources remain outside the numerical kernel identity's src/include roots.

The registered unpack hooks import the exact supplied packed bytes using the
original public numerical API, hash the existing loader buffer, bind the new
handle to that identity and publish it only after registration succeeds.
Never use the original destroy functions on these owned handles. Use the
matching diagnostic destroy hook or release the complete context. Evaluation
may borrow a handle only while its context/registration remains live. Context
destruction cannot race any operation or numerical use of its owned handles.

The registry requires an absolute qualified file reference; it records the
supplied bytes' identity rather than claiming to have authenticated that file.
A future replay reader must verify file content, size and kernel identity
before importing it. Registration is not interpolation or physical validation.
Unpack rejection and registration failures have separate diagnostic status
codes; the actual original unpack result is reported separately.

With BUILD_TESTING enabled, SHA tests run without table assets. To register
the real-table lifecycle test, explicitly supply
`PB11_DIAGNOSTIC_THERMAL_TEST_FILE` and optionally
`PB11_DIAGNOSTIC_BEAM_TEST_FILE`. These are qualified local test assets, not
automatically generated production caches. Linux GNU/Clang test linking wraps
the original table destroy functions to observe exactly one matching call.
This proves destruction dispatch, not complete heap leak freedom.

Installed consumers of the exported diagnostics target must resolve
`Threads::Threads` (`find_package(Threads REQUIRED)`) before loading
`pb11Targets.cmake`, in addition to the ordinary core target requirements.

Full trial-input serialization, exclusive snapshot publication, independent
coupled-call replay and the BALDUR Fortran hook are not implemented yet. These
foundations must not be used to claim reconstruction or correction of the
recorded DD numerical failure.
