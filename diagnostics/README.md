# Optional local-trial diagnostics

Enable `PB11_BUILD_DIAGNOSTICS=ON` to build `pb11::diagnostics`. This separate
target implements streaming SHA256, an explicitly owned table registry,
local-trial snapshot capture, and independent replay with C interfaces. It has no reaction, collision or burn formulas.
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
Replay verifies file content, size and kernel identity before importing it. Registration is not interpolation or physical validation.
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

## Trial capture and replay

`fusion_capture_trial.h` defines a borrowed C view of one local covered or
packet trial. A caller supplies its original numerical status and exact input
buffers before they change. Required pointers and dimensions must satisfy the
view contract. The optional thermal pointer is significant: nullptr selects
direct evaluation; an all-null nonnull handle array requests table mode and can
produce an original invalid-argument status for active thermal channels.
All optional pointer-presence distinctions are retained. Capture rejects
unsupported overlapping input spans; it does not reproduce arbitrary invalid
output-pointer/alias layouts. Outputs are reconstructed using valid, distinct
buffers, with the floor-output pointer paired to the recorded floor limits.

The binary file uses versioned, ordered fields, explicit little-endian integer
and IEEE binary64 representations, and a SHA256 checksum. Addresses and C ABI
padding are never written. All controls, reservoir populations and energies,
fast S/T arrays, grid edges, collision logs, inert baths, external birth/escape,
floor limits and table-domain policy are included. Thermal/beam tables are
references to registered packed bytes, not embedded or repacked payloads.
Capture accepts budgets of 1 through 128 MiB including checksum; a bounded
stream writes a same-directory temporary file. On POSIX systems publication
is atomic and exclusive, with mode 0600. Existing files are never overwritten.
The writer reports unsupported on Windows; the wire format is portable.
Registered inputs must remain live and unchanged throughout capture.

`fusion_trial_replay snapshot.pbtrial` verifies record shape/checksum/current
kernel, then imports only exact referenced files (regular nonsymlink, expected
size, digest, kernel), deduplicating repeated table references. There is no
cache construction, fallback policy change, threshold relaxation or host-state
mutation. Missing files report diagnostic IO (106), changed identities report
107, and original numerical status is separate. The JSON output always says
`physics_accepted:false`. Exit 0 means the original/replayed numerical statuses
match, including matched numerical failures; it is not physical acceptance.
The output digest explicitly encodes all output members, not ABI padding.

Snapshot tests compare exact inputs and complete outputs against independent
calls to the original covered and packet APIs. Status 0/3/4 fixture tests are
local arithmetic parity, not the actual failed DD host trial. Real-table
identity tests exercise registry references and file changes with reactions
disabled; they do not claim cache use, interpolation validation or physical
acceptance. Test outputs use exclusive fresh directories; repeated CTest runs
preserve previous evidence. Table assets are optional and explicitly supplied.

The BALDUR capture hook is not implemented yet. No actual DD rejected tuple
has been reconstructed or its cause established by this component alone.
