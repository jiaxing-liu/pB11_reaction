# Validation provenance

Coordinator report and evidence archive are in the paired BALDUR checkout at
`docs/exl50u-program/THERMAL_TABLE_CACHE.md` and `thermal-table-cache/`.
The design named below resides in that coordinator directory.

# D146: reusable thermal-table persistence

Host parent5068990; external parent e0b7881. New revisions are pending final
archive verification/paired commits. This increment implements the immutable
thermal-table cache designed in THERMAL_TABLE_CACHE_DESIGN.md.

## Library and host behavior

Additive C and Fortran APIs expose kernel identity, exact complete constructor
request matching and bounded pack-size/pack/unpack. Canonical little-endian
binary64/integer encoding preserves all coefficient and dense spectrum bits.
The format is distinct from beam tables, with exact lengths,512MiB cap, original
constructor storage budget, shared source/data fingerprint and FNV1a checksum.
The checksum detects corruption; it is not authenticity or proof of sampled
accuracy. Import validates metadata, source/control domains, knot/edge ordering,
per-node stoichiometry and energy, and uses the existing thermal parent-domain
gate (including pB12MeV limit). It performs no quadrature or filesystem access.
A shared internal preflight helper preserves the original direct-source parent
construction and checks, avoiding a second implementation of that model bound.

Fortran checks array extents before C_LOC, clears failed result counts/handles
and keeps caller-owned lifetime. BALDUR probe gains optional --thermal-cache-dir.
Host stream I/O creates missing files with status=new, never overwrites an
existing cache, and requires complete request equality after unpack. Each
channel records loaded flag, knots/direct-evaluation metadata and preparation
CPU seconds. The selected core libpb11.a now comes from the same --ifx-build as
the selected Fortran bindings, preventing stale/mixed static-core selection.
The host owns file paths; the independent C++ library remains host-neutral.

## Verification and preserved failures

GNU and Intel-backed Fortran tests pass actual DT table roundtrip, exact knot
and off-knot spectra, exact repack, full request matching and invalid buffers.
Checksum-repaired corruptions exercise semantic checks rather than only the
checksum: model identity/version/length, dimensions/budgets, nonfinite metadata,
invalid source domains, repeated edges, missing endpoints and negative spectra.
A real narrow-domain pB table roundtrip and a source-bound mutation check the pB
path. That tiny fixture deliberately uses coarse quadrature and loose gates for
serialization coverage; it is not physical-model convergence evidence.
A synthetic minimum-subnormal alpha-cell perturbation preserves its exact bit
pattern through pack/unpack and exact-knot evaluation. No tail clipping added.
Original thermal table and beam-cache C++ tests and original thermal/beam-byte
Fortran bindings pass. GNU host file I/O tests verify exact load, missing file,
rejected overwrite with unchanged bytes, request mismatch and corrupt file.

First root C++ test collided with the original fixture's TableOwner name;
renamed BytesTableOwner. First Fortran edge mutation1e-30 rounded away; NEAREST
now makes an actual one-ULP mismatch. Missing TARGET on C_LOC(answer) was fixed
in source review. First GNU compile-only attempt exhausted /tmp; subsequent
builds use existing workspace TMPDIR. Failed logs are retained. The first
preflight edit's source anchor assertion failed before source modification;
its successful build did not test the later helper. Corrected build and tests
are identified separately. Final ASan+UBSan3131 configure/build/runtime all passed after the helper
and pB/subnormal additions; the earlier sanitizer success is preserved separately.

## Actual BALDUR cache parity

The controlled DT network run requests100us and accepts9steps through110us,
with moving geometry and handoff enabled. DD/DT/DHe3 thermal-channel tables are
prepared; this is not the final matched0.6nG four-fuel operating comparison.

Creation46443 and cached83523 are both terminal0 and complete the requested
window. Seven physics CSVs (channels, species, nuclear history, host energy,
plasma history/profiles and handoff diagnostics) plus jobxdat are byte-identical.
All four channel tables report create0 then load1; knot/evaluation metadata are
unchanged. Direct-evaluation counts in a loaded table describe its original
construction, not fresh evaluations during loading.

| Measurement | Create/save | Load |
|---|---:|---:|
| Total table preparation CPU seconds |155.443405|0.090315|
| Whole probe wall seconds |159.31169045|3.87053895|

Preparation includes I/O and validation. Wall time includes initialization and
integration, excludes compilation/calibration. These timings concern this
short controlled run under concurrent workloads; they do not predict pB
long-discharge speedup. parity.json records hashes, scope and exact file list.

## Remaining full-program scope

This table cache neither supplies missing transport/fast-loss closures nor
completes restart support or production configuration. Loaded tables retain
sampled interpolation guarantees only, and changing source/data revision can
invalidate previous caches. Keep pinned live builds; never bypass identity
checks. Full pB/DD windows, all twelve main cases, physical convergence and
final reporting remain required. No new ignition or device-performance claim.
