# Reusable beam-table bytes

The immutable beam source table can be saved and loaded with
`fusion_c_beam_birth_table_pack_size`, `fusion_c_beam_birth_table_pack` and
`fusion_c_beam_birth_table_unpack`. The caller owns buffers and file I/O. No
BALDUR state, mesh coordinates in real space, density, host time advance or
filesystem path is stored by the library.

A loaded table uses the same `info`, `evaluate`, `destroy` and coupled sparse
beam-entry APIs as a newly built table. It retains the complete source
spectrum, including every nonzero subnormal entry, source options, energy
grid, sampled temperature knots and construction diagnostics. It is not a
moment-only replacement or a new interpolation model.

## Compatibility and provenance

`fusion_c_beam_birth_table_kernel_identity()` returns a process-lifetime
64-character lowercase SHA256 string. CMake generates it from sorted relative
paths and file SHA256 values for `.cpp`, `.h`, `.hpp` and `.inc` files under
`src`/`include`, plus the Boost version. Edits and additions/removals trigger
regeneration on the next CMake build. Absolute paths, git state and compiler
flags are excluded. This is a conservative source/data compatibility key:
unrelated source edits may invalidate a cache. It is not binary provenance.
New source/data file types outside this list must be added to the identity
inputs before their values can affect the beam kernel.

Import rejects a different identity. Moving an unchanged source tree does not
change it. The checksum detects accidental byte corruption; it is not an
authentication mechanism. Keep caches with trusted build/input provenance.
Import validates structure and per-knot conservation but does not recompute
the direct quadratures, rediscover sampled error envelopes or prove unsampled
accuracy. Stored envelopes remain the original sampled evidence. Independent
spectral/debit checks and coupled convergence are still required.

## Version1 wire format

All integers are unsigned64 little-endian; signed API option fields must be
nonnegative and fit their declared integer types. All real values use the
IEEE754 binary64 bit pattern in little-endian order. No C structure padding,
pointers or native Fortran storage is serialized.

1. Magic `0x3142544e4f495346`, format version1, total byte length.
2. The64 ASCII bytes of kernel identity, without its NUL terminator.
3. The42 metadata words: the nine real fields and seven integer fields of
   `fusion_beam_birth_table_info_v1`; then source options (eight real fields,
   nine integer fields); then controls (six real fields, three integer
   fields), each group in declaration order in the public headers.
4. `cells+1` energy-grid edges.
5. Each knot: temperature,31 source-coefficient real fields in declaration
   order, number of sparse entries, then `(index,value)` pairs. Indices are
   in the dense seven-species-major grid and strictly increasing within a
   knot. Every stored value is strictly positive, including subnormals.
6. Trailing FNV-1a64 checksum over every preceding byte, using offset basis
   `14695981039346656037` and multiplier `1099511628211`.

Length is `432 + 8*(cells+1) + 264*knots + 16*stored_entries` bytes.
No trailing bytes are accepted. The format is bounded by256MiB,100,000cells,
100,000knots and12,500,000 cumulative sampled spectral entries; stored entries
cannot exceed that sampled count. Individual sparse knots have at most
`7*cells` entries. These are validation/resource bounds, not accuracy targets.

## Failure and lifetime behavior

Size/written outputs clear to zero on failure. A short pack buffer is left
untouched. Unpack clears the output handle before checking input and publishes
only after complete validation. It bounds lengths/counts before allocations,
checks identity/checksum, source/control domains, ordered grid/temperatures,
sampled-envelope consistency, finite nonnegative coefficients and sparse
indices/values. Every knot must satisfy the same particle/energy balance
criterion as table construction; sparse validation visits stored entries
without allocating a dense grid per knot.

The importer does not impose an extra relative projectile-moment constraint:
that is not an existing constructor invariant when a tiny moment rounds to
zero. Valid serialized tables must survive exact round trips, including
representable subnormal source values.

The caller must destroy each successful handle once. Importing into a variable
that already owns a table does not destroy the old table automatically. Arrays
must not overlap the table, each other or output control values.

Build through CMake so the private generated identity header is available.

Fortran exposes `fusion_beam_birth_table_pack_size/pack/unpack` with explicit
`c_size_t` lengths and caller-owned contiguous `integer(c_int8_t)` byte arrays.
Pack uses `intent(inout)` to preserve buffers on early rejection. Wrappers
reject lengths beyond the actual array before forming a C pointer.


## Exact request matching for caller-owned caches

`fusion_c_beam_birth_table_matches_request` compares every constructor input:
channel, projectile slot/energy, full energy grid, temperature interval, named
source options and all accuracy/resource controls. It returns OK with matches0
for a different request, including invalid numeric request values; null inputs
and invalid cell extents are argument errors. NaNs never match a valid table.
Equality is numerical, without a tolerance; signed zeros compare equal.
The Fortran wrapper derives cells from the supplied edge array and returns a
logical result. Neither interface integrates, mutates nor selects a fallback.

A successful unpack alone does not establish suitability for another request.
Callers should require a matching request before reuse. Import still verifies
the current kernel identity; a source change may intentionally invalidate an
older cache even when the reaction equations themselves are unchanged.
