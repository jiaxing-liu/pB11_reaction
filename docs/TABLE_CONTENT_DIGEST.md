# Immutable table content identity

Thermal and beam table C APIs now expose `content_digest(table,digest,length)`;
matching ISO_C_BINDING wrappers return32 raw bytes and `integer(c_size_t)`.
The digest is SHA256 of the complete existing canonical pack, including its
format, current kernel identity and checksum. Length is the exact packed size.
This does not change table evaluation or validate a physical model. It does not
perform file I/O, authentication, host restart or persist beam routing identity.

All nonnull outputs are zeroed before validation, including failures caused by
another null output. One temporary packed buffer is bounded by existing limits:
512MiB thermal,256MiB beam. No new OpenSSL/system dependency. The private helper
implements NIST FIPS180-4 sections4–6; this is not a certification claim:
https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.180-4.pdf
Chat supplied the implementation and review; executor changed invalid private
null input to throw, so it cannot silently represent a successful zero digest.

Tests: standard empty/abc/multiblock/million-a answers;15 deterministic Python
hashlib comparisons including55/56/63/64/65-byte padding boundaries; two actual
DT thermal/beam tables matched Python. Pack/unpack identity, unchanged pack after
query, changed table metadata and C failure-output contracts pass. Intel Fortran
default and-r8 digests independently match Python and clear failed outputs.
GNU CMake CTest passes C++, Python and Fortran fixtures. Normal Intel BALDUR build
passes and contains both new ABI symbols. Raw logs/bytes/hashes in
validation/table-digest/. This test does not claim all-channel physics coverage.

Reproduce with a fresh CMake build, build targets test_fusion_table_digest and
(optional Fortran) test_fusion_table_digest_fortran, then `ctest -R
fusion_table_digest --output-on-failure`. Python tests are registered if the
interpreter is found. The C++ fixture generates the files required by Fortran.

Kernel identity intentionally includes all source/header files, so this code
change invalidates prior cache identities even though no physics changed.
Do not rewrite old cache headers to bypass this guard. Existing running host
processes retain their linked library and handles. Full restart must record
capture-time content identity plus beam channel/slot/cell, validate candidate
handles, then retain those same handles through publication (no later reload).
