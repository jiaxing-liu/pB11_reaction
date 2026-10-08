// Standalone diagnostic test; compile with -I/path/to/pB11_reaction/src.
#include "../detail/sha256_stream.hpp"
#include "sha256_internal.hpp"
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using pb11_diagnostics::detail::Sha256;
using Digest = std::array<unsigned char, 32>;

static void require(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(what);
}
static std::string hex(const Digest& bytes) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string out;
    for (auto byte : bytes) {
        out += digits[byte >> 4]; out += digits[byte & 15u];
    }
    return out;
}
static Digest hash(const void* bytes, std::size_t size) {
    Sha256 sha; sha.update(bytes, size); return sha.digest();
}
static void vectors() {
    require(hex(hash(nullptr, 0)) ==
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "empty vector");
    require(hex(hash("abc", 3)) ==
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "abc vector");
    Sha256 million;
    std::array<unsigned char, 1000> chunk{}; chunk.fill('a');
    for (unsigned i = 0; i < 1000; ++i) million.update(chunk.data(), chunk.size());
    require(hex(million.digest()) ==
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0", "million-a vector");
}
static void splits_and_binary() {
    for (std::size_t size : {55u, 56u, 63u, 64u, 65u, 127u, 128u, 129u}) {
        std::vector<unsigned char> bytes(size);
        for (std::size_t i = 0; i < size; ++i)
            bytes[i] = static_cast<unsigned char>((i * 197u) & 255u);
        const Digest expected = fusion_detail::sha256(bytes.data(), bytes.size());
        require(hash(bytes.data(), bytes.size()) == expected, "binary whole buffer");
        for (std::size_t split = 0; split <= size; ++split) {
            Sha256 sha;
            sha.update(bytes.data(), split);
            sha.update(nullptr, 0);
            sha.update(bytes.data() + split, size - split);
            require(sha.digest() == expected, "every binary split");
        }
        Sha256 singles;
        for (auto byte : bytes) singles.update(&byte, 1);
        require(singles.digest() == expected, "single-byte updates");
    }
    std::array<unsigned char, 256> all{};
    for (unsigned i = 0; i < all.size(); ++i) all[i] = static_cast<unsigned char>(i);
    require(hash(all.data(), all.size()) == fusion_detail::sha256(all.data(), all.size()),
            "all byte values including zero and high bits");
}
static void independent_state_and_contracts() {
    Sha256 sha; sha.update("a", 1);
    const auto first = sha.digest();
    require(first == hash("a", 1) && sha.digest() == first, "repeat digest preserves state");
    Sha256 copy = sha;
    sha.update("bc", 2); copy.update("xy", 2);
    require(sha.digest() == hash("abc", 3), "digest then update");
    require(copy.digest() == hash("axy", 3), "copy independently updates");
    const auto before = sha.digest();
    bool rejected = false;
    try { sha.update(nullptr, 1); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected && sha.digest() == before, "nonzero null rejected without mutation");
    sha.update(nullptr, 0);
    require(sha.digest() == before, "zero null permitted");
    // Oversize input must be rejected before dereferencing the valid tiny pointer.
    if (std::numeric_limits<std::size_t>::max() > std::numeric_limits<std::uint64_t>::max() / 8u) {
        const unsigned char byte = 0;
        rejected = false;
        try { sha.update(&byte, std::numeric_limits<std::size_t>::max()); }
        catch (const std::overflow_error&) { rejected = true; }
        require(rejected && sha.digest() == before, "bit length overflow rejected without mutation");
    }
}
int main() {
    try { vectors(); splits_and_binary(); independent_state_and_contracts(); }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
    std::cout << "SHA256 stream tests passed\n";
    return 0;
}
