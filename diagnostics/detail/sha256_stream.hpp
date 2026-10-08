// Incremental SHA-256, FIPS 180-4 sections 4-6.
// Compression arithmetic adapted from pB11_reaction/src/sha256_internal.hpp.
// This diagnostic helper is separate from the physics cache identity closure.
#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace pb11_diagnostics { namespace detail {
class Sha256 {
    using u32 = std::uint32_t;
    using u64 = std::uint64_t;
public:
    void update(const void* bytes, std::size_t size) {
        if (size != 0 && bytes == nullptr)
            throw std::invalid_argument("null SHA256 bytes");
        constexpr u64 max_bytes = std::numeric_limits<u64>::max() / 8u;
        if (size > max_bytes - bytes_)
            throw std::overflow_error("SHA256 bit length overflow");
        if (size == 0) return;
        bytes_ += static_cast<u64>(size);
        const auto* data = static_cast<const unsigned char*>(bytes);
        if (used_ != 0) {
            const std::size_t take = size < 64u - used_ ? size : 64u - used_;
            std::memcpy(buffer_.data() + used_, data, take);
            used_ += take; data += take; size -= take;
            if (used_ == 64u) { compress(buffer_.data()); used_ = 0; }
        }
        while (size >= 64u) {
            compress(data); data += 64u; size -= 64u;
        }
        if (size != 0) {
            std::memcpy(buffer_.data(), data, size);
            used_ = size;
        }
    }

    std::array<unsigned char, 32> digest() const {
        Sha256 copy = *this;
        copy.buffer_[copy.used_++] = 0x80u;
        if (copy.used_ > 56u) {
            std::memset(copy.buffer_.data() + copy.used_, 0, 64u - copy.used_);
            copy.compress(copy.buffer_.data());
            copy.used_ = 0;
        }
        std::memset(copy.buffer_.data() + copy.used_, 0, 56u - copy.used_);
        const u64 bits = copy.bytes_ * 8u;
        for (unsigned i = 0; i < 8; ++i)
            copy.buffer_[56u + i] = static_cast<unsigned char>(bits >> (56u - 8u * i));
        copy.compress(copy.buffer_.data());
        std::array<unsigned char, 32> out{};
        for (unsigned i = 0; i < 8; ++i)
            for (unsigned j = 0; j < 4; ++j)
                out[4u * i + j] = static_cast<unsigned char>(copy.state_[i] >> (24u - 8u * j));
        return out;
    }
private:
    static u32 rotr(u32 x, unsigned n) { return (x >> n) | (x << (32u - n)); }
    static u32 load_be32(const unsigned char* p) {
        return (static_cast<u32>(p[0]) << 24) | (static_cast<u32>(p[1]) << 16) |
               (static_cast<u32>(p[2]) << 8) | static_cast<u32>(p[3]);
    }
    void compress(const unsigned char* block) {
    static constexpr u32 k[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
        0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
        0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
        0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
        0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
        0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
        0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
        0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
    };

        u32 w[64];

        for (unsigned i = 0; i < 16; ++i) {
            w[i] = load_be32(block + 4u * i);
        }

        for (unsigned i = 16; i < 64; ++i) {
            const u32 x = w[i - 15];
            const u32 y = w[i - 2];

            const u32 s0 =
                rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3);

            const u32 s1 =
                rotr(y, 17) ^ rotr(y, 19) ^ (y >> 10);

            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }

        u32 a = state_[0];
        u32 b = state_[1];
        u32 c = state_[2];
        u32 d = state_[3];
        u32 e = state_[4];
        u32 f = state_[5];
        u32 g = state_[6];
        u32 hh = state_[7];

        for (unsigned i = 0; i < 64; ++i) {
            const u32 S1 =
                rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);

            const u32 ch =
                (e & f) ^ ((~e) & g);

            const u32 t1 =
                hh + S1 + ch + k[i] + w[i];

            const u32 S0 =
                rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);

            const u32 maj =
                (a & b) ^ (a & c) ^ (b & c);

            const u32 t2 = S0 + maj;

            hh = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }

        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
        state_[4] += e;
        state_[5] += f;
        state_[6] += g;
        state_[7] += hh;
    }
    std::array<u32, 8> state_{{
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u}};
    std::array<unsigned char, 64> buffer_{};
    u64 bytes_ = 0;
    std::size_t used_ = 0;
};
}} // namespace pb11_diagnostics::detail
