// SHA-256, NIST FIPS 180-4 sections 4-6.
// https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.180-4.pdf
// Private helper; caller guarantees valid bytes and <=512 MiB.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace fusion_detail {

inline std::array<unsigned char, 32>
sha256(const unsigned char* data, std::size_t size)
{
    using u32 = std::uint32_t;
    using u64 = std::uint64_t;

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

    auto rotr = [](u32 x, unsigned n) -> u32 {
        return (x >> n) | (x << (32u - n));
    };

    auto load_be32 = [](const unsigned char* p) -> u32 {
        return (static_cast<u32>(p[0]) << 24) |
               (static_cast<u32>(p[1]) << 16) |
               (static_cast<u32>(p[2]) <<  8) |
                static_cast<u32>(p[3]);
    };

    u32 h[8] = {
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
    };

    auto compress = [&](const unsigned char block[64]) {
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

        u32 a = h[0];
        u32 b = h[1];
        u32 c = h[2];
        u32 d = h[3];
        u32 e = h[4];
        u32 f = h[5];
        u32 g = h[6];
        u32 hh = h[7];

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

        h[0] += a;
        h[1] += b;
        h[2] += c;
        h[3] += d;
        h[4] += e;
        h[5] += f;
        h[6] += g;
        h[7] += hh;
    };

    // Defensive contract for this private helper:
    // nullptr is permitted only for an empty message.
    if (size != 0 && data == nullptr) {
        throw std::invalid_argument("null SHA256 bytes");
    }

    const std::size_t full_blocks = size / 64u;

    for (std::size_t i = 0; i < full_blocks; ++i) {
        compress(data + i * 64u);
    }

    const std::size_t rem = size % 64u;

    // SHA-256 padding needs either one or two final 64-byte blocks.
    unsigned char tail[128] = {};
    if (rem != 0) {
        std::memcpy(tail, data + full_blocks * 64u, rem);
    }

    tail[rem] = 0x80u;

    const std::size_t tail_size = (rem <= 55u) ? 64u : 128u;

    // FIPS 180-4 requires a message length below 2^64 bits.
    // The project-level <=512 MiB limit is far below overflow here.
    const u64 bit_length = static_cast<u64>(size) * 8u;

    unsigned char* len = tail + tail_size - 8u;
    len[0] = static_cast<unsigned char>(bit_length >> 56);
    len[1] = static_cast<unsigned char>(bit_length >> 48);
    len[2] = static_cast<unsigned char>(bit_length >> 40);
    len[3] = static_cast<unsigned char>(bit_length >> 32);
    len[4] = static_cast<unsigned char>(bit_length >> 24);
    len[5] = static_cast<unsigned char>(bit_length >> 16);
    len[6] = static_cast<unsigned char>(bit_length >>  8);
    len[7] = static_cast<unsigned char>(bit_length);

    compress(tail);
    if (tail_size == 128u) {
        compress(tail + 64u);
    }

    std::array<unsigned char, 32> out{};

    for (unsigned i = 0; i < 8; ++i) {
        out[4u * i + 0u] = static_cast<unsigned char>(h[i] >> 24);
        out[4u * i + 1u] = static_cast<unsigned char>(h[i] >> 16);
        out[4u * i + 2u] = static_cast<unsigned char>(h[i] >>  8);
        out[4u * i + 3u] = static_cast<unsigned char>(h[i]);
    }

    return out;
}

} // namespace fusion_detail
