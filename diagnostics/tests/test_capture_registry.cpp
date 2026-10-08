// Real packed-table lifecycle diagnostic. Optional GNU ld destruction wrappers
// observe calls to the actual library destructors, not full heap leak freedom.
// No physics gate.
#include "fusion_capture_registry.h"
#include "sha256_internal.hpp"
#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
static size_t thermal_destroy_calls = 0, beam_destroy_calls = 0;
extern "C" void __real_fusion_c_birth_table_destroy(fusion_birth_table_v1*);
extern "C" void __real_fusion_c_beam_birth_table_destroy(fusion_beam_birth_table_v1*);
extern "C" void __wrap_fusion_c_birth_table_destroy(fusion_birth_table_v1* handle) {
    ++thermal_destroy_calls;
    __real_fusion_c_birth_table_destroy(handle);
}
extern "C" void __wrap_fusion_c_beam_birth_table_destroy(fusion_beam_birth_table_v1* handle) {
    ++beam_destroy_calls;
    __real_fusion_c_beam_birth_table_destroy(handle);
}
#endif

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
using Context = std::unique_ptr<fusion_capture_context_v1,
    decltype(&fusion_capture_context_destroy_v1)>;
static Context context(uint64_t capacity) {
    fusion_capture_context_v1* out = nullptr;
    require(fusion_capture_context_create_v1(capacity, &out) == FUSION_CAPTURE_OK && out,
            "create context");
    return Context(out, fusion_capture_context_destroy_v1);
}
struct Input { std::string path; std::vector<unsigned char> bytes; std::string digest; };
static std::string hex(const std::array<unsigned char, 32>& digest) {
    constexpr char digits[] = "0123456789abcdef";
    std::string out;
    for (auto b : digest) { out += digits[b >> 4]; out += digits[b & 15]; }
    return out;
}
static Input load(const char* name, size_t limit) {
    Input in;
    in.path = std::filesystem::absolute(name).lexically_normal().string();
    require(in.path.size() <= 4096, "input path too long");
    std::ifstream file(in.path, std::ios::binary | std::ios::ate);
    require(bool(file), "open packed file");
    const auto end = file.tellg();
    require(end > 0 && static_cast<uint64_t>(end) <= limit, "packed file size bound");
    in.bytes.resize(static_cast<size_t>(end));
    file.seekg(0);
    require(bool(file.read(reinterpret_cast<char*>(in.bytes.data()), in.bytes.size())),
            "read packed file exactly once");
    in.digest = hex(fusion_detail::sha256(in.bytes.data(), in.bytes.size()));
    return in;
}
static void unchanged(const Input& in) {
    require(hex(fusion_detail::sha256(in.bytes.data(), in.bytes.size())) == in.digest,
            "loader bytes mutated");
}
static bool cleared(const fusion_capture_table_identity_v1& out) {
    const auto* bytes = reinterpret_cast<const unsigned char*>(&out);
    for (size_t i = 0; i < sizeof(out); ++i) if (bytes[i]) return false;
    return true;
}
static void identity(fusion_capture_context_v1* ctx, int kind, const void* handle,
                     const Input& in) {
    fusion_capture_table_identity_v1 out;
    std::memset(&out, 0xa5, sizeof(out));
    require(fusion_capture_table_identity_v1_get(ctx, kind, handle, &out) == FUSION_CAPTURE_OK,
            "identity lookup");
    require(out.kind == kind && out.packed_bytes == in.bytes.size(), "kind and packed length");
    require(out.content_sha256[64] == 0 && in.digest == out.content_sha256,
            "existing SHA256 reference matches");
    require(out.kernel_identity[64] == 0 &&
            std::string(out.kernel_identity) == fusion_c_beam_birth_table_kernel_identity(),
            "kernel identity matches actual library");
    require(std::string(out.qualified_path) == in.path, "qualified path matches");
}
static void missing(fusion_capture_context_v1* ctx, int kind, const void* handle) {
    fusion_capture_table_identity_v1 out;
    std::memset(&out, 0xa5, sizeof(out));
    require(fusion_capture_table_identity_v1_get(ctx, kind, handle, &out) ==
            FUSION_CAPTURE_UNKNOWN_HANDLE && cleared(out), "unknown/wrong kind clears identity");
}
static fusion_birth_table_v1* thermal(fusion_capture_context_v1* ctx, const Input& in) {
    fusion_birth_table_v1* handle = nullptr;
    int original = -9;
    require(fusion_capture_thermal_unpack_v1(ctx, in.bytes.data(), in.bytes.size(),
            in.path.c_str(), &handle, &original) == FUSION_CAPTURE_OK && original == 0 && handle,
            "real thermal unpack publishes with original success");
    return handle;
}
static void thermal_tests(const Input& in) {
    fusion_capture_context_v1* invalid = reinterpret_cast<fusion_capture_context_v1*>(uintptr_t(1));
    require(fusion_capture_context_create_v1(0, &invalid) == FUSION_CAPTURE_INVALID && !invalid,
            "zero capacity clears context output");
    invalid = reinterpret_cast<fusion_capture_context_v1*>(uintptr_t(1));
    require(fusion_capture_context_create_v1(1000001, &invalid) == FUSION_CAPTURE_INVALID && !invalid,
            "excess capacity clears context output");
    require(fusion_capture_context_create_v1(1, nullptr) == FUSION_CAPTURE_INVALID,
            "null context output rejected");
    auto ctx = context(1);
    auto* first = thermal(ctx.get(), in);
    identity(ctx.get(), FUSION_CAPTURE_THERMAL_TABLE, first, in);
    missing(ctx.get(), FUSION_CAPTURE_BEAM_TABLE, first);
    int unknown_token = 0;
    missing(ctx.get(), FUSION_CAPTURE_THERMAL_TABLE, &unknown_token);
    auto* second = first;
    int original = 55;
    require(fusion_capture_thermal_unpack_v1(ctx.get(), in.bytes.data(), in.bytes.size(),
            in.path.c_str(), &second, &original) == FUSION_CAPTURE_CAPACITY && !second && original == -1,
            "capacity rejects second unpack before original call");
    identity(ctx.get(), FUSION_CAPTURE_THERMAL_TABLE, first, in);
    auto* wrong = reinterpret_cast<fusion_beam_birth_table_v1*>(first);
    require(fusion_capture_beam_destroy_v1(ctx.get(), &wrong) == FUSION_CAPTURE_UNKNOWN_HANDLE &&
            reinterpret_cast<void*>(wrong) == first, "wrong destroy preserves handle");
    auto* unknown = reinterpret_cast<fusion_birth_table_v1*>(&unknown_token);
    require(fusion_capture_thermal_destroy_v1(ctx.get(), &unknown) == FUSION_CAPTURE_UNKNOWN_HANDLE &&
            reinterpret_cast<void*>(unknown) == &unknown_token, "unknown destroy preserves handle");
    const void* removed = first;
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    size_t thermal_before = thermal_destroy_calls, beam_before = beam_destroy_calls;
#endif
    require(fusion_capture_thermal_destroy_v1(ctx.get(), &first) == FUSION_CAPTURE_OK && !first,
            "correct destroy clears handle");
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    require(thermal_destroy_calls == thermal_before + 1 && beam_destroy_calls == beam_before,
            "explicit thermal release calls matching original destructor exactly once");
#endif
    missing(ctx.get(), FUSION_CAPTURE_THERMAL_TABLE, removed);
    // Corrupt only a copy's magic byte; preserve the caller's original buffer.
    auto corrupt = in.bytes;
    corrupt[0] ^= 0xff;
    second = nullptr; original = 55;
    require(fusion_capture_thermal_unpack_v1(ctx.get(), corrupt.data(), corrupt.size(),
            in.path.c_str(), &second, &original) == FUSION_CAPTURE_UNPACK_REJECTED &&
            !second && original != 0 && original != 55,
            "corrupt bytes return separate diagnostic 104 and actual original rejection");
    first = thermal(ctx.get(), in); // Failed unpack did not consume capacity.
    identity(ctx.get(), FUSION_CAPTURE_THERMAL_TABLE, first, in);
    removed = first;
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    thermal_before = thermal_destroy_calls; beam_before = beam_destroy_calls;
#endif
    ctx.reset(); // Deliberately leave one owned thermal table.
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    require(thermal_destroy_calls == thermal_before + 1 && beam_destroy_calls == beam_before,
            "context calls matching original thermal destructor exactly once");
#endif
    auto fresh = context(1);
    missing(fresh.get(), FUSION_CAPTURE_THERMAL_TABLE, removed);
    second = reinterpret_cast<fusion_birth_table_v1*>(uintptr_t(1)); original = 55;
    require(fusion_capture_thermal_unpack_v1(fresh.get(), in.bytes.data(), in.bytes.size(),
            "relative/path", &second, &original) == FUSION_CAPTURE_INVALID && !second && original == -1,
            "invalid unpack clears outputs");
    unchanged(in);
}
static void beam_tests(const Input& in) {
    auto ctx = context(1);
    fusion_beam_birth_table_v1* handle = nullptr;
    int original = -9;
    require(fusion_capture_beam_unpack_v1(ctx.get(), in.bytes.data(), in.bytes.size(),
            in.path.c_str(), &handle, &original) == FUSION_CAPTURE_OK && handle && original == 0,
            "real beam unpack publishes with original success");
    identity(ctx.get(), FUSION_CAPTURE_BEAM_TABLE, handle, in);
    missing(ctx.get(), FUSION_CAPTURE_THERMAL_TABLE, handle);
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    size_t thermal_before = thermal_destroy_calls, beam_before = beam_destroy_calls;
#endif
    require(fusion_capture_beam_destroy_v1(ctx.get(), &handle) == FUSION_CAPTURE_OK && !handle,
            "beam release clears handle");
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    require(beam_destroy_calls == beam_before + 1 && thermal_destroy_calls == thermal_before,
            "explicit beam release calls matching original destructor exactly once");
#endif
    require(fusion_capture_beam_unpack_v1(ctx.get(), in.bytes.data(), in.bytes.size(),
            in.path.c_str(), &handle, &original) == FUSION_CAPTURE_OK && handle && original == 0,
            "beam capacity reusable after release");
    const void* stale = handle;
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    thermal_before = thermal_destroy_calls; beam_before = beam_destroy_calls;
#endif
    ctx.reset(); // Deliberately leave one owned beam table.
#ifdef PB11_DIAGNOSTIC_DESTROY_OBSERVATION
    require(beam_destroy_calls == beam_before + 1 && thermal_destroy_calls == thermal_before,
            "context calls matching original beam destructor exactly once");
#endif
    auto fresh = context(1);
    missing(fresh.get(), FUSION_CAPTURE_BEAM_TABLE, stale);
    unchanged(in);
}
int main(int argc, char** argv) {
    try {
        require(argc == 2 || argc == 3, "usage: test_capture_registry thermal-packed-file [beam-packed-file]");
        const Input thermal_input = load(argv[1], 512ULL * 1024 * 1024);
        thermal_tests(thermal_input);
        if (argc == 3) {
            const Input beam_input = load(argv[2], 256ULL * 1024 * 1024);
            beam_tests(beam_input);
        }
        std::cout << "capture registry lifecycle diagnostic passed; physics acceptance not assessed\n";
    } catch (const std::exception& error) {
        std::cerr << "capture registry test: " << error.what() << '\n'; return 1;
    }
    return 0;
}
