#include "util/simd_dispatch.hpp"

#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <new>
#include <string>

#if defined(__x86_64__) || defined(__i386__)
  #include <cpuid.h>
  #define IKAFSSN_ARCH_X86 1
#endif

#if defined(__aarch64__)
  #include <sys/auxv.h>
  #include <asm/hwcap.h>
  #define IKAFSSN_ARCH_ARM 1
#endif

#ifndef IKAFSSN_ENABLE_SIMD
  #define IKAFSSN_ENABLE_SIMD 1
#endif

namespace ikafssn {

namespace {

// Persistent state populated by do_init_once() and queried via the public
// accessors. Stored in plain atomics because writes happen exactly once.
std::once_flag       g_once;
std::atomic<SimdCap> g_auto_cap   { SimdCap::Scalar };
std::atomic<SimdCap> g_current    { SimdCap::Scalar };
std::atomic<bool>    g_slow_bmi2  { false };
std::atomic<bool>    g_avx512f_only { false };
std::atomic<bool>    g_initialized { false };

// Normalize a force-simd token: lowercase, strip '-', '_', '.', whitespace.
std::string normalize_token(const char* raw) {
    std::string out;
    if (!raw) return out;
    out.reserve(std::strlen(raw));
    for (const char* p = raw; *p; ++p) {
        unsigned char c = static_cast<unsigned char>(*p);
        if (c == '-' || c == '_' || c == '.' || std::isspace(c)) continue;
        out.push_back(static_cast<char>(std::tolower(c)));
    }
    return out;
}

#if IKAFSSN_ARCH_X86
struct X86Info {
    char vendor[13]      = {0};
    unsigned family      = 0;
    unsigned model       = 0;
    bool     has_avx2    = false;
    bool     has_bmi2    = false;
    bool     has_avx512f = false;
    bool     has_avx512bw = false;
    bool     has_avx512vbmi  = false;
    bool     has_avx512vbmi2 = false;
    bool     has_sse42   = false;
};

X86Info detect_x86() noexcept {
    X86Info info;

    __builtin_cpu_init();

    // Vendor + family/model via raw CPUID for AMD Zen1/Zen2 detection.
    unsigned eax = 0, ebx = 0, ecx = 0, edx = 0;
    if (__get_cpuid(0u, &eax, &ebx, &ecx, &edx)) {
        std::memcpy(info.vendor + 0, &ebx, 4);
        std::memcpy(info.vendor + 4, &edx, 4);
        std::memcpy(info.vendor + 8, &ecx, 4);
        info.vendor[12] = '\0';
    }
    if (__get_cpuid(1u, &eax, &ebx, &ecx, &edx)) {
        unsigned base_family = (eax >> 8) & 0xF;
        unsigned ext_family  = (eax >> 20) & 0xFF;
        unsigned base_model  = (eax >> 4) & 0xF;
        unsigned ext_model   = (eax >> 16) & 0xF;
        info.family = (base_family == 0xF) ? (base_family + ext_family) : base_family;
        info.model  = (base_family == 0x6 || base_family == 0xF)
                        ? ((ext_model << 4) | base_model)
                        : base_model;
    }

    info.has_sse42        = __builtin_cpu_supports("sse4.2");
    info.has_avx2         = __builtin_cpu_supports("avx2");
    info.has_bmi2         = __builtin_cpu_supports("bmi2");
    info.has_avx512f      = __builtin_cpu_supports("avx512f");
    info.has_avx512bw     = __builtin_cpu_supports("avx512bw");
    info.has_avx512vbmi   = __builtin_cpu_supports("avx512vbmi");
    info.has_avx512vbmi2  = __builtin_cpu_supports("avx512vbmi2");

    return info;
}
#endif // IKAFSSN_ARCH_X86

void do_init_once(Logger* logger) noexcept {
    SimdCap      auto_cap = SimdCap::Scalar;
    SimdAuxFlags flags;
    bool         build_disabled = (IKAFSSN_ENABLE_SIMD == 0);

#if IKAFSSN_ARCH_X86
    X86Info x = build_disabled ? X86Info{} : detect_x86();
    if (!build_disabled) {
        // AVX-512 BW is the minimum AVX-512 tier we support. AVX-512F-only
        // (KNL/KNM) is demoted to AVX2.
        if (x.has_avx512vbmi2 && x.has_avx512vbmi && x.has_avx512bw && x.has_avx512f) {
            auto_cap = SimdCap::AVX512VBMI2;
        } else if (x.has_avx512vbmi && x.has_avx512bw && x.has_avx512f) {
            auto_cap = SimdCap::AVX512VBMI;
        } else if (x.has_avx512bw && x.has_avx512f) {
            auto_cap = SimdCap::AVX512BW;
        } else if (x.has_avx2) {
            auto_cap = SimdCap::AVX2;
            if (x.has_avx512f) flags.has_avx512f_only = true;
        } else if (x.has_sse42) {
            auto_cap = SimdCap::SSE42;
        } else {
            auto_cap = SimdCap::Scalar;
        }

        // AMD Zen1/Zen2 (and older Bulldozer/Jaguar): pdep/pext are microcoded
        // and ~10x slower than scalar. Mark the BMI2-slow auxiliary flag so
        // callers can opt out of those code paths while keeping the AVX2 tier.
        if (std::strcmp(x.vendor, "AuthenticAMD") == 0) {
            if (x.family == 0x15 || x.family == 0x16 || x.family == 0x17) {
                flags.slow_bmi2 = true;
            }
        }
    }
#elif IKAFSSN_ARCH_ARM
    if (!build_disabled) {
        unsigned long hwcap  = getauxval(AT_HWCAP);
        unsigned long hwcap2 = getauxval(AT_HWCAP2);
        bool has_neon = (hwcap & HWCAP_ASIMD) != 0;
        bool has_sve  = (hwcap & HWCAP_SVE)   != 0;
        bool has_sve2 = false;
        bool has_sme  = false;
        bool has_sme2 = false;
        #ifdef HWCAP2_SVE2
            has_sve2 = (hwcap2 & HWCAP2_SVE2) != 0;
        #endif
        #ifdef HWCAP2_SME
            has_sme  = (hwcap2 & HWCAP2_SME)  != 0;
        #endif
        #ifdef HWCAP2_SME2
            has_sme2 = (hwcap2 & HWCAP2_SME2) != 0;
        #endif
        if (has_sme2)      auto_cap = SimdCap::SME2;
        else if (has_sme)  auto_cap = SimdCap::SME;
        else if (has_sve2) auto_cap = SimdCap::SVE2;
        else if (has_sve)  auto_cap = SimdCap::SVE;
        else if (has_neon) auto_cap = SimdCap::NEON;
        else               auto_cap = SimdCap::Scalar;
        // SME / SME2 are detected only; the present kernels do not use them.
        // Cap the auto-tier at SVE2 to avoid surprising user-visible labels.
        if (auto_cap == SimdCap::SME || auto_cap == SimdCap::SME2) {
            auto_cap = has_sve2 ? SimdCap::SVE2 : (has_sve ? SimdCap::SVE
                                                            : SimdCap::NEON);
        }
        (void)has_neon;
    }
#else
    (void)build_disabled;
#endif

    g_auto_cap.store(auto_cap, std::memory_order_release);
    g_slow_bmi2.store(flags.slow_bmi2, std::memory_order_release);
    g_avx512f_only.store(flags.has_avx512f_only, std::memory_order_release);

    // IKAFSSN_FORCE_SIMD honors only downgrade requests; over-requests collapse
    // to Scalar with a WARN.
    SimdCap current = auto_cap;
    const char* env = std::getenv("IKAFSSN_FORCE_SIMD");
    bool forced = false;
    bool force_warned = false;

    if (build_disabled) {
        current = SimdCap::Scalar;
        if (logger) logger->info("simd: scalar (build-time disabled)");
    } else if (env && *env) {
        ParseForceResult parsed = parse_force_simd_env(env);
        if (!parsed.explicit_value) {
            // unrecognized token
            if (logger) logger->warn(
                "simd: IKAFSSN_FORCE_SIMD='%s' not recognized; using scalar",
                env);
            current = SimdCap::Scalar;
            forced = true;
            force_warned = true;
        } else if (static_cast<int>(parsed.cap) >
                   static_cast<int>(auto_cap)) {
            // over-request rejected
            if (logger) logger->warn(
                "simd: IKAFSSN_FORCE_SIMD=%s exceeds detected capability %s; "
                "using scalar",
                env, simd_cap_name(auto_cap).data());
            current = SimdCap::Scalar;
            forced = true;
            force_warned = true;
        } else {
            current = parsed.cap;
            forced  = true;
        }
    }

    g_current.store(current, std::memory_order_release);
    g_initialized.store(true, std::memory_order_release);

    if (!build_disabled && logger) {
        if (forced && !force_warned) {
            logger->info(
                "simd: %s (forced via IKAFSSN_FORCE_SIMD=%s; auto-detected: %s)",
                simd_cap_name(current).data(),
                env ? env : "",
                simd_cap_name(auto_cap).data());
        } else if (!forced) {
            logger->info("simd: %s (auto)", simd_cap_name(current).data());
        }

#if IKAFSSN_ARCH_X86
        if (logger->verbose()) {
            logger->debug("simd cpu: vendor=%s family=0x%x model=0x%x "
                          "sse42=%d avx2=%d bmi2=%d avx512f=%d "
                          "avx512bw=%d avx512vbmi=%d avx512vbmi2=%d "
                          "slow_bmi2=%d avx512f_only=%d",
                          x.vendor, x.family, x.model,
                          (int)x.has_sse42, (int)x.has_avx2,
                          (int)x.has_bmi2, (int)x.has_avx512f,
                          (int)x.has_avx512bw, (int)x.has_avx512vbmi,
                          (int)x.has_avx512vbmi2,
                          (int)flags.slow_bmi2, (int)flags.has_avx512f_only);
        }
#endif
    }
}

} // namespace

void init_simd_dispatch(Logger* logger) noexcept {
    std::call_once(g_once, [&] { do_init_once(logger); });
}

SimdCap current_simd_cap() noexcept {
    return g_current.load(std::memory_order_acquire);
}

SimdAuxFlags current_simd_flags() noexcept {
    SimdAuxFlags flags;
    flags.slow_bmi2        = g_slow_bmi2.load(std::memory_order_acquire);
    flags.has_avx512f_only = g_avx512f_only.load(std::memory_order_acquire);
    return flags;
}

SimdCap auto_detected_simd_cap() noexcept {
    return g_auto_cap.load(std::memory_order_acquire);
}

std::string_view simd_cap_name(SimdCap cap) noexcept {
    switch (cap) {
        case SimdCap::Scalar:       return "scalar";
        case SimdCap::SSE42:        return "sse42";
        case SimdCap::AVX2:         return "avx2";
        case SimdCap::AVX512BW:     return "avx512bw";
        case SimdCap::AVX512VBMI:   return "avx512vbmi";
        case SimdCap::AVX512VBMI2:  return "avx512vbmi2";
        case SimdCap::NEON:         return "neon";
        case SimdCap::SVE:          return "sve";
        case SimdCap::SVE2:         return "sve2";
        case SimdCap::SME:          return "sme";
        case SimdCap::SME2:         return "sme2";
    }
    return "unknown";
}

SimdCap force_simd_cap(SimdCap requested) noexcept {
    SimdCap auto_cap = auto_detected_simd_cap();
    SimdCap applied  = (static_cast<int>(requested) <= static_cast<int>(auto_cap))
                         ? requested
                         : SimdCap::Scalar;
    g_current.store(applied, std::memory_order_release);
    return applied;
}

ParseForceResult parse_force_simd_env(const char* env_value) noexcept {
    ParseForceResult result{ SimdCap::Scalar, false, "" };
    if (!env_value) return result;
    std::string tok = normalize_token(env_value);
    result.token = env_value;
    if (tok.empty()) {
        return result;
    }
    if (tok == "auto") {
        // "auto" means: use the auto-detected value. We model this as the
        // maximum SimdCap; force_simd_cap() then clamps it to auto_cap.
        // Use a per-arch sentinel so cross-arch comparisons stay sound.
#if IKAFSSN_ARCH_X86
        result.cap = SimdCap::AVX512VBMI2;
#elif IKAFSSN_ARCH_ARM
        result.cap = SimdCap::SVE2;
#else
        result.cap = SimdCap::Scalar;
#endif
        result.explicit_value = true;
        return result;
    }

    struct Entry { const char* name; SimdCap cap; };
    static constexpr Entry table[] = {
        {"scalar",       SimdCap::Scalar},
        {"sse42",        SimdCap::SSE42},
        {"sse4",         SimdCap::SSE42},
        {"avx2",         SimdCap::AVX2},
        {"avx512bw",     SimdCap::AVX512BW},
        {"avx512vbmi",   SimdCap::AVX512VBMI},
        {"avx512vbmi2",  SimdCap::AVX512VBMI2},
        {"neon",         SimdCap::NEON},
        {"asimd",        SimdCap::NEON},
        {"sve",          SimdCap::SVE},
        {"sve2",         SimdCap::SVE2},
        {"sme",          SimdCap::SME},
        {"sme2",         SimdCap::SME2},
    };
    for (const Entry& e : table) {
        if (tok == e.name) {
            result.cap = e.cap;
            result.explicit_value = true;
            return result;
        }
    }
    // Unrecognized → keep explicit_value = false.
    return result;
}

void check_required_tier_or_skip() noexcept {
    const char* req = std::getenv("IKAFSSN_TEST_REQUIRE_TIER");
    if (!req || *req == '\0') return;
    ParseForceResult parsed = parse_force_simd_env(req);
    if (!parsed.explicit_value) {
        // Unrecognized token: do not silently skip — let the caller observe a
        // genuine bug in the parser rather than a phantom skip.
        return;
    }
    SimdCap auto_cap = auto_detected_simd_cap();
    if (static_cast<int>(parsed.cap) > static_cast<int>(auto_cap)) {
        std::fprintf(stderr,
            "SKIP: CPU does not support %s (auto-detected=%s)\n",
            req, simd_cap_name(auto_cap).data());
        std::exit(77);  // CTest standard skip exit code
    }
}

// Always compiled (linkage available to test binaries that define
// IKAFSSN_TESTING). The header gates only the declaration so production code
// cannot accidentally call it.
void reset_simd_dispatch_for_testing() noexcept {
    g_once.~once_flag();
    new (&g_once) std::once_flag();
    g_auto_cap.store(SimdCap::Scalar, std::memory_order_release);
    g_current.store(SimdCap::Scalar, std::memory_order_release);
    g_slow_bmi2.store(false, std::memory_order_release);
    g_avx512f_only.store(false, std::memory_order_release);
    g_initialized.store(false, std::memory_order_release);
}

} // namespace ikafssn
