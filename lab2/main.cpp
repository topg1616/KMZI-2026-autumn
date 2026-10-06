#include "number_theory.h"
#include "rsa_core.h"
#include "oaep.h"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <x86intrin.h>

using namespace rsa;

namespace {

std::mt19937_64& rng() {
    static std::mt19937_64 gen{std::random_device{}()};
    return gen;
}

Message random_message(std::size_t len) {
    Message m{};
    m.len = len;
    for (std::size_t i = 0; i < len; ++i) m.data[i] = static_cast<std::uint8_t>(rng());
    return m;
}

void random_seed(std::uint8_t seed[H_LEN]) {
    for (std::size_t i = 0; i < H_LEN; ++i) seed[i] = static_cast<std::uint8_t>(rng());
}

BigInt bytes_to_int(const std::array<std::uint8_t, K_LEN>& b) {
    BigInt x = 0;
    for (std::uint8_t c : b) x = (x << 8) | c;
    return x;
}

std::array<std::uint8_t, K_LEN> int_to_bytes(const BigInt& x) {
    std::array<std::uint8_t, K_LEN> b{};
    BigInt t = x;
    for (int i = static_cast<int>(K_LEN) - 1; i >= 0; --i) {
        b[i] = static_cast<std::uint8_t>((t & 0xFF).convert_to<unsigned>());
        t >>= 8;
    }
    return b;
}

// алгоритм Гарнера на малых числах
void demo_garner_small() {
    const BigInt p = 11, q = 13, n = p * q;  // N = 143
    const BigInt e = 7, d = 43;              // lambda = lcm(10, 12) = 60
    const BigInt m = 9;
    const BigInt c = mod_pow(m, e, n);       // 48
    const BigInt dp = d % (p - 1), dq = d % (q - 1), qinv = mod_inverse(q % p, p);
    const BigInt m1 = mod_pow(c % p, dp, p);
    const BigInt m2 = mod_pow(c % q, dq, q);
    BigInt diff = (m1 - m2) % p;
    if (diff < 0) diff += p;
    const BigInt h = (qinv * diff) % p;
    const BigInt restored = m2 + h * q;
    std::cout << "Garner demo: c=" << c << " m1=" << m1 << " m2=" << m2
              << " h=" << h << " restored=" << restored << "\n";
    if (restored != m) throw std::runtime_error("Garner demo failed");
}

}

int main() {
    try {
        demo_garner_small();

        std::cout << "Keygen 2048 bit...\n";
        const auto t0 = std::chrono::steady_clock::now();
        const RsaKeys key = generate_keys(2048);
        const auto t1 = std::chrono::steady_clock::now();
        std::cout << "Keygen: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                  << " ms\n";

        // граничные длины сообщения 0 и MAX_MSG
        for (const std::size_t len : {std::size_t(0), MAX_MSG}) {
            const Message msg = random_message(len);
            std::uint8_t seed[H_LEN];
            random_seed(seed);
            const EncodedBlock em = oaep_encode(msg, seed);
            const BigInt c = rsa_encrypt(bytes_to_int(em.data), key);
            const EncodedBlock em2{int_to_bytes(rsa_decrypt_crt(c, key))};
            const DecodeResult res = oaep_decode(em2);
            bool same = res.ok && res.msg.len == msg.len;
            for (std::size_t i = 0; same && i < len; ++i)
                same = (res.msg.data[i] == msg.data[i]);
            if (!same) throw std::runtime_error("Round-trip failed for len=" + std::to_string(len));
            std::cout << "Round-trip OK, len=" << len << "\n";
        }

        // такты процессора, валидный vs намеренно повреждённый шифротекст
        const Message msg = random_message(32);
        std::uint8_t seed[H_LEN];
        random_seed(seed);
        const EncodedBlock em = oaep_encode(msg, seed);
        std::array<std::uint8_t, K_LEN> good = int_to_bytes(rsa_encrypt(bytes_to_int(em.data), key));
        std::array<std::uint8_t, K_LEN> tampered = good;
        tampered[100] ^= 0x01;

        auto measure = [&](const std::array<std::uint8_t, K_LEN>& ct) {
            const BigInt c = bytes_to_int(ct);
            for (int i = 0; i < 3; ++i)
                oaep_decode(int_to_bytes(rsa_decrypt_crt(c, key)));
            const std::uint64_t t = __rdtsc();
            const int reps = 50;
            for (int i = 0; i < reps; ++i)
                oaep_decode(int_to_bytes(rsa_decrypt_crt(c, key)));
            return (__rdtsc() - t) / reps;
        };
        const std::uint64_t ticks_good = measure(good);
        const std::uint64_t ticks_bad  = measure(tampered);
        std::cout << "Ticks valid=" << ticks_good << " tampered=" << ticks_bad
                  << " diff=" << (100.0 * static_cast<double>(ticks_bad) /
                                  static_cast<double>(ticks_good) - 100.0) << "%\n";

        // Повреждённый шифротекст обязан отклоняться
        {
            const DecodeResult res =
                oaep_decode(int_to_bytes(rsa_decrypt_crt(bytes_to_int(tampered), key)));
            if (res.ok) throw std::runtime_error("Tampered ciphertext accepted!");
            std::cout << "Tampered ciphertext rejected\n";
        }

        // эмпирическое ускорение CRT (~4x)
        {
            const BigInt c = bytes_to_int(good);
            std::uint64_t t = __rdtsc();
            for (int i = 0; i < 10; ++i) rsa_decrypt_crt(c, key);
            const std::uint64_t tc = (__rdtsc() - t) / 10;
            t = __rdtsc();
            for (int i = 0; i < 10; ++i) rsa_decrypt_plain(c, key);
            const std::uint64_t tp = (__rdtsc() - t) / 10;
            std::cout << "Decrypt ticks: CRT=" << tc << " plain=" << tp
                      << " speedup=" << static_cast<double>(tp) / static_cast<double>(tc) << "x\n";
        }

        std::cout << "Self-tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}