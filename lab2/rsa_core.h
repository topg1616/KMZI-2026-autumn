#pragma once
#include "number_theory.h"

namespace rsa {

BigInt rsa_encrypt(const BigInt& m, const RsaKeys& k);       // c = m^e mod N
BigInt rsa_decrypt_crt(const BigInt& c, const RsaKeys& k);   // алгоритм Гарнера
BigInt rsa_decrypt_plain(const BigInt& c, const RsaKeys& k); // m = c^d mod N

}