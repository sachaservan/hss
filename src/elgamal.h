#pragma once
#include "config.h"

#include <gmp.h>

#include "entropy.h"

typedef struct elgamal_key {
  mpz_t sk;
  mpz_t pk;
} elgamal_key_t;

typedef struct elgamal_cipher {
  mpz_t c1;
  mpz_t c2;
} elgamal_cipher_t;

elgamal_key_t elgamal_keygen();
elgamal_cipher_t elgamal_encrypt(elgamal_key_t k, const mpz_t m);
void elgamal_decrypt(mpz_t rop, elgamal_key_t k, elgamal_cipher_t c);
