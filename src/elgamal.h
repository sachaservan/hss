#pragma once
#include "config.h"

#include <gmp.h>

#include "entropy.h"
#include "fbase.h"

typedef struct elgamal_key {
  mpz_t sk;
  mpz_t pk;
} elgamal_key_t[1];


typedef struct elgamal_cipher {
  /* NOTE: the first part is already inverted! */
  mpz_t c1;
  mpz_t c2;

  mpz_t c1e64;
  mpz_t c1e128;

  fbase_ptr fb_c1;
  fbase_ptr fb_c1e64;
  fbase_ptr fb_c1e128;
  fbase_ptr fb_c2;

} elgamal_cipher_t[1];

void elgamal_encrypt(elgamal_cipher_t c, const elgamal_key_t k, const mpz_t m);
void elgamal_decrypt(mpz_t rop, const elgamal_key_t k, const elgamal_cipher_t c);

void elgamal_keygen(elgamal_key_t key);

void elgamal_cipher_init(elgamal_cipher_t c);
void elgamal_cipher_clear(elgamal_cipher_t c);
void elgamal_cipher_set(elgamal_cipher_t rop, const elgamal_cipher_t op1);

#define ELGAMAL_KEY(func, k)                   \
  mpz_ ## func ## s(k->pk, k->sk, NULL);

/* we're not instantiating just mpz_t anymore :( */
#define ELGAMAL_CIPHER(func, c)                \
  elgamal_cipher_ ## func(c)



#define elgamal_encrypt_shares(share1, share2, key, plaintext) \
  elgamal_encrypt(share1, key, plaintext);                     \
  elgamal_cipher_set(share2, share1)
