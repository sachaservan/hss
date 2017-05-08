#pragma once

#include <gmp.h>
#include <stdint.h>

#include "elgamal.h"

/**
 * p is our prime modulus, and is 2^n - g
 * where g is referred to as "gamma" (built-in function in C, so transliterated)
 */
extern const char* p_str;
extern mpz_t p, q;
extern uint32_t gg;

void hss_init();
void hss_del();


/** A level-1 share is the El-Gamal encryption of a secretly-shared value,
 *  plus the encryption of the product of each bit.
 */
typedef struct ssl1 {
  elgamal_cipher_t w;
  elgamal_cipher_t cw[160];
} ssl1_t[1];


void ssl1_init(ssl1_t s);
void ssl1_clear(ssl1_t s);

/** A level-2 share are subractive shares.
    This shares have at most 192 bits.
 */
typedef struct ssl2 {
  mpz_t x;
  mpz_t cx;
} ssl2_t[1];

#define ssl2_add(rop, a, b) mpz_add(rop, a, b)
void ssl2_init(ssl2_t s);
void ssl2_clear(ssl2_t s);
