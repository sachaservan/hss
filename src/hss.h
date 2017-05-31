#pragma once

#include <gmp.h>
#include <stdint.h>

#include "elgamal.h"

#define SK_SIZE   160
#define SK_BLOCKS ((SK_SIZE)/(SS_BASE))

void hss_init();
void hss_del();


/** A level-1 share is the El-Gamal encryption of a secretly-shared value,
 *  plus the encryption of the product of each bit.
 */

typedef struct ssl1 {
  elgamal_cipher_t w;
  elgamal_cipher_t cw[SK_BLOCKS];
} ssl1_t[1];


void ssl1_init(ssl1_t s);
void ssl1_clear(ssl1_t s);
void ssl1_share(ssl1_t r1, ssl1_t r2, const mpz_t v, const elgamal_key_t key);
void ssl1_open(mpz_t rop, const ssl1_t r1, const ssl1_t r2, const elgamal_key_t key);

/** A level-2 share are subractive shares.
    This shares have at most 192 bits.
 */
typedef struct ssl2 {
  uint32_t x;
  mpz_t cx;
} ssl2_t[1];

#define ssl2_add(rop, a, b) mpz_add(rop, a, b)
void ssl2_init(ssl2_t s);
void ssl2_clear(ssl2_t s);
void ssl2_share(ssl2_t s1, ssl2_t s2, const mpz_t v, const mpz_t sk);
void ssl2_open(mpz_t rop, const ssl2_t s1, const ssl2_t s2);


void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2);
