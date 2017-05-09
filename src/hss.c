#include "config.h"

#include <assert.h>
#include <stdlib.h>

#include "entropy.h"
#include "hss.h"

/**
 * p is our prime modulus, and is 2^n - g
 * where g is referred to as "gamma" (built-in function in C, so transliterated)
 */
const char* p_str =
  "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "505CAF";

mpz_t p, q;
uint32_t gg = 11510609;

void hss_init()
{
  mpz_init_set_str(p, p_str, 0);

  mpz_init_set(q, p);
  mpz_sub_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);
}


void hss_del()
{
  mpz_clear(p);
}

void ssl1_init(ssl1_t s)
{
  elgamal_cipher_init(s->w);

  for (size_t t = 0; t < 160; t++) {
    elgamal_cipher_init(s->cw[t]);
  }
}

void ssl1_clear(ssl1_t s)
{
  elgamal_cipher_clear(s->w);

  for (size_t t = 0; t < 160; t++) {
    elgamal_cipher_clear(s->cw[t]);
  }
}

void ssl1_share(ssl1_t r1, ssl1_t r2, const mpz_t v, const elgamal_key_t key)
{
  mpz_t zero;
  mpz_init_set_ui(zero, 0);

  elgamal_encrypt_shares(r1->w, r2->w, key, v);

  for (size_t t = 0; t < 160; t++) {
    if (mpz_tstbit(key->sk, 159-t)) {
      elgamal_encrypt_shares(r1->cw[t], r2->cw[t], key, v);
    } else {
      elgamal_encrypt_shares(r1->cw[t], r2->cw[t], key, zero);
    }
  }

  mpz_clear(zero);
}

void ssl1_open(mpz_t rop, const ssl1_t r1, const ssl1_t r2, const elgamal_key_t key)
{
  mpz_t rop1, rop2;
  mpz_inits(rop1, rop2, NULL);

  elgamal_decrypt(rop1, key, r1->w);
  elgamal_decrypt(rop2, key, r2->w);

  assert(!mpz_cmp(rop1, rop2));
  mpz_set(rop, rop1);

  mpz_clears(rop1, rop2, NULL);
}


void ssl2_init(ssl2_t s)
{
  s->x = 0;
  mpz_inits(s->cx, NULL);
}

void ssl2_clear(ssl2_t s)
{
  mpz_clear(s->cx);
}


void ssl2_share(ssl2_t s1, ssl2_t s2, const mpz_t v, const mpz_t sk)
{
  /* sampling one byte here is already sufficient.
   * However, the purpose of this function is testing,
   * so here we go sampling over the whole space */
  getrandom(&s1->x, 3, GRND_NONBLOCK);
  //mpz_urandomb(s1->x, _rstate, 192);
  //mpz_add(s2->x, v, s1->x);

  const uint32_t _v = (uint32_t) mpz_get_ui(v);
  s2->x = s1->x + _v;

  mpz_urandomb(s1->cx, _rstate, 192);
  mpz_mul(s2->cx, sk, v);
  mpz_add(s2->cx, s2->cx, s1->cx);
}


void ssl2_open(mpz_t rop, const ssl2_t s1, const ssl2_t s2)
{
  if (s1->x > s2->x) mpz_set_ui(rop, s1->x - s2->x);
  else               mpz_set_ui(rop, s2->x - s1->x);
  //mpz_sub(rop, s2->x, s1->x);
  //mpz_abs(rop, rop);
}
