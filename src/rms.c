#include "config.h"

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include <gmp.h>

#include "elgamal.h"
#include "entropy.h"
#include "rms.h"
#include "hss.h"
#include "timeit.h"

/** this function is only for testing purposes. */
void ssl2_share(ssl2_t s1, ssl2_t s2, const mpz_t v, const mpz_t sk)
{
  mpz_rrandomb(s1->x, _rstate, 192);
  mpz_add(s2->x, v, s1->x);

  mpz_rrandomb(s1->cx, _rstate, 192);
  mpz_mul(s2->cx, sk, v);
  mpz_add(s2->cx, s2->cx, s1->cx);
}


void ssl2_merge(mpz_t rop, const ssl2_t s1, const ssl2_t s2)
{

  mpz_sub(rop, s2->x, s1->x);
}


void ssl1_share(ssl1_t r1, ssl1_t r2, const mpz_t v, const elgamal_key_t key)
{
  mpz_t x, zero;
  mpz_init(x);
  mpz_init_set_ui(zero, 0);
  /* first share is random */
  mpz_urandomm(x, _rstate, q);

  r1->w = elgamal_encrypt(key, x);
  for (size_t t = 0; t < 160; t++) {
    if (mpz_tstbit(key.sk, t)) {
      r1->cw[t] = elgamal_encrypt(key, x);
    } else {
      r1->cw[t] = elgamal_encrypt(key, zero);
    }
  }

  mpz_add(x, x, v);
  r2->w = elgamal_encrypt(key, x);
  for (size_t t = 0; t < 160; t++) {
    if (mpz_tstbit(key.sk, t)) {
      r2->cw[t] = elgamal_encrypt(key, x);
    } else {
      r2->cw[t] = elgamal_encrypt(key, zero);
    }
  }

  mpz_clear(x);
}

void ssl1_merge(mpz_t rop, const ssl1_t r1, const ssl1_t r2, const elgamal_key_t key)
{
  mpz_t rop1, rop2;
  mpz_inits(rop1, rop2, NULL);

  elgamal_decrypt(rop1, key, r1->w);
  elgamal_decrypt(rop2, key, r2->w);

  mpz_invert(rop, rop1, p);
  mpz_mul(rop, rop, rop2);
  mpz_mod(rop, rop, p);

  mpz_clears(rop1, rop2, NULL);
}


#define strip_size 16
uint32_t naif_convert(mpz_t n)
{
  uint32_t i;
  mpz_t t;
  mpz_init_set_ui(t, 1);
  mpz_mul_2exp(t, t, 1536-strip_size);


  for (i = 0; mpz_cmp(n, t) > -1; i++) {
    mpz_mul_2exp(n, n, 1);
    mpz_mod(n, n, p);
  }

  mpz_clear(t);
  return i;
}

void hss_mul(mpz_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  mpz_t tmp;
  mpz_init(tmp);

  mpz_sub(tmp, q, sl2->cx);
  mpz_powm(tmp, sl1->w.c1, tmp, p);

  mpz_powm(rop, sl1->w.c2, sl2->x, p);
  mpz_mul(rop, rop, tmp);
  mpz_mod(rop, rop, p);

  const uint32_t converted = naif_convert(rop);
  mpz_set_ui(rop, converted);

  mpz_clear(tmp);
}

int main()
{
  /* set up entropy, prime modulus etc. */
  mpz_entropy_init();
  hss_init();

  mpz_t test;
  mpz_init(test);

  /* test elgamal */
  elgamal_key_t key = elgamal_keygen();
  //elgamal_cipher_t c = elgamal_encrypt(key, sk);
  //elgamal_decrypt(test, key, c);
  //assert(!mpz_cmp(sk, test));


  /* end test */

  mpz_t x;
  mpz_t y;

  mpz_init_set_ui(x, 0);
  //mpz_rrandomb(x, _rstate, 1);
  mpz_init_set_ui(y, 1);
  //mpz_urandomb(y, _rstate, 1);

  /* attempt to create a secret share */
  ssl1_t r1, r2;
  ssl2_t s1, s2;
  mpz_t t1, t2;

  ssl1_init(r1);
  ssl1_init(r2);
  ssl2_init(s1);
  ssl2_init(s2);
  mpz_inits(t1, t2, NULL);


  INIT_TIMEIT();

  for (int i = 0; i < 1e1; i++) {

    ssl1_share(r1, r2, y, key);
    ssl2_share(s1, s2, x, key.sk);

    ssl2_merge(test, s1, s2);
    assert(!mpz_cmp(test, x));
    ssl1_merge(test, r1, r2, key);
    assert(!mpz_cmp_ui(test, 2));

    START_TIMEIT();
    hss_mul(t1, r1, s1);
    END_TIMEIT();
    hss_mul(t2, r2, s2);
    ssl2_merge(test, s1, s2);
    //gmp_printf("%Zd %Zd\n", s1, s2);
    assert(!mpz_cmp_ui(test, 0));
  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
  ssl2_clear(s1);
  ssl2_clear(s2);
  ssl1_clear(r1);
  ssl1_clear(r2);

  mpz_clears(t1, t2, NULL);

  mpz_clears(x, y, NULL);
  mpz_clears(key.sk, key.pk, NULL);
  hss_del();
  return 0;
}
