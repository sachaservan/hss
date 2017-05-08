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
  mpz_abs(rop, rop);
}


void ssl1_share(ssl1_t r1, ssl1_t r2, const mpz_t v, const elgamal_key_t key)
{
  mpz_t zero;
  mpz_init_set_ui(zero, 0);

  r1->w = elgamal_encrypt(key, v);
  for (size_t t = 0; t < 160; t++) {
    if (mpz_tstbit(key.sk, t)) {
      r1->cw[t] = elgamal_encrypt(key, v);
    } else {
      r1->cw[t] = elgamal_encrypt(key, zero);
    }
  }

  mpz_set(r2->w.c1, r1->w.c1);
  mpz_set(r2->w.c2, r1->w.c2);
  //r2->w = elgamal_encrypt(key, v);
  for (size_t t = 0; t < 160; t++) {
    if (mpz_tstbit(key.sk, t)) {
      r2->cw[t] = elgamal_encrypt(key, v);
    } else {
      r2->cw[t] = elgamal_encrypt(key, zero);
    }
  }
}

void ssl1_merge(mpz_t rop, const ssl1_t r1, const ssl1_t r2, const elgamal_key_t key)
{
  mpz_t rop1, rop2;
  mpz_inits(rop1, rop2, NULL);

  elgamal_decrypt(rop1, key, r1->w);
  elgamal_decrypt(rop2, key, r2->w);

  assert(!mpz_cmp(rop1, rop2));
  mpz_set(rop, rop1);

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
    mpz_mul_ui(n, n, 2);
    mpz_mod(n, n, p);
  }

  mpz_clear(t);
  return i;
}

elgamal_key_t key;

void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  mpz_t tmp;
  mpz_init(tmp);

  mpz_powm(tmp, sl1->w.c1, sl2->cx, p);
  mpz_invert(tmp, tmp, p);

  mpz_powm(rop->x, sl1->w.c2, sl2->x, p);
  mpz_mul(rop->x, rop->x, tmp);
  mpz_mod(rop->x, rop->x, p);

  const uint32_t converted = naif_convert(rop->x);
  mpz_set_ui(rop->x, converted);

  mpz_clear(tmp);
}

int main()
{
  /* set up entropy, prime modulus etc. */
  mpz_entropy_init();
  hss_init();

  mpz_t test;
  mpz_init(test);

  mpz_t x;
  mpz_t y;
  mpz_inits(x, y, NULL);

  uint64_t expected;

  /* test elgamal */
  key = elgamal_keygen();
  mpz_urandomb(x, _rstate, 128);
  //elgamal_cipher_t c = elgamal_encrypt(key, x);
  //elgamal_decrypt(test, key, c);
  //assert(!mpz_cmp(x, test));


  ssl1_t r1, r2;
  ssl2_t s1, s2;
  ssl2_t t1, t2;

  ssl1_init(r1);
  ssl1_init(r2);
  ssl2_init(s1);
  ssl2_init(s2);
  ssl2_init(t1);
  ssl2_init(t2);


  INIT_TIMEIT();
  for (int i = 0; i <  (int) 1e1; i++) {

    mpz_rrandomb(x, _rstate, 1);
    mpz_urandomb(y, _rstate, 1);
    //mpz_set_ui(x, 0);
    //mpz_set_ui(y, 1);

    ssl2_share(s1, s2, x, key.sk);
    ssl2_merge(test, s1, s2);
    assert(!mpz_cmp(test, x));

    ssl1_share(r1, r2, y, key);
    ssl1_merge(test, r1, r2, key);
    expected = mpz_cmp_ui(x, 0) && mpz_cmp_ui(y, 0) ? 2 : 1;
    assert(!mpz_cmp_ui(test, expected));

    START_TIMEIT();
    hss_mul(t1, r1, s1);
    END_TIMEIT();
    hss_mul(t2, r2, s2);
    //gmp_printf("resulting shares: %Zx %Zx\n", t1->x, t2->x);
    //gmp_printf("x, y: %Zx %Zx\n", x, y);

    ssl2_merge(test, t2, t1);
    //gmp_printf("result: %Zx\n", test);
    expected = (!mpz_cmp_ui(x, 0) || !mpz_cmp_ui(y, 0)) ? 0 : 1;
    assert(!mpz_cmp_ui(test, expected));

  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
  ssl2_clear(s1);
  ssl2_clear(s2);
  ssl1_clear(r1);
  ssl1_clear(r2);
  ssl2_clear(t1);
  ssl2_clear(t2);

  mpz_clears(x, y, NULL);
  mpz_clears(key.sk, key.pk, NULL);
  hss_del();
  return 0;
}
