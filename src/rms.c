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
  mpz_urandomb(s1->x, _rstate, 192);
  mpz_add(s2->x, v, s1->x);

  mpz_urandomb(s1->cx, _rstate, 192);
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

  elgamal_encrypt_shares(r1->w, r2->w, key, v);
  for (size_t t = 0; t < 160; t++) {
    if (mpz_tstbit(key->sk, t)) {
      elgamal_encrypt_shares(r1->cw[t], r2->cw[t], key, v);
    } else {
      elgamal_encrypt_shares(r1->cw[t], r2->cw[t], key, zero);
    }
  }

  mpz_clear(zero);
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

static inline
uint32_t __mul_single(mpz_t op1,
                      mpz_t op2,
                      const mpz_t c1,
                      const mpz_t c2,
                      const mpz_t x,
                      const mpz_t cx)
{

  mpz_powm(op1, c1, cx, p);
  mpz_invert(op1, op1, p);

  mpz_powm(op2, c2, x, p);
  mpz_mul(op2, op2, op1);
  mpz_mod(op2, op2, p);
  return naif_convert(op2);
}

void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  mpz_t op1, op2;
  uint32_t converted;
  mpz_inits(op1, op2, NULL);

  converted = __mul_single(op1, op2,
                           sl1->w->c1, sl1->w->c2, sl2->x, sl2->cx);
  mpz_set_ui(rop->x, converted);

  mpz_set_ui(rop->cx, 0);
  for (ssize_t t = 159; t >= 0; t--) {
    converted = __mul_single(op1, op2,
                             sl1->cw[t]->c1, sl1->cw[t]->c2, sl2->x, sl2->cx);
    mpz_add_ui(rop->cx, rop->cx, converted);
    mpz_mul_2exp(rop->cx, rop->cx, 1);
  }
  mpz_div_2exp(rop->cx, rop->cx, 1);

  mpz_clears(op1, op2, NULL);
}

int main()
{
  /* set up entropy, prime modulus etc. */
  mpz_entropy_init();
  hss_init();

  mpz_t test;
  mpz_init(test);

  mpz_t x, y, xy;
  mpz_inits(x, y, xy, NULL);

  elgamal_key_t key;
  elgamal_key_init(key);

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

    mpz_urandomb(y, _rstate, 1);
    mpz_urandomb(x, _rstate, 1);
    /* mpz_set_ui(x, 1); */
    /* mpz_set_ui(y, 1); */

    ssl2_share(s1, s2, x, key->sk);
    ssl2_merge(test, s1, s2);
    assert(!mpz_cmp(test, x));

    ssl1_share(r1, r2, y, key);
    ssl1_merge(test, r1, r2, key);
    assert(!mpz_cmp_ui(test, mpz_cmp_ui(y, 0) ? 2 : 1));

    START_TIMEIT();
    hss_mul(t1, r1, s1);
    END_TIMEIT();
     hss_mul(t2, r2, s2);
#ifndef NDEBUG
    gmp_printf("%Zx %Zx\n", x, y);
    gmp_printf("%Zx %Zx\n", s1->x, s2->x);
#endif

    mpz_mul(xy, x, y);
    ssl2_merge(test, t2, t1);
    assert(!mpz_cmp(test, xy));

    mpz_sub(test, t2->cx, t1->cx);
    mpz_abs(test, test);
    assert(((!mpz_cmp_ui(xy, 1) && !mpz_cmp(test, key->sk))) ||
           ((!mpz_cmp_ui(xy, 0)) && !mpz_cmp_ui(test, 0)));
  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
  ssl2_clear(s1);
  ssl2_clear(s2);
  ssl1_clear(r1);
  ssl1_clear(r2);
  ssl2_clear(t1);
  ssl2_clear(t2);

  mpz_clears(x, y, NULL);
  elgamal_key_clear(key);
  hss_del();
  return 0;
}
