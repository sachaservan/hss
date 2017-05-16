#include "config.h"

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include <gmp.h>

#include "ddlog.h"
#include "elgamal.h"
#include "entropy.h"
#include "group.h"
#include "hss.h"
#include "rms.h"
#include "timeit.h"

INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);


static inline
uint32_t mul_single(const elgamal_cipher_t c,
                    const uint32_t x,
                    const mpz_t cx)
{
  mpz_t op1, op2;
  mpz_inits(op1, op2, NULL);

  /* c1: first block */
  fb_powmp_ui(op1, c->fb_c1, cx->_mp_d[0]);
  /* c1: second block */
  fb_powmp_ui(op2, c->fb_c1e64, cx->_mp_d[1]);
  mpz_mul_modp(op1, op2, op1);
  /* c1: third block */
  fb_powmp_ui(op2, c->fb_c1e128, cx->_mp_d[2]);
  mpz_mul_modp(op1, op2, op1);
  /* c2 */
  fb_powmp_ui(op2, c->fb_c2, x);
  mpz_mul_modp(op2, op2, op1);

  START_TIMEIT();
  const uint32_t converted = convert(PTR(op2));
  END_TIMEIT();

  mpz_clears(op1, op2, NULL);
  return converted;
}

void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  uint32_t converted;
  rop->x = mul_single(sl1->w, sl2->x, sl2->cx);

  mpz_set_ui(rop->cx, 0);
  for (size_t t = 0; t < 160; t++) {
    mpz_mul_2exp(rop->cx, rop->cx, 1);
    converted = mul_single(sl1->cw[t], sl2->x, sl2->cx);
    mpz_add_ui(rop->cx, rop->cx, converted);
  }
}

int main()
{
  group_init();
  mpz_entropy_init();
  dlog_precompute();

  mpz_t test;
  mpz_init(test);

  mpz_t x, y, xy;
  mpz_inits(x, y, xy, NULL);

  elgamal_key_t key;
  ELGAMAL_KEY(init, key);
  elgamal_keygen(key);


  mpz_t expected_mod, base;
  mpz_inits(expected_mod, base, NULL);
  mpz_urandomm(base, _rstate, p);
  mpz_powm_ui(expected_mod, base, 2, p);
  mpz_pow_ui(test, base, 2);
  remp(test);
  // gmp_printf("%Zx\n%Zx\n", test, expected_mod);
  assert(!mpz_cmp(test, expected_mod));
  mpz_clear(expected_mod);

  ssl1_t r1, r2;
  ssl2_t s1, s2;
  ssl2_t t1, t2;

  ssl1_init(r1);
  ssl1_init(r2);
  ssl2_init(s1);
  ssl2_init(s2);
  ssl2_init(t1);
  ssl2_init(t2);

  mpz_urandomb(y, _rstate, 1);
  mpz_urandomb(x, _rstate, 1);
  /* mpz_set_ui(x, 1); */
  /* mpz_set_ui(y, 1); */

  ssl2_share(s1, s2, x, key->sk);
  ssl2_open(test, s1, s2);
  assert(!mpz_cmp(test, x));

  ssl1_share(r1, r2, y, key);

  ssl1_open(test, r1, r2, key);
  assert(!mpz_cmp_ui(test, mpz_cmp_ui(y, 0) ? 2 : 1));

  for (int i = 0; i <  (int) 1e2; i++) {
    hss_mul(t1, r1, s1);
    hss_mul(t2, r2, s2);

#ifndef NDEBUG
    gmp_printf("%Zx %Zx\n", x, y);
    gmp_printf("%d %d\n", s1->x, s2->x);
#endif

    mpz_mul(xy, x, y);
    ssl2_open(test, t2, t1);
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
  ELGAMAL_KEY(clear, key);
  group_clear();
  return 0;
}
