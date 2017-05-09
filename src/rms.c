#include "config.h"

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include <gmp.h>

#include "ddlog.h"
#include "elgamal.h"
#include "entropy.h"
#include "rms.h"
#include "hss.h"
#include "timeit.h"


static inline
void fbpowm(mpz_t rop, const_fbptable_t T, const uint32_t exp)
{
  const uint8_t *e = (uint8_t *) &exp;

  mpz_mul(rop, T[0][e[0]], T[1][e[1]]);
  mpz_mod(rop, rop, p);
  mpz_mul(rop, rop, T[2][e[2]]);
  mpz_mod(rop, rop, p);
  mpz_mul(rop, rop, T[3][e[3]]);
  mpz_mod(rop, rop, p);
}


static inline
uint32_t __mul_single(mpz_t op1,
                      mpz_t op2,
                      const mpz_t c1,
                      const mpz_t c2,
                      const_fbptable_t T,
                      const uint32_t x,
                      const mpz_t cx)
{

  mpz_powm(op1, c1, cx, p);
  mpz_invert(op1, op1, p);

  //mpz_t test; mpz_init(test);
  //mpz_powm_ui(test, c2, x, p);
  fbpowm(op2, T, x);
  //if (mpz_cmp(test, op2)) gmp_printf("base: %Zx\nexp: %x\npcomp: %Zx\nreal: %Zd\n", c2, x, op2, test);
  //mpz_clear(test);

  mpz_powm_ui(op2, c2, x, p);
  mpz_mul(op2, op2, op1);
  mpz_mod(op2, op2, p);
  const uint32_t converted = convert(op2->_mp_d);
  return converted;
}

void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  mpz_t op1, op2;
  uint32_t converted;
  mpz_inits(op1, op2, NULL);

  converted = __mul_single(op1, op2,
                           sl1->w->c1,
                           sl1->w->c2,
                           sl1->T,
                           sl2->x,
                           sl2->cx);
  rop->x = converted;

  mpz_set_ui(rop->cx, 0);
  for (size_t t = 0; t < 160; t++) {
    converted = __mul_single(op1, op2,
                             sl1->cw[t]->c1,
                             sl1->cw[t]->c2,
                             sl1->T,
                             sl2->x,
                             sl2->cx);
    mpz_add_ui(rop->cx, rop->cx, converted);
    mpz_mul_2exp(rop->cx, rop->cx, 1);
  }
  mpz_div_2exp(rop->cx, rop->cx, 1);

  mpz_clears(op1, op2, NULL);
}

int main()
{
  mpz_entropy_init();
  hss_init();
  dlog_precompute();

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
    ssl2_open(test, s1, s2);
    assert(!mpz_cmp(test, x));

    ssl1_share(r1, r2, y, key);
    ssl1_open(test, r1, r2, key);
    assert(!mpz_cmp_ui(test, mpz_cmp_ui(y, 0) ? 2 : 1));

    START_TIMEIT();
    hss_mul(t1, r1, s1);
    END_TIMEIT();
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
  elgamal_key_clear(key);
  hss_del();
  return 0;
}
