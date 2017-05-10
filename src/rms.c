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
void modp(mpz_t rop)
{
  const int limbs = rop->_mp_size - 24;

  if (limbs < 0) return;
  else if (limbs == 0 && mpz_cmp(rop, p) < 0) return;
  else if (limbs == 0 && mpz_cmp(rop, p) >= 0) {
    mpz_sub(rop, rop, p);
    return;
  }

  mp_limb_t a[limbs + 10];
  uint64_t gamma = gg;
  mpn_mul(a, rop->_mp_d + 24, limbs, &gamma, 1);
  for (int i = 24; i < rop->_mp_size; i++) rop->_mp_d[i] = 0;
  mpn_add(rop->_mp_d, rop->_mp_d, 24, a, limbs+1);
  mpz_mul_2exp(rop, rop, 1);
  mpz_div_2exp(rop, rop, 1);

  modp(rop);
}

static inline
void fbpowm(mpz_t rop, const mpz_t T[4][256], const uint32_t exp)
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
                      const mpz_t c1e64,
                      const mpz_t c1e128,
                      const mpz_t c2,
                      const uint32_t x,
                      const mpz_t cx)
{

  /* first block */
  mpz_powm_ui(op1, c1, cx->_mp_d[0], p);
  /* second block */
  mpz_powm_ui(op2, c1e64, cx->_mp_d[1], p);
  mpz_mul(op1, op2, op1);
  modp(op1);
  /* third block */
  mpz_powm_ui(op2, c1e128, cx->_mp_d[2], p);
  mpz_mul(op1, op2, op1);
  modp(op1);

  //mpz_powm(op1, c1, cx, p);
  mpz_powm_ui(op2, c2, x, p);
  mpz_mul(op2, op2, op1);
  modp(op2);
  const uint32_t converted = convert(op2->_mp_d);
  return converted;
}

void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  mpz_t op1, op2;
  mpz_inits(op1, op2, NULL);

  rop->x = __mul_single(op1, op2,
                        sl1->w->c1,
                        sl1->w->c1e64,
                        sl1->w->c1e128,
                        sl1->w->c2,
                        sl2->x,
                        sl2->cx);

  mpz_set_ui(rop->cx, 0);
  for (size_t t = 0; t < 160; t++) {
    const uint32_t converted =
      __mul_single(op1, op2,
                   sl1->cw[t]->c1,
                   sl1->cw[t]->c1e64,
                   sl1->cw[t]->c1e128,
                   sl1->cw[t]->c2,
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
  ELGAMAL_KEY(init, key);
  elgamal_keygen(key);


  mpz_t expected_mod, base;
  mpz_inits(expected_mod, base, NULL);
  mpz_urandomm(base, _rstate, p);
  mpz_powm_ui(expected_mod, base, 2, p);
  mpz_pow_ui(test, base, 2);
  modp(test);
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
  ELGAMAL_KEY(clear, key);
  hss_del();
  return 0;
}
