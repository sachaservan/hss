#include "config.h"

#include <assert.h>
#include <stdio.h>

#include <gmp.h>

#include "ddlog.h"
#include "hss.h"
#include "timeit.h"


int main()
{
  hss_init();

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
    INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);

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

    printf("%d" SEP  "%d" SEP "%d" SEP TIMEIT_FORMAT "\n",
           FAILURE, FB_BASE, SS_BASE, GET_TIMEIT());
  }

  ssl2_clear(s1);
  ssl2_clear(s2);
  ssl1_clear(r1);
  ssl1_clear(r2);
  ssl2_clear(t1);
  ssl2_clear(t2);

  mpz_clears(x, y, NULL);
  ELGAMAL_KEY(clear, key);
  hss_clear();
  return 0;
}
