#include "config.h"

#include <assert.h>
#include <stdio.h>

#include <gmp.h>

#include "entropy.h"
#include "fbase.h"
#include "group.h"
#include "timeit.h"

int main()
{
  mpz_entropy_init();
  group_init();
  mpz_t test, expected, base;
  uint64_t exp;
  mpz_inits(test, expected, base, NULL);

  mpz_t expected_mod;
  mpz_inits(expected_mod, NULL);


  mpz_urandomm(base, _rstate, p);
  mpz_powm_ui(expected_mod, base, 2, p);
  mpz_pow_ui(test, base, 2);
  remp(test);
  assert(!mpz_cmp(test, expected_mod));
  mpz_clear(expected_mod);

  mpz_urandomm(base, _rstate, p);
  fbase_ptr pbase = fb_init();
  fb_set(pbase, base);

  for (int i = 0; i < (int) 1e4; i++) {
    INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);

    getrandom(&exp, 8, GRND_NONBLOCK);

    mpz_powm_ui(expected, base, exp, p);
    START_TIMEIT();
    fb_powmp_ui(test, pbase, exp);
    END_TIMEIT();

    assert(!mpz_cmp(test, expected));
    printf("%d" SEP TIMEIT_FORMAT "\n", FB_BASE, GET_TIMEIT());
  }


  fb_clear(pbase);
  mpz_clears(test, expected, base, NULL);
}
