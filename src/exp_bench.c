#include <assert.h>
#include <stdio.h>

#include <gmp.h>

#include "entropy.h"
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
  // gmp_printf("%Zx\n%Zx\n", test, expected_mod);
  assert(!mpz_cmp(test, expected_mod));
  mpz_clear(expected_mod);

  INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);
  for (int i = 0; i < (int) 1e4; i++) {
    getrandom(&exp, 4, GRND_NONBLOCK);
    mpz_urandomm(base, _rstate, p);

    mpz_powm_ui(expected, base, exp, p);
    mpz_set(test, base);
    START_TIMEIT();
    powmp_ui(test, test, exp);
    END_TIMEIT();

    assert(!mpz_cmp(test, expected));
  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
  mpz_clears(test, expected, base, NULL);
}
