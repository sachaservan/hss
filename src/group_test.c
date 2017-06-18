#include <assert.h>

#include <gmp.h>

#include "group.h"
#include "entropy.h"

mpz_t test, expected, x, y;


int main()
{
  group_init();
  mpz_entropy_init();
  mpz_inits(test, expected, x, y, NULL);

  for (int i = 0; i < 1e4; i++) {
    mpz_urandomm(x, _rstate, p);
    mpz_urandomm(y, _rstate, p);

    mpz_mul(expected, x, y);
    mpz_mod(expected, expected, p);
    mul_modp(test, y, x);
    if (mpz_cmp(test, expected)) {

      gmp_printf("%lu %lu\n%lu %lu\n",
                 PTR(test)[1], PTR(test)[0],
                 PTR(expected)[1], PTR(expected)[0]);
      assert(0);
    }
  }


  mpz_clears(test, expected, x, y, NULL);
  group_clear();

  return 0;
}
