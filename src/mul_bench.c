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
  mpz_t test, expected, x, y;
  mpz_inits(test, expected, x, y, NULL);

  printf("e6time\n");
  for (int n = 0; n < 0x01 << 16; n++) {
    INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);

    for (int i = 0; i < (int) 1e6; i++) {
      mpz_urandomm(x, _rstate, p);
      mpz_urandomm(y, _rstate, p);

      START_TIMEIT();
      mul_modp(expected, x, y);
      END_TIMEIT();
    }
    printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
  }

  mpz_clears(test, expected, x, y, NULL);
}
