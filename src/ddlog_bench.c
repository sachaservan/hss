#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <strings.h>

#include "ddlog.h"
#include "entropy.h"
#include "group.h"
#include "timeit.h"

int main()
{
  group_init();
  dlog_precompute();
  mpz_entropy_init();

  INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);
  uint32_t converted;
  for (int i=0; i < (int) (1e9 / (0x01 << strip_size) * 3.5); i++) {
    mpz_t n, n0;
    mpz_inits(n, n0, NULL);

    mpz_urandomm(n0, _rstate, p);
    mpz_set(n, n0);
    START_TIMEIT();
    converted = convert(n->_mp_d);
    END_TIMEIT();
    mpz_set(n, n0);
#ifndef NDEBUG
    uint32_t expected = naif_convert(n);
    printf("%d %d\n", converted, expected);
    assert(converted == expected);
#endif
    mpz_clears(n, n0, NULL);
  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());

  group_clear();
  return 0;

}
