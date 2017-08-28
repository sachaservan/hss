#include "config.h"

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

  const int repeat = (0x01 << FAILURE) * 3;
  for (int i=0; i < repeat; i++) {
    mpz_t n;
    mpz_init(n);

    mpz_urandomm(n, _rstate, p);
    INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);

    START_TIMEIT();
    uint32_t steps = convert(n->_mp_d);
    END_TIMEIT();

    printf("%d" SEP "%u" SEP TIMEIT_FORMAT "\n",
           FAILURE, steps, GET_TIMEIT());

    mpz_clear(n);
  }

  group_clear();
  return 0;

}
