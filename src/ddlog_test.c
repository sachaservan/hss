#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <strings.h>

#include "ddlog.h"
#include "entropy.h"
#include "group.h"
#include "timeit.h"

#ifdef NDEBUG
#error "running test without debug?"
#endif

int main()
{

  group_init();
  dlog_precompute();
  mpz_entropy_init();


  for (int i=0; i < (int) (1e3); i++) {
    mpz_t n, n0;
    mpz_inits(n, n0, NULL);

    mpz_urandomm(n0, _rstate, p);
    mpz_set(n, n0);
    uint32_t expected = naif_convert(n);
    mpz_set(n, n0);
    uint32_t converted = convert(n->_mp_d);
    printf("%d %d\n", converted, expected);
    assert(converted == expected);
    mpz_clears(n, n0, NULL);
  }

  group_clear();
  return 0;

}
