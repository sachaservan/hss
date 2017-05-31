#include "config.h"

#include <assert.h>
#include <stdio.h>

#include <gmp.h>

#include "ddlog.h"
#include "hss.h"

int main()
{
  hss_init();


  mpz_t x, y;
  uint32_t xc, yc;

  mpz_inits(x, y, NULL);

  mpz_urandomm(x, _rstate, q);
  mpz_set_ui(y, 2);

  /* y is a uniformly random group element */
  mpz_powm(y, y, x, p);
  /* x = 2y  (mod p) */
  mpz_mul_2exp(x, y, 1);
  mpz_mod(x, x, p);

  xc = convert(x->_mp_d);
  yc = convert(y->_mp_d);
  printf("%d %d\n", xc, yc);

  hss_clear();
  return 0;

}
