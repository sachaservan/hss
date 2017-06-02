#include "config.h"

#include <assert.h>
#include <stdio.h>

#include <gmp.h>

#include "ddlog.h"
#include "hss.h"

int main()
{
  hss_init();


  mpz_t x, y, z;
  uint32_t xc, yc;

  mpz_inits(x, y, z, NULL);

  /* z is a uniformly random group element */
  mpz_urandomm(x, _rstate, q);
  mpz_set_ui(z, 2);
  mpz_powm(z, z, x, p);
  /* y = 2z  (mod p) */
  mpz_mul_2exp(y, z, 1);
  mpz_mod(y, y, p);
  /* x = 2y  (mod p) */
  mpz_mul_2exp(x, y, 1);
  mpz_mod(x, z, p);

  printf("%lx %lx\n", x->_mp_d[23], y->_mp_d[23]);
  xc = convert(x->_mp_d);
  yc = convert(y->_mp_d);

  printf("%d %d    %d\n", xc, yc, distinguished(z));


  hss_clear();
  mpz_clears(x, y, z, NULL);
  return 0;

}
