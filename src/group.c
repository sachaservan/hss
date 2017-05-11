#include "group.h"

const char* p_str =
  "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "505CAF";

mpz_t p, q;
const uint64_t gg = 11510609;

void group_init()
{
  mpz_init_set_str(p, p_str, 0);

  mpz_init_set(q, p);
  mpz_sub_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);

}


void group_clear()
{
  mpz_clears(p, q, NULL);
}

void powmp_ui(mpz_t x, const mpz_t base, uint64_t exp)
{
  mpz_t y;
  mpz_init_set_ui(y, 1);
  mpz_set(x, base);

  while (exp > 1) {
    if (exp &  1) mpz_mul_modp(y, x, y);
    mpz_mul_modp(x, x, x);
    exp >>= 1;
  }

  mpz_mul_modp(x, x, y);
  mpz_clear(y);
}
