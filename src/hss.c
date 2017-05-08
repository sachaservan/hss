#include "config.h"

#include "hss.h"

/**
 * p is our prime modulus, and is 2^n - g
 * where g is referred to as "gamma" (built-in function in C, so transliterated)
 */
const char* p_str =
  "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "505CAF";

mpz_t p, q;
uint32_t gg = 11510609;

void hss_init()
{
  mpz_init_set_str(p, p_str, 0);

  mpz_init_set(q, p);
  mpz_sub_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);
}


void hss_del()
{
  mpz_clear(p);
}

void ssl1_init(ssl1_t s)
{
  /* mpz_init(s->w); */
  /* for (size_t t = 0; t < 160; t++) { */
  /*   mpz_init(s->cw[t]); */
  /* } */
}

void ssl1_clear(ssl1_t s)
{
  mpz_clears(s->w.c1, s->w.c2, NULL);

  for (size_t t = 0; t < 160; t++) {
    mpz_clears(s->cw[t].c1, s->cw[t].c2, NULL);
  }
}

void ssl2_init(ssl2_t s)
{
  mpz_inits(s->x, s->cx, NULL);
}

void ssl2_clear(ssl2_t s)
{
  mpz_clear(s->x);
  mpz_clear(s->cx);
}
