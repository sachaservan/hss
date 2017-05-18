#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <gmp.h>

#include "fbase.h"
#include "group.h"


fbase_ptr fb_init()
{
  fbase_ptr pb = (fbase_ptr) calloc(FB_FRAMES, sizeof(fbase_unit));

  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= FB_MASK; i++) {
      mpz_init(pb[j][i]);
    }
  }
  return pb;
}

void fb_set_small(fbase_t pb, const mpz_t n)
{
  mpz_t e;
  mpz_init(e);
  for (size_t j = 0; j < (FB_FRAMES)/2; j++) {
    for (size_t i = 0; i <= FB_MASK; i++) {
      uint64_t e =  (0x01 <<  (FB_BASE)*j) * i;
      powmp_ui(pb[j][i], n, e);

      /* force size to be constant. */
      _mpz_realloc(pb[j][i], 24);
      SIZ(pb[j][i]) = 24;

    }
  }
  mpz_clear(e);
}

void fb_copy(fbase_t dst, fbase_t source)
{
  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= FB_MASK; i++) {
      mpz_set(dst[j][i], source[j][i]);

      /* force size to be constant */
      _mpz_realloc(dst[j][i], 24);
      SIZ(dst[j][i]) = 24;
    }
  }
}

void fb_set(fbase_t pb, const mpz_t n)
{
  mpz_t e;
  mpz_init(e);
  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= FB_MASK; i++) {
      mpz_set_ui(e, 1);
      mpz_mul_2exp(e, e, (FB_BASE)*j);
      mpz_mul_ui(e, e, i);

      mpz_powm(pb[j][i], n, e, p);
    }
  }
  mpz_clear(e);
}

void fb_clear(fbase_t pb)
{

  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= FB_MASK; i++) {
      mpz_clear(pb[j][i]);
    }
  }
  free(pb);
}


void __attribute__((optimize("unroll-loops")))
fb_powmp_ui(mpz_t rop, fbase_t pb, const uint64_t exp)
{
#define e(i) ((exp >> (i * (FB_BASE))) & (FB_MASK))

  mpz_mul_modp(rop, pb[0][e(0)], pb[1][e(1)]);
  for (size_t j = 2; j < FB_FRAMES; j++) {
    if (e(j) != 0) {
      mpz_mul_modp(rop, rop,  pb[j][e(j)]);
    }
  }
}
