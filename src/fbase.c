#include <stdint.h>
#include <stdlib.h>

#include <gmp.h>

#include "fbase.h"
#include "group.h"


fbase_ptr fb_init()
{
  fbase_ptr pb = (fbase_ptr) calloc(FB_FRAMES, sizeof(fbase_unit));

  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= 0xFF; i++) {
      mpz_init(pb[j][i]);
    }
  }
  return pb;
}

void fb_set_small(fbase_t pb, const mpz_t n)
{
  mpz_t e;
  mpz_init(e);
  for (size_t j = 0; j < FB_FRAMES/2; j++) {
    for (size_t i = 0; i <= 0xFF; i++) {
      uint64_t e =  (0x01 <<  8*j) * i;
      powmp_ui(pb[j][i], n, e);
    }
  }
  mpz_clear(e);
}

void fb_copy(fbase_t dst, fbase_t source)
{
  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= 0xFF; i++) {
      mpz_set(dst[j][i], source[j][i]);
    }
  }
}

void fb_set(fbase_t pb, const mpz_t n)
{
  mpz_t e;
  mpz_init(e);
  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= 0xFF; i++) {
      mpz_set_ui(e, 1);
      mpz_mul_2exp(e, e, 8*j);
      mpz_mul_ui(e, e, i);

      mpz_powm(pb[j][i], n, e, p);
    }
  }
  mpz_clear(e);
}

void fb_clear(fbase_t pb)
{

  for (size_t j = 0; j < FB_FRAMES; j++) {
    for (size_t i = 0; i <= 0xFF; i++) {
      mpz_clear(pb[j][i]);
    }
  }
  free(pb);
}
