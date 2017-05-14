#pragma once

#include <gmp.h>
#include <stdint.h>

#include "group.h"

#define FB_FRAMES 8
typedef mpz_t fbase_unit[256];
typedef fbase_unit   fbase_t[FB_FRAMES];
typedef fbase_unit  *fbase_ptr;

static inline void __attribute__((optimize("unroll-loops")))
fb_powmp_ui(mpz_t rop, fbase_t pb, const uint64_t exp)
{
  const uint8_t *e = (uint8_t *) &exp;

  mpz_mul_modp(rop, pb[0][e[0]], pb[1][e[1]]);
  for (size_t j = 2; j < FB_FRAMES; j++) {
    const size_t exp = e[j];
    if (exp != 0) {
      mpz_mul_modp(rop, rop,  pb[j][exp]);
    }
  }
}

fbase_ptr fb_init();
void fb_set(fbase_t pb, const mpz_t n);
void fb_set_small(fbase_t pb, const mpz_t n);
void fb_clear(fbase_t pb);
void fb_copy(fbase_t source, fbase_t dst);

#define fb_init_set(pb, base) fb_init(pb); fb_set(pb, base)
