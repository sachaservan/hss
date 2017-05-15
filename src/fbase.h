#pragma once

#include <gmp.h>
#include <stdint.h>

#include "group.h"

#define FB_FRAMES 8
typedef mpz_t fbase_unit[256];
typedef fbase_unit   fbase_t[FB_FRAMES];
typedef fbase_unit  *fbase_ptr;

void fb_powmp_ui(mpz_t rop, fbase_t pb, const uint64_t exp);

fbase_ptr fb_init();
void fb_set(fbase_t pb, const mpz_t n);
void fb_set_small(fbase_t pb, const mpz_t n);
void fb_clear(fbase_t pb);
void fb_copy(fbase_t source, fbase_t dst);

#define fb_init_set(pb, base) fb_init(pb); fb_set(pb, base)
