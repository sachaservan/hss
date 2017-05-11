#pragma once

#include <stdint.h>
#include <gmp.h>

/**
 * p is our prime modulus, and is 2^n - g
 * where g is referred to as "gamma" (built-in function in C, so transliterated)
 */
extern const char* p_str;
extern mpz_t p, q;
extern const uint64_t gg;

void group_init();
void group_clear();
