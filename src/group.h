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



/* some gmp internal funcitons to speed up modulus… */

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)
#define MPN_NORMALIZE(DST, NLIMBS)				\
  do {									\
    while (1)								\
      {									\
	if ((DST)[(NLIMBS) - 1] != 0)					\
	  break;							\
	(NLIMBS)--;							\
      }									\
  } while (0)



static inline
void remp(mpz_t rop)
{
  int32_t limbs = SIZ(rop) - 24;

  while (limbs > 0) {
    /* note: this declarations MUST happen after checking
     * for positivity of limbs. */
    uint64_t a[limbs+1];
    /* copy the most significant part of rop into a,
     * then set it to zero */
    for (int i = 24; i < SIZ(rop); i++) {
      a[i-24] = PTR(rop)[i];
      PTR(rop)[i] = 0;
    }
    a[limbs] = 0;

    mpn_addmul_1(PTR(rop), a, limbs+1, gg);
    MPN_NORMALIZE(PTR(rop), SIZ(rop));
    limbs = SIZ(rop) - 24;
  }

  if (mpn_cmp(PTR(rop), PTR(p), SIZ(rop)) >= 0) {
    mpn_sub(PTR(rop), PTR(rop), 24, PTR(p), 24);
    MPN_NORMALIZE(PTR(rop), SIZ(rop));
  }

}

void powmp_ui(mpz_t rop, const mpz_t base, uint64_t exp);

#define mpz_mul_modp(rop, op1, op2) mpz_mul(rop, op1, op2); remp(rop);
