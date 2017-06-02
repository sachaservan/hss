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

typedef __uint128_t uint128_t;

void group_init();
void group_clear();



/* some gmp internal funcitons to speed up modulus… */

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)
#define MPN_COPY(dst, src, l) \
  do {                                                  \
    for (int i = 0; i < l; i++) (dst)[i] = (src)[i]; \
  } while(0)

#define MPN_ZERO(dst, l)                                \
  do {                                                  \
    for (size_t i = 0; i < l; i++) (dst)[i] = 0;        \
  } while(0)

#define MPN_NORMALIZE(DST, NLIMBS)                                      \
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
void mul_modp(mpz_t rop, mpz_t op1, mpz_t op2);

#define mpz_mul_modp(rop, op1, op2)             \
  do {                                          \
    mpz_mul(rop, op1, op2);                     \
    remp(rop);                                  \
  } while (0)

#define mpz_mul_ui_modp(rop, op1, op2)          \
  do {                                          \
    mpz_mul_ui(rop, op1, op2);                  \
    remp(rop);                                  \
  } while (0)
