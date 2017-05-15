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

  mpz_init(q);
  mpz_sub_ui(q, p, 1);
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

void mul_modp(mpz_t rop, mpz_t op1, mpz_t op2)
{
  mp_limb_t z2[24+2];
  mp_limb_t z0[24+2];
  mp_limb_t z1[24+12]= {0};
  mp_limb_t z1f[12+1];
  mp_limb_t z1s[12+1];
  mp_limb_t z3[24];
  const mp_limb_t *op1_hi = PTR(op1) + 12;
  const mp_limb_t *op2_hi = PTR(op2) + 12;
  const mp_limb_t *op1_lo = PTR(op1);
  const mp_limb_t *op2_lo = PTR(op2);
  const mp_limb_t *z1_hi = z1 + 12;
  const mp_limb_t *z1_lo = z1;
  mp_limb_t carry;

  /* make all operands of the same width */
  _mpz_realloc(rop, 24);
  _mpz_realloc(op1, 24);
  _mpz_realloc(op2, 24);
  SIZ(rop) = SIZ(op1) = SIZ(op2) = 24;

  mpn_mul_n(z2, op1_hi, op2_hi, 12);
  z2[24] = z2[25] = 0;
  mpn_mul_n(z0, op1_lo, op2_lo, 12);
  z0[24] = z0[25] = 0;

  carry = mpn_add_n(z1f, op1_lo, op1_hi, 12);
  z1f[12] = carry;
  carry = mpn_add_n(z1s, op2_lo, op2_hi, 12);
  z1s[12] = carry;
  mpn_mul_n(z1, z1f, z1s, 12+1);
  mpn_sub_n(z1, z1, z2, 24+2);
  mpn_sub_n(z1, z1, z0, 24+2);

  mpn_add_n(z3, z2, z1_hi, 24);
  carry = mpn_mul_1(PTR(rop), z3, 24, gg);
  carry += mpn_add_n(PTR(rop), PTR(rop), z0, 24);
  carry += mpn_add_n(PTR(rop)+12, PTR(rop)+12, z1_lo, 12);
  mpn_add_1(PTR(rop), PTR(rop), 24, carry*gg);

  MPN_NORMALIZE(PTR(rop), SIZ(rop));
}
