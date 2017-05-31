#include "config.h"

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include <gmp.h>

#include "ddlog.h"
#include "elgamal.h"
#include "entropy.h"
#include "group.h"
#include "hss.h"
#include "timeit.h"

static inline
uint32_t mul_single(const elgamal_cipher_t c,
                    const uint32_t x,
                    const mpz_t cx)
{
  mpz_t op1, op2;
  mpz_inits(op1, op2, NULL);

  /* c1: first block */
  fb_powmp_ui(op1, c->fb_c1, cx->_mp_d[0]);
  /* c1: second block */
  fb_powmp_ui(op2, c->fb_c1e64, cx->_mp_d[1]);
  mpz_mul_modp(op1, op2, op1);
  /* c1: third block */
  fb_powmp_ui(op2, c->fb_c1e128, cx->_mp_d[2]);
  mpz_mul_modp(op1, op2, op1);
  /* c2 */
  fb_powmp_ui(op2, c->fb_c2, x);
  mpz_mul_modp(op2, op2, op1);

  const uint32_t converted = convert(PTR(op2));

  mpz_clears(op1, op2, NULL);
  return converted;
}

void hss_mul(ssl2_t rop, const ssl1_t sl1, const ssl2_t sl2)
{
  uint32_t converted;
  rop->x = mul_single(sl1->w, sl2->x, sl2->cx);

  mpz_set_ui(rop->cx, 0);
  for (size_t t = 0; t < SK_BLOCKS; t++) {
    mpz_mul_2exp(rop->cx, rop->cx, SS_BASE);
    converted = mul_single(sl1->cw[t], sl2->x, sl2->cx);
    mpz_add_ui(rop->cx, rop->cx, converted);
  }
}
