#include "config.h"

#include <assert.h>
#include <stdlib.h>

#include "entropy.h"
#include "elgamal.h"
#include "hss.h"

void ssl1_init(ssl1_t s)
{
  ELGAMAL_CIPHER(init, s->w);

  for (size_t t = 0; t < SK_BLOCKS; t++) {
    ELGAMAL_CIPHER(init, s->cw[t]);
  }
}

void ssl1_clear(ssl1_t s)
{
  ELGAMAL_CIPHER(clear, s->w);

  for (size_t t = 0; t < SK_BLOCKS; t++) {
    ELGAMAL_CIPHER(clear, s->cw[t]);
  }
}



void ssl1_share(ssl1_t r1, ssl1_t r2, const mpz_t v, const elgamal_key_t key)
{
  mpz_t q, r, x;

  mpz_init_set(q, key->sk);
  mpz_inits(r, x, NULL);

  elgamal_encrypt_shares(r1->w, r2->w, key, v);

  for (size_t t = 0; t < SK_BLOCKS; t++) {
    mpz_fdiv_r_2exp(r, q, SS_BASE);
    mpz_fdiv_q_2exp(q, q, SS_BASE);
    mpz_mul(x, v, r);
    /* do it in reverse so that when computing it's just incremental */
    elgamal_encrypt_shares(r1->cw[SK_BLOCKS - 1 - t],
                           r2->cw[SK_BLOCKS - 1 - t],
                           key, x);
  }

  mpz_clears(q, r, x, NULL);
}

void ssl1_open(mpz_t rop, const ssl1_t r1, const ssl1_t r2, const elgamal_key_t key)
{
  mpz_t rop1, rop2;
  mpz_inits(rop1, rop2, NULL);

  elgamal_decrypt(rop1, key, r1->w);
  elgamal_decrypt(rop2, key, r2->w);

  assert(!mpz_cmp(rop1, rop2));
  mpz_set(rop, rop1);

  mpz_clears(rop1, rop2, NULL);
}


void ssl2_init(ssl2_t s)
{
  s->x = 0;
  mpz_inits(s->cx, NULL);
}

void ssl2_clear(ssl2_t s)
{
  mpz_clear(s->cx);
}


void ssl2_share(ssl2_t s1, ssl2_t s2, const mpz_t v, const mpz_t sk)
{
  /* sampling one byte here is already sufficient.
   * However, the purpose of this function is testing,
   * so here we go sampling over the whole space */
  getrandom(&s1->x, 3, GRND_NONBLOCK);
  //mpz_urandomb(s1->x, _rstate, 192);
  //mpz_add(s2->x, v, s1->x);

  const uint32_t _v = (uint32_t) mpz_get_ui(v);
  s2->x = s1->x + _v;

  mpz_urandomb(s1->cx, _rstate, 192);
  mpz_mul(s2->cx, sk, v);
  mpz_add(s2->cx, s2->cx, s1->cx);
}


void ssl2_open(mpz_t rop, const ssl2_t s1, const ssl2_t s2)
{
  if (s1->x > s2->x) mpz_set_ui(rop, s1->x - s2->x);
  else               mpz_set_ui(rop, s2->x - s1->x);
  //mpz_sub(rop, s2->x, s1->x);
  //mpz_abs(rop, rop);
}
