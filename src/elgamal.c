#include "config.h"

#include <gmp.h>

#include "elgamal.h"
#include "entropy.h"
#include "group.h"
#include "hss.h"


void elgamal_keygen(elgamal_key_t rop)
{
  mpz_set_ui(rop->pk, 2);

  //  mpz_urandomm(rop->sk, _rstate, q);
  mpz_urandomb(rop->sk, _rstate, 160);
  mpz_powm(rop->pk, rop->pk, rop->sk, p);
}

void elgamal_encrypt(elgamal_cipher_t rop, const elgamal_key_t k, const mpz_t m)
{
  mpz_t x;

  mpz_init(x);
  mpz_urandomm(x, _rstate, q);

  mpz_set_ui(rop->c1, 2);
  mpz_powm(rop->c1, rop->c1, x, p);
  mpz_invert(rop->c1, rop->c1, p);

  /* cached intermediate values */
  mpz_t e;
  mpz_init_set_ui(e, 1);
  mpz_mul_2exp(e, e, 64);
  mpz_powm(rop->c1e64, rop->c1, e, p);
  mpz_powm(rop->c1e128, rop->c1e64, e, p);
  mpz_clear(e);

  mpz_set(rop->c2, k->pk);
  mpz_powm(rop->c2, rop->c2, x, p);
  mpz_set_ui(x, 2);
  mpz_powm(x, x, m, p);
  mpz_mul(rop->c2, rop->c2, x);
  mpz_mod(rop->c2, rop->c2, p);

  mpz_clear(x);

  /* fixed bases */
  //fb_set(rop->fb_c1, rop->c1);
  /* fb_set(rop->fb_c1e64, rop->c1e64); */
  /* fb_set(rop->fb_c1e128, rop->c1e128); */
  /* fb_set(rop->fb_c2, rop->c2); */


}


void elgamal_decrypt(mpz_t rop, const elgamal_key_t k, const elgamal_cipher_t c)
{
  mpz_powm(rop, c->c1, k->sk, p);
  mpz_mul(rop, rop, c->c2);
  mpz_mod(rop, rop, p);
}


void elgamal_cipher_init(elgamal_cipher_t c)
{
  mpz_inits(c->c1, c->c2, c->c1e64, c->c1e128, NULL);   \
  c->fb_c1 = fb_init();
  c->fb_c1e64 = fb_init();
  c->fb_c1e128 = fb_init();
  c->fb_c1e128 = fb_init();
  c->fb_c2 = fb_init();
}

void elgamal_cipher_clear(elgamal_cipher_t c)
{
  mpz_clears(c->c1, c->c2, c->c1e64, c->c1e128, NULL);   \
  fb_clear(c->fb_c1);
  fb_clear(c->fb_c1e64);
  fb_clear(c->fb_c1e128);
  fb_clear(c->fb_c1e128);
  fb_clear(c->fb_c2);
}


void elgamal_cipher_set(elgamal_cipher_t rop, const elgamal_cipher_t op1)
{
  mpz_set(rop->c1, op1->c1);
  mpz_set(rop->c2, op1->c2);
  mpz_set(rop->c1e64, op1->c1e64);
  mpz_set(rop->c1e128, op1->c1e128);

  /* fixed bases */
  //fb_set(rop->fb_c1, rop->c1);
  //fb_set(rop->fb_c1e64, rop->c1e64);
  //fb_set(rop->fb_c1e128, rop->c1e128);
  //fb_set(rop->fb_c2, rop->c2);

}
