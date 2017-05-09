#include "config.h"

#include <gmp.h>

#include "elgamal.h"
#include "entropy.h"
#include "hss.h"


void elgamal_key_init(elgamal_key_t rop)
{
  mpz_init(rop->sk);
  mpz_init_set_ui(rop->pk, 2);

  //  mpz_urandomm(rop->sk, _rstate, q);
  mpz_urandomb(rop->sk, _rstate, 160);
  mpz_powm(rop->pk, rop->pk, rop->sk, p);
}

void elgamal_key_clear(elgamal_key_t key)
{
  mpz_clears(key->sk, key->pk, NULL);
}


void elgamal_encrypt(elgamal_cipher_t rop, const elgamal_key_t k, const mpz_t m)
{
  mpz_t x;

  mpz_init(x);
  mpz_urandomm(x, _rstate, q);

  mpz_set_ui(rop->c1, 2);
  mpz_powm(rop->c1, rop->c1, x, p);
  mpz_invert(rop->c1, rop->c1, p);

  mpz_set(rop->c2, k->pk);
  mpz_powm(rop->c2, rop->c2, x, p);
  mpz_set_ui(x, 2);
  mpz_powm(x, x, m, p);
  mpz_mul(rop->c2, rop->c2, x);
  mpz_mod(rop->c2, rop->c2, p);

  mpz_clear(x);
}


void elgamal_decrypt(mpz_t rop, const elgamal_key_t k, const elgamal_cipher_t c)
{
  mpz_powm(rop, c->c1, k->sk, p);
  mpz_mul(rop, rop, c->c2);
  mpz_mod(rop, rop, p);
}

void elgamal_cipher_init(elgamal_cipher_t rop)
{
  mpz_inits(rop->c1, rop->c2, NULL);
}
void elgamal_cipher_set(elgamal_cipher_t rop, const elgamal_cipher_t op1)
{
  mpz_set(rop->c1, op1->c1);
  mpz_set(rop->c2, op1->c2);
}

void elgamal_cipher_clear(elgamal_cipher_t op)
{
  mpz_clears(op->c1, op->c2, NULL);
}
