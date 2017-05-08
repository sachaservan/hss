#include "config.h"

#include <gmp.h>

#include "elgamal.h"
#include "entropy.h"
#include "hss.h"


elgamal_key_t elgamal_keygen()
{
  elgamal_key_t rop;
  mpz_init(rop.sk);
  mpz_init_set_ui(rop.pk, 2);

  mpz_urandomm(rop.sk, _rstate, q);
  mpz_powm(rop.pk, rop.pk, rop.sk, p);
  return rop;
}

elgamal_cipher_t elgamal_encrypt(elgamal_key_t k, mpz_t m)
{
  mpz_t x;
  elgamal_cipher_t rop;

  mpz_init(x);
  mpz_urandomm(x, _rstate, q);

  mpz_init_set_ui(rop.c1, 2);
  mpz_powm(rop.c1, rop.c1, x, p);

  mpz_init_set(rop.c2, k.pk);
  mpz_powm(rop.c2, rop.c2, x, p);
  mpz_set_ui(x, 2);
  mpz_powm(x, x, m, p);
  mpz_mul(rop.c2, rop.c2, x);
  mpz_mod(rop.c2, rop.c2, p);

  mpz_clear(x);
  return rop;
}


void elgamal_decrypt(mpz_t rop, elgamal_key_t k, elgamal_cipher_t c)
{
  mpz_powm(rop, c.c1, k.sk, p);
  mpz_invert(rop, rop, p);
  mpz_mul(rop, rop, c.c2);
  mpz_mod(rop, rop, p);
}
