#include "config.h"

#include <stdio.h>

#include <gmp.h>
#include <openssl/ossl_typ.h>
#include <openssl/bio.h>
#include <openssl/bn.h>
#include <openssl/evp.h>
#include <openssl/ec.h>

#include "entropy.h"
#include "group.h"
#include "timeit.h"

void EC_POINT_get_random(const EC_GROUP *group, EC_POINT *r, BN_CTX *ctx) {
    BIGNUM *k = NULL;
    k = BN_new();

    if (!EC_GROUP_get_order(group, k, ctx)) goto err;
    if (!BN_pseudo_rand(k, BN_num_bits(k), 0, 0)) goto err;
    if (!EC_POINT_mul(group, r, k, NULL, NULL, ctx)) goto err;
    if (!EC_POINT_is_on_curve(group, r, ctx)) goto err;

 err:
    if (k) BN_free(k);
}

int main()
{
  mpz_entropy_init();
  group_init();
  mpz_t x, y, xy;
  mpz_inits(x, y, xy, NULL);

  BN_CTX *ctx;
  EC_GROUP *group;
  EC_POINT *P, *Q;
  ctx = BN_CTX_new();

  group = EC_GROUP_new_by_curve_name(NID_secp521r1);


  for (int n = 0; n < (int) 1e3; n++) {
    /* block for Z_p with our reminder */
    {
      INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);
      for (int i = 0; i < (int) 1e5; i++) {
        mpz_urandomm(x, _rstate, p);
        mpz_urandomm(y, _rstate, p);

        START_TIMEIT();
        mul_modp(xy, x, y);
        END_TIMEIT();
      }
      printf(TIMEIT_FORMAT SEP, GET_TIMEIT());
    }
    /* block for Z_p with vanilla reminder */
    {
      INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);
      for (int i = 0; i < (int) 1e5; i++) {
        mpz_urandomm(x, _rstate, p);
        mpz_urandomm(y, _rstate, p);

        START_TIMEIT();
        mpz_mul(xy, x, y);
        mpz_mod(xy, xy, p);
        END_TIMEIT();
      }
      printf(TIMEIT_FORMAT SEP, GET_TIMEIT());
    }
    /* EC  addition */
    {
      INIT_TIMEIT(CLOCK_PROCESS_CPUTIME_ID);
      for (int i = 0; i < (int) 1e5; i++) {
        P = EC_POINT_new(group);
        Q = EC_POINT_new(group);
        EC_POINT_get_random(group, Q, ctx);

        START_TIMEIT();
        EC_POINT_add(group, Q, Q, P, ctx);
        END_TIMEIT();
      }
      printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
    }
  }

  mpz_clears(x, y, xy, NULL);
}
