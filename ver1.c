#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <linux/random.h>
#include <sys/syscall.h>
#include <sys/time.h>

#include <gmp.h>


#define INIT_TIMEIT() \
  struct timeval __start, __end;                \
  double __sdiff = 0, __udiff = 0

#define START_TIMEIT()                          \
  gettimeofday(&__start, NULL)

#define END_TIMEIT()                                                    \
  gettimeofday(&__end, NULL);                                           \
  __sdiff += (__end.tv_sec - __start.tv_sec);                           \
  __udiff += (__end.tv_usec - __start.tv_usec)

#define GET_TIMEIT()                            \
  __sdiff + __udiff * 1e-6

#define TIMEIT_FORMAT "%lf"


/**
 * p is our prime modulus, and is 2^n - g
 * where g is referred to as "gamma" (built-in function in C, so transliterated)
 */
const static char* p_str =
  "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "505CAF";
const static char* g_str =
  "11510609";

mpz_t p, g;
uint32_t gg = 11510609;

static inline ssize_t
getrandom(void *buffer, size_t length, unsigned int flags)
{
  return syscall(SYS_getrandom, buffer, length, flags);
}



uint32_t convert(mpz_t n)
{
  uint32_t steps = 0;
  size_t i = 0;

  /**
   * Here we start making a bunch of assumptions.
   * First, we assume the "w" here is 64 bits, which should be the size
   * (in bits) of a mp_limb_t.
   * Secondly, the amount of zeros to check, "d" here is 8.
   */
  assert(sizeof(mp_limb_t) == 8);
  assert(n->_mp_size == 24);
#define mp_size 24

  uint64_t * restrict nn = n->_mp_d;

#define distinguished(x) ((x[23] & (ULLONG_MAX << (64-8))) == 0)

  while (!distinguished(nn)) {
    for (i = 0; i < 64; i+=4) {
      if (((nn[23] | nn[22]) & (0x0F << i)) != 0) break;
    }
    if (i == 64) {
      /**
       * We found no distinguished point for the next 64 steps.
       * Boost it!
       */
      const __int128_t a = nn[23] * gg;
      mpn_lshift(nn, nn, 64, 0);
      mpn_add_n(nn, nn, (mp_limb_t *) &a, 24); // YOLO
      steps += 64;
    } else {
      for (; i < 64; i+=4) {
        if ((nn[23] & (0x0F << i)) != 0)  break;
      }
    }
    if (i == 64) {
      /**
       * We found no distinguished point for the next 32 steps.
       * Boost it!
       */
      const uint64_t a = nn[23] >> 32 * gg;
      mpn_lshift(nn, nn, 32, 0);
      mpn_add_1(nn, nn, 24, a);
      steps += 32;
    } else if (i >= 32) {
      /**
       * We found no distinguished point for the next 16 steps.
       * Boost it!
       */
      const uint64_t a = (nn[23] & (ULLONG_MAX << 32)) * gg;
      mpn_lshift(nn, nn, 16, 0);
      mpn_add_1(nn, nn, 24, a);
      steps += 16;
    } else {
      /**
       * If there is nothing else to do, then just multiply by two.
       */
      if (mpn_lshift(nn, nn, 24, 1)) {
        mpn_add_1(nn, nn, 24, gg);
      }
      steps++;
    }
  }
  return steps;
}


uint32_t naif_convert(mpz_t n)
{
  uint32_t i;
  mpz_t t;
  mpz_init_set_ui(t, 1);
  mpz_mul_2exp(t, t, 1536-8);


  for (i = 0; mpz_cmp(n, t) > -1; i++) {
    mpz_mul_2exp(n, n, 1);
    mpz_mod(n, n, p);
  }

  mpz_clear(t);
  return i;
}

int main()
{
  mpz_init_set_str(p, p_str, 0);
  mpz_init_set_str(g, g_str, 10);

  gmp_randstate_t _rstate;
  unsigned long int _rseed;

  gmp_randinit_default(_rstate);
  getrandom(&_rseed, sizeof(unsigned long int), GRND_RANDOM);
  gmp_randseed_ui(_rstate, _rseed);

  mpz_t n, n0;
  mpz_inits(n, n0, NULL);

  INIT_TIMEIT();
  uint32_t converted;
  for (int i=0; i < 1e4; i++) {
    mpz_urandomm(n0, _rstate, p);
    mpz_set(n, n0);
    START_TIMEIT();
    converted = convert(n);
    END_TIMEIT();
    mpz_set(n, n0);
    assert(converted == naif_convert(n));
  }
  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());

  mpz_clears(n, n0, p, g, NULL);
  return 0;

}
