#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
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



uint32_t convert(uint64_t * nn)
{
  uint32_t steps = 0;

  /**
   * Here we start making a bunch of assumptions.
   * First, we assume the "w" here is 64 bits, which should be the size
   * (in bits) of a mp_limb_t.
   * Secondly, the amount of zeros to check, "d" here is 8.
   */

#define strip_size 8
#define distinguished(x) ((x[23] & (~(ULLONG_MAX >> strip_size)))) == 0
#define lbindex(x) __builtin_clzll(x);

  while (!distinguished(nn)) {
    /**
     * Here we try to find a strip of zeros for "w2" bits.
     * When we find one (up to w2 = 64), then we jump of w = w/2.
     * I tried to optimize this code:
     * - by integrating the if statement above with the for loop invariant;
     * - by making the loop algebraic (i.e. no if-s), given that in the
     *   generated assembly I read a lot of jumps.
     * Unfortunately, both approaches actually lead to a slow down in the code.
     */
    uint64_t x = nn[23];
    if (x == 0) return steps;
    uint32_t first_bit = lbindex(x);
    uint32_t second_bit = 0;

    while (x != 0) {
      /* clear that bit */
      x &= ~(0x8000000000000000 >> first_bit);

      if (x == 0) {
        const uint32_t w = 64 - first_bit;
        if (w > strip_size) {
          return steps + first_bit + 1;
        } else {
          /**
           * We found no distinguished point.
           */
          const uint64_t a = mpn_lshift(nn, nn, 24, 36) * gg;
          mpn_add_1(nn, nn, 24, a);
          steps += 36;
        }
      } else {
        second_bit = lbindex(x);
        if (second_bit - first_bit > strip_size) {
          return steps + first_bit + 1;
        } else {
          first_bit = second_bit;
        }
      }
    }
  }
  return steps;
}


uint32_t naif_convert(mpz_t n)
{
  uint32_t i;
  mpz_t t;
  mpz_init_set_ui(t, 1);
  mpz_mul_2exp(t, t, 1536-strip_size);


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
  getrandom(&_rseed, sizeof(unsigned long int), GRND_NONBLOCK); //GRND_RANDOM
  gmp_randseed_ui(_rstate, _rseed);

  mpz_t n, n0;
  mpz_inits(n, n0, NULL);

  INIT_TIMEIT();
  uint32_t converted;
  for (int i=0; i < 1e5; i++) {
    mpz_urandomm(n0, _rstate, p);
    mpz_set(n, n0);
    START_TIMEIT();
    converted = convert(n->_mp_d);
    END_TIMEIT();
    mpz_set(n, n0);
    assert(converted == naif_convert(n));
  }
  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());


  /* memset(n->_mp_d, 0, 24*8); */
  /* memset(n0->_mp_d, 0, 24*8); */
  /* n0->_mp_d[0] = 13423523; */
  /* n0->_mp_d[1] = 1; */
  /* uint64_t v[64] = {0}; */
  /* unpack(v, n->_mp_d); */
  /* pack(n0, v); */

  mpz_clears(n, n0, p, g, NULL);
  return 0;

}
