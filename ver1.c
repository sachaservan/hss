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


INIT_TIMEIT();
static const uint32_t strip_size = 8;

uint32_t convert(uint64_t * nn)
{
  uint32_t steps = 0;

  /**
   * Here we start making a bunch of assumptions.
   * First, we assume the "w" here is 64 bits, which should be the size
   * (in bits) of a mp_limb_t.
   * Secondly, the amount of zeros to check, "d" here is 8.
   */
  static const uint32_t window = 32;
  static const uint32_t strip_size = 16;
  static const uint32_t halfstrip_size = strip_size/2;
  static const uint64_t topmask = ~(ULLONG_MAX >> halfstrip_size);
  static const uint64_t topbigmask = ~(ULLONG_MAX >> strip_size);
  //static const uint64_t tailmask = (1 << halfstrip_size) - 1;

  /** Check if x has a halfstrip starting from start */
#define halfstrip(x, start) (((x)) << (start)) & topmask)
#define distinguished(x) (((x)[23] & (~(ULLONG_MAX >> strip_size)))) == 0
#define lbindex(x) __builtin_clzll(x);
#define tbindex(x) __builtin_ctzll(x);

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
    uint32_t w;

    START_TIMEIT();
    for (uint32_t w2 = halfstrip_size; w2 < window*2-halfstrip_size; w2 += halfstrip_size) {
      if (!(x & (topmask >> w2))) {
        for (w = w2-1; !(x & (topmask >> w)); w--);
        ++w;
        if (!(x & (topbigmask >> w))) return steps + w;
      }
    }
    END_TIMEIT();
    /**
     * We found no distinguished point.
     */
    const uint64_t a = mpn_lshift(nn, nn, 24, window) * gg;
    mpn_add_1(nn, nn, 24, a);
    steps += window;
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
  getrandom(&_rseed, sizeof(unsigned long int), GRND_NONBLOCK);
  gmp_randseed_ui(_rstate, _rseed);

  mpz_t n, n0;
  mpz_inits(n, n0, NULL);

  uint32_t converted;
  for (int i=0; i < 1e4; i++) {
    mpz_urandomm(n0, _rstate, p);
    mpz_set(n, n0);
    converted = convert(n->_mp_d);
    mpz_set(n, n0);
#ifndef NDEBUG
    uint32_t expected = naif_convert(n);
    if (converted != expected) printf("%d %d\n", converted, expected);
    assert(converted == expected);
#endif
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
