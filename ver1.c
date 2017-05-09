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
#include <immintrin.h>

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
static const uint32_t strip_size = 16;
#define halfstrip_size (strip_size/2)
static uint8_t lookup[256];
static uint8_t offset[256];

uint32_t convert(uint64_t *nn)
{
  assert(strip_size == 16);
  uint32_t steps;
  static const uint32_t window = 7;

#define distinguished(x) ((x)[23] & ~(ULLONG_MAX >> strip_size)) == 0
  const __m128i rotmask = _mm_set_epi8(0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);
  for (steps = 0; !distinguished(nn); steps += window*8) {
    START_TIMEIT();

    __m128i x = _mm_lddqu_si128((__m128i *) (nn + 22));
    __m128i mask = _mm_set_epi8(0xFF, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    for (int32_t i = 14; i > 0; --i) {
      mask = _mm_shuffle_epi8(mask, rotmask);
        const bool zero = _mm_testz_si128(mask, x);
        if (zero) {
          const uint8_t previous = _mm_extract_epi8(x, i-1);
          const uint8_t next = _mm_extract_epi8(x, i+1);
          if (previous <= lookup[next]) {
            END_TIMEIT();
            return steps + (15-i)*8 - offset[next];
          }
        }
    }
    END_TIMEIT();

    /**
     * We found no distinguished point.
     */
    const uint64_t a = mpn_lshift(nn, nn, 24, window) * gg;
    mpn_add_1(nn, nn, 24, a);
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

  for (uint32_t i = 0; i <= 0xFF; i++) {
    uint32_t j = ffs(i) ? ffs(i) - 1 : 8;
    lookup[i] = 0xFF >> (8-j);
    offset[i] = j;
  }

  mpz_t n, n0;
  mpz_inits(n, n0, NULL);

  uint32_t converted;
  for (uint64_t i=0; i < 1e3; i++) {
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

  mpz_clears(n, n0, NULL);
  mpz_clears(p, g, NULL);
  return 0;

}
