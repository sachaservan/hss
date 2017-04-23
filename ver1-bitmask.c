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

typedef __uint128_t uint128_t;

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

static uint8_t lookup[256];
static uint8_t offset[256];

static inline ssize_t
getrandom(void *buffer, size_t length, unsigned int flags)
{
  return syscall(SYS_getrandom, buffer, length, flags);
}


INIT_TIMEIT();
#define strip_size 16
#define halfstrip_size ((strip_size)/2)

static inline
void add_1(uint64_t *b, const size_t start, const size_t len, uint64_t a)
{
  for (size_t i = start; a != 0; i = (i+1)%len) {
    const uint64_t r = b[i] + a;
    a = (r < a);
    b[i] = r;
  }
  /*
     we don't check for overflows in "i".
     If it happens, buy a lottery ticket and retry.
   */
}

uint32_t convert(uint64_t * nn)
{
  static const uint64_t topmask = ~(ULLONG_MAX >> halfstrip_size);
  static const uint64_t topbigmask = ~(ULLONG_MAX >> strip_size);
  static const uint64_t bottommask = (0x01  << halfstrip_size) -1;

  uint32_t w;
  uint32_t steps;
  size_t head = 23;
#define next_head  ((head + 23) % 24)
#define tail       ((head + 1)  % 24)
#define next_tail  ((head + 2)  % 24)

#define distinguished(x) (((x)[head] & topbigmask)) == 0

  /** Check the most significant block */
  const uint64_t x = nn[head];
  for (uint32_t w2 = halfstrip_size; w2 < 64-halfstrip_size; w2 += halfstrip_size) {
    if (!(x & (topmask >> w2))) {
      for (w = w2-1; !(x & (topmask >> w)); w--);
      ++w;
      if (!(x & (topbigmask >> w))) return w;
    }
  }

  for (steps = 64; !distinguished(nn); steps += 64) {
    const uint64_t x = nn[head];
    const uint64_t y = nn[next_head];

    if (!(x & bottommask)) {
      const size_t previous = (x >> halfstrip_size) & bottommask;
      const uint8_t next = y >> (64 - halfstrip_size);
      if (next <= lookup[previous]) return steps - halfstrip_size - offset[previous];
    }

    if (!(y & topmask)) {
      const size_t previous = x & bottommask;
      const uint8_t next = (y >> (64 - 2*halfstrip_size)) & bottommask;
      if (next <= lookup[previous]) return steps - offset[previous];
    }

    for (uint32_t w2 = halfstrip_size; w2 < 64-halfstrip_size; w2 += halfstrip_size) {
      if (!(y & (topmask >> w2))) {
        for (w = w2-1; !(y & (topmask >> w)); w--);
        ++w;
        if (!(y & (topbigmask >> w)))  return steps + w;
      }
    }

    /**
     * We found no distinguished point.
     */
    const uint128_t a = (uint128_t) x * gg;
    const uint64_t al = (uint64_t) a;
    const uint64_t ah = (a >> 64);
    head = next_head;
    nn[tail] = al;
    add_1(nn, next_tail, 24, ah);

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


  for (size_t i = 0; i <= 0xFF; i++) {
    uint32_t j = ffs(i) ? ffs(i) - 1 : 8;
    lookup[i] = 0xFF >> (8-j);
    offset[i] = j;
  }

  uint32_t converted;
  for (int i=0; i < (int) 15258; i++) {
    mpz_t n, n0;
    mpz_inits(n, n0, NULL);

    mpz_urandomm(n0, _rstate, p);
    mpz_set(n, n0);
    START_TIMEIT();
    converted = convert(n->_mp_d);
    END_TIMEIT();
    mpz_set(n, n0);
#ifndef NDEBUG
    uint32_t expected = naif_convert(n);
    //printf("%d %d\n", converted, expected);
    assert(converted == expected);
#endif
    mpz_clears(n, n0, NULL);
  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());

  mpz_clears(p, g, NULL);
  return 0;

}
