#include <stdint.h>
#include <limits.h>
#include <strings.h>

#include <gmp.h>

#include "ddlog.h"
#include "group.h"
#include "hss.h"

typedef __uint128_t uint128_t;

uint8_t lookup[256];
uint8_t offset[256];

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

uint32_t __attribute__((optimize("unroll-loops"))) convert(uint64_t * nn)
{
  static const uint64_t topmask = ~(ULLONG_MAX >> halfstrip_size);
  static const uint64_t topbigmask = ~(ULLONG_MAX >> strip_size);
  static const uint64_t bottommask = (0x01  << halfstrip_size) -1;
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
      const size_t previous = (x >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
      const uint8_t next =    (x >> (64 - halfstrip_size - w2 - halfstrip_size)) & bottommask;
      if (next <= lookup[previous]) return w2 - offset[previous];
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
        const size_t previous = (y >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
        const uint8_t next =    (y >> (64 - halfstrip_size - w2 - halfstrip_size)) & bottommask;
        if (next <= lookup[previous]) return steps + w2 - offset[previous];
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
    mpz_mul_ui(n, n, 2);
    mpz_mod(n, n, p);
  }

  mpz_clear(t);
  return i;
}


void dlog_precompute()
{
  for (size_t i = 0; i <= 0xFF; i++) {
    uint32_t j = ffs(i) ? ffs(i) - 1 : 8;
    lookup[i] = 0xFF >> (8-j);
    offset[i] = j;
  }
}
