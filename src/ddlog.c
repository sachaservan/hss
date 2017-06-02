#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include <strings.h>

#include <gmp.h>

#include "ddlog.h"
#include "group.h"
#include "hss.h"

uint32_t lookup[256];
uint32_t offset[256];

static const uint64_t topmask = ~(ULLONG_MAX >> halfstrip_size);
//static const uint64_t topbigmask = ~(ULLONG_MAX >> strip_size);
static const uint64_t bottommask = (0x01  << halfstrip_size) -1;

static const uint64_t failuremask = ~(ULLONG_MAX >> FAILURE);
#define distinguished_limb 0x8000000000000000
#define clz(x) __builtin_clzll(x)

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


#define zdistinguished(x) (((x)[head] & topbigmask)) == 0

uint32_t __attribute__((optimize("unroll-loops"))) convert(uint64_t * nn)
{
  uint32_t steps;
  size_t head = 23;
#define next_head  ((head + 23) % 24)
#define tail       ((head + 1)  % 24)
#define next_tail  ((head + 2)  % 24)
#define x (nn[head])
#define y (nn[next_head])

  /** Check the most significant block */
  for (uint32_t w2 = clz(x) + halfstrip_size - (clz(x) % halfstrip_size);
       w2 < 64 - halfstrip_size;
       w2 += halfstrip_size) {
    if (!(x & (topmask >> w2))) {
      const size_t previous = (x >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
      const uint32_t next   = (x >> (64 - halfstrip_size - w2 - halfstrip_size)) & bottommask;
      if (next <= lookup[previous]) return w2 - offset[previous] - 1;
    }
  }

  for (steps = 64; true; steps += 64) {
    if (!(x & bottommask)) {
      const size_t previous = (x >> halfstrip_size) & bottommask;
      const uint32_t next   =  y >> (64 - halfstrip_size);
      if (next <= lookup[previous]) return steps - halfstrip_size - offset[previous] - 1;
    }

    if (!(y & topmask)) {
      const size_t previous =  x & bottommask;
      const uint32_t next   = (y >> (64 - 2*halfstrip_size)) & bottommask;
      if (next <= lookup[previous]) return steps - offset[previous] - 1;
    }

    for (uint32_t w2 = halfstrip_size; w2 < 64-halfstrip_size; w2 += halfstrip_size) {
      if (!(y & (topmask >> w2))) {
        const size_t previous = (y >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
        const uint32_t next   = (y >> (64 - halfstrip_size - w2 - halfstrip_size)) & bottommask;
        if (next <= lookup[previous]) return steps + w2 - offset[previous] - 1;
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
}

bool distinguished(mpz_t n)
{
  return (n->_mp_d[23] & failuremask) == distinguished_limb;
}

uint32_t naif_convert(mpz_t n)
{
  uint32_t steps;

  for (steps = 0; !distinguished(n); steps++) {
    mpz_mul_ui_modp(n, n, 2);
  }

  return steps;
}


void dlog_precompute()
{
  for (size_t i = 0; i <= bottommask; i++) {
    uint32_t j = ffs(i) ? ffs(i) - 1 : halfstrip_size;
    lookup[i] = bottommask >> (halfstrip_size - j);
    offset[i] = j;
  }
}
