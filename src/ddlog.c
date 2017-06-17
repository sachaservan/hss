#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include <strings.h>

#include <gmp.h>

#include "ddlog.h"
#include "group.h"
#include "hss.h"

uint64_t lookup[0x01 << halfstrip_size];
uint64_t offset[0x01 << halfstrip_size];

static const uint64_t topmask = ~(ULLONG_MAX >> halfstrip_size);
static const uint64_t bottommask = (0x01  << halfstrip_size) - 1;
static const uint64_t topbigmask = ~(ULLONG_MAX >> strip_size);

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



#define offset_xy_l (64 % halfstrip_size)
#define offset_xy_r (halfstrip_size - offset_xy_l)
#define mask_xy_l   ((0x01 << offset_xy_l) - 1)
#define mask_xy_r   (~(ULLONG_MAX >> offset_xy_r))
#define strip_xy    (((x & mask_xy_l) << offset_xy_r) | (y >> (64 - offset_xy_r)))
#define offset_yz_l ((64 - offset_xy_r) % halfstrip_size)
#define offset_yz_r (halfstrip_size - offset_yz_l)
#define mask_yz_l   ((0x01 << offset_yz_l) - 1)
#define strip_yz    (((y & mask_yz_l) << offset_yz_r) | (z >> (64 - offset_yz_r)))

#define next_head      ((head + 23) % 24)
#define next_next_head ((head + 23*2) % 24)
#define tail           ((head + 1)  % 24)
#define next_tail      ((head + 2)  % 24)


uint32_t __attribute__((optimize("unroll-loops"))) convert_lookup(uint64_t * nn)
{
  uint32_t steps;
  size_t head = 23;
  uint32_t w2;

  /** Check the most significant block */
  const uint64_t x = nn[head];
  for (w2 = clz(x) + halfstrip_size - (clz(x) % halfstrip_size);
       w2 < 64 - offset_xy_l - halfstrip_size;
       w2 += halfstrip_size) {
    if (!(x & (topmask >> w2))) {
      const size_t previous = (x >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
      const uint32_t next   = (x >> (64 - halfstrip_size - w2 - halfstrip_size)) & bottommask;
      if (next <= lookup[previous]) {
        return w2 - offset[previous] - 1;
      }
    }
  }

  for (steps = 64; true; steps += 64) {
    const uint64_t x = nn[head];
    const uint64_t y = nn[next_head];
    const uint64_t z = nn[next_next_head];

    if (!((x >> offset_xy_l) & bottommask)) {
      const size_t previous = (x >> (halfstrip_size + offset_xy_l)) & bottommask;
      const uint64_t next   =  strip_xy;
      if (next <= lookup[previous]) {
        return steps - offset_xy_l - halfstrip_size - offset[previous] - 1;
      }
    }

    if (!strip_xy) {
      const size_t previous = (x >> offset_xy_l) & bottommask;
      const uint64_t next   = ((y << offset_xy_r) >> (64 - halfstrip_size));
      if (next <= lookup[previous]) {
        return steps - offset_xy_l - offset[previous] - 1;
      }
    }

    if (!((y << offset_xy_r) & topmask)) {
      const size_t previous = strip_xy;
      const uint64_t next   = ((y << (offset_xy_r + halfstrip_size)) >> (64 - halfstrip_size));
      if (next <= lookup[previous]) {
        return steps + offset_xy_r - offset[previous] - 1;
      }
    }

    for (w2 = halfstrip_size + offset_xy_r;
         w2 < 64 - halfstrip_size - offset_yz_l;
         w2 += halfstrip_size) {
      if (!(y & (topmask >> w2))) {
        const size_t previous = (y >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
        const uint64_t next   = (y >> (64 - halfstrip_size - w2 - halfstrip_size)) & bottommask;
        if (next <= lookup[previous]) {
          return steps + w2 - offset[previous] - 1;
        }
      }
    }

    if (!(y & (topmask >> w2))) {
      const size_t previous = (y >> (64 - halfstrip_size - w2 + halfstrip_size)) & bottommask;
      const uint64_t next   = strip_yz;
      if (next <= lookup[previous]) {
        return steps + w2 - offset[previous] - 1;
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


uint32_t __attribute__((optimize("unroll-loops"))) convert_ec17(uint64_t * nn)
{
  uint32_t steps;
  size_t head = 23;
  uint32_t w, w2;

  /** Check the most significant block */
  const uint64_t x = nn[head];
  for (w2 = clz(x) + halfstrip_size - (clz(x) % halfstrip_size);
       w2 < 64 - offset_xy_l - halfstrip_size;
       w2 += halfstrip_size) {
    if (!(x & (topmask >> w2))) {
        for (w = w2-1; !(x & (topmask >> w)); w--);
        ++w;
        if (!(x & (topbigmask >> w)))  return w - 1;
    }
  }

  for (steps = 64; true; steps += 64) {
    const uint64_t x = nn[head];
    const uint64_t y = nn[next_head];
    const uint64_t z = nn[next_next_head];

    if (!((x >> offset_xy_l) & bottommask)) {
      for (w=0; !(x >> (halfstrip_size + offset_xy_l + w) & 1); w++);
      if (!(strip_xy >> w)) {
        return steps - offset_xy_l - halfstrip_size - w - 1;
      }
    }

    if (!strip_xy) {
      for (w = 0; !((x >> (offset_xy_l + w)) & 1); w++);
      if (!(y << offset_xy_r & (topmask << w)))
        return steps - offset_xy_l - w - 1;
    }

    if (!((y << offset_xy_r) & topmask)) {
      for (w = 0; !((strip_xy >> w) & 1); w++);
      if (!((y << (offset_xy_r + halfstrip_size)) & (topmask << w)))
        return steps + offset_xy_r - w - 1;
    }

    for (w2 = halfstrip_size + offset_xy_r;
         w2 < 64 - halfstrip_size - offset_yz_l;
         w2 += halfstrip_size) {
      if (!(y & (topmask >> w2))) {
        for (w = w2-1; !(y & (topmask >> w)); w--);
        ++w;
        if (!(y & (topbigmask >> w)))  return steps + w - 1;
      }
    }

    if (!(y & (topmask >> w2))) {
      for (w = w2-1; !(y & (topmask >> w)); w--);
      w++;
      if (!(strip_yz >> (w2 - w))) {
        return steps + w - 1;
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

void dlog_precompute()
{
  for (size_t i = 0; i <= bottommask; i++) {
    uint32_t j = ffs(i) ? ffs(i) - 1 : halfstrip_size;
    lookup[i] = bottommask >> (halfstrip_size - j);
    offset[i] = j;
  }
}
