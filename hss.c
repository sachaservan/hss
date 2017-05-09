#define _GNU_SOURCE
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdio.h>
#include <sys/time.h>

#include <gmp.h>

const static char* p_str =
  "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "505CAF";
mpz_t p;

mpz_t modxp;

/* this goes in a separate file, hold on */
gmp_randstate_t _rstate;
unsigned long int _rseed;

#include <linux/random.h>
#include <sys/syscall.h>


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


static inline ssize_t
getrandom(void *buffer, size_t length, unsigned int flags)
{
  return syscall(SYS_getrandom, buffer, length, flags);
}

void mpz_entropy_init()
{
  gmp_randinit_default(_rstate);
  getrandom(&_rseed, sizeof(unsigned long int), GRND_NONBLOCK);
  gmp_randseed_ui(_rstate, _rseed);
}

typedef struct pbase  {
  mpz_t base;
  mpz_t T[4][256];
} pbase_t;

void naif_fbpowm(mpz_t rop, const pbase_t pb, uint32_t exp)
{
  mpz_set(rop, pb.base);
  mpz_powm_ui(rop, rop, exp, p);
}


void fbprecompute(pbase_t *pb)
{
  for (size_t j = 0; j < 4; j++) {
    for (size_t i = 0; i <= 0xFF; i++) {
      uint64_t e = (0x01 << 8*j) * i;
      mpz_t *t = pb->T[j];
      mpz_init(t[i]);
      mpz_powm_ui(t[i], pb->base, e, p);
    }
  }
}

static inline
void fbpowm(mpz_t rop, const pbase_t * const pb, uint32_t exp)
{
  const uint8_t *e = (uint8_t *) &exp;

  mpz_mul(rop, pb->T[0][e[0]], pb->T[1][e[1]]);
  mpz_mod(rop, rop, p);
  mpz_mul(rop, rop, pb->T[2][e[2]]);
  mpz_mod(rop, rop, p);
  mpz_mul(rop, rop, pb->T[3][e[3]]);
  mpz_mod(rop, rop, p);

}

int main()
{

  mpz_t x;
  uint32_t exp;
  mpz_init_set_str(p, p_str, 0);

  mpz_entropy_init();

  mpz_init(x);

  /* the size of the cyclic subgrup <2>
   * is #\GG / gcd(#\GG, 2) */
  mpz_init_set(modxp, p);
  mpz_sub_ui(modxp, modxp, 1);
  mpz_div_2exp(modxp, modxp, 1);

  pbase_t pb;
  mpz_init_set_ui(pb.base, 2);
  //  mpz_urandomm(pb.base, _rstate, p);
  fbprecompute(&pb);

  mpz_t expected;
  mpz_init(expected);

  INIT_TIMEIT();
  for (int i = 0; i < (int) 1e4; i++) {
    getrandom(&exp, sizeof(exp), GRND_NONBLOCK);
    START_TIMEIT();
    fbpowm(x, &pb, exp);
    END_TIMEIT();
    //naif_fbpowm(expected, pb, exp);
    //assert(!mpz_cmp(expected, x));
  }

  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());
  mpz_clear(x);
  mpz_clear(pb.base);
  return 0;
}
