#define _GNU_SOURCE
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <linux/random.h>
#include <sys/syscall.h>
#include <sys/time.h>

#include <gmp.h>


#define START_TIMEIT()                                          \
  struct timeval __start, __end; gettimeofday(&__start, NULL)

#define END_TIMEIT()                                                    \
  gettimeofday(&__end, NULL);                                           \
  double __sdiff = (__end.tv_sec - __start.tv_sec), __udiff = (__end.tv_usec - __start.tv_usec)


#define GET_TIMEIT()                            \
  __sdiff + __udiff * 1e-6

#define TIMEIT_FORMAT "%lf"

const static char *p_str =
  "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
  "505CAF";

static inline ssize_t
getrandom(void *buffer, size_t length, unsigned int flags)
{
  return syscall(SYS_getrandom, buffer, length, flags);
}


mpz_t p, t;

uint32_t convert(mpz_t n)
{
  uint32_t i;

  for (i = 0; mpz_cmp(n, t) > -1; i++) {
    mpz_mul_2exp(n, n, 1);
    mpz_mod(n, n, p);
  }
  return i;
}

int main()
{
  // mpz_t p;
  mpz_init_set_str(p, p_str, 0);

  gmp_randstate_t _rstate;
  unsigned long int _rseed;

  gmp_randinit_default(_rstate);
  getrandom(&_rseed, sizeof(unsigned long int), GRND_RANDOM);
  gmp_randseed_ui(_rstate, _rseed);

  // mpz_t t;
  mpz_init_set_str(t, "1", 0);
  mpz_mul_2exp(t, t, 1536-8);

  mpz_t n;
  mpz_init(n);

  START_TIMEIT();
  for (int i=0; i < 5e3; i++) {
    mpz_urandomm(n, _rstate, p);
    convert(n);
  }
  END_TIMEIT();
  printf(TIMEIT_FORMAT "\n", GET_TIMEIT());

  return 0;

}
