#pragma once
#include <linux/random.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <gmp.h>

extern gmp_randstate_t _rstate;
extern unsigned long int _rseed;

static inline ssize_t
getrandom(void *buffer, size_t length, unsigned int flags)
{
  return syscall(SYS_getrandom, buffer, length, flags);
}


void mpz_entropy_init();
