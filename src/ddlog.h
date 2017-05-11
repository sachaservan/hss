#pragma once

#include <stdint.h>
#include <gmp.h>

#define strip_size 16
#define halfstrip_size ((strip_size)/2)

extern uint32_t lookup[256];
extern uint32_t offset[256];


uint32_t convert(uint64_t *nn);
uint32_t naif_convert(mpz_t n);

void dlog_precompute();
