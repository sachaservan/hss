#pragma once

#include "config.h"
#include <stdbool.h>
#include <stdint.h>

#include <gmp.h>

#define strip_size (FAILURE - 1)
#define halfstrip_size ((strip_size)/2)

extern uint64_t lookup[0x01 << strip_size];
extern uint64_t offset[0x01 << strip_size];


uint32_t convert(uint64_t *nn);
uint32_t naif_convert(mpz_t n);
bool distinguished(mpz_t n);

void dlog_precompute();
