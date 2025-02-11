#pragma once

#include "config.h"
#include <stdbool.h>
#include <stdint.h>

#include <gmp.h>

#define strip_size (FAILURE - 1)
#define halfstrip_size ((strip_size)/2)

extern uint64_t lookup[0x01 << halfstrip_size];
extern uint64_t offset[0x01 << halfstrip_size];

/** default function used for conversion */
#define convert convert_lookup

uint32_t convert_lookup(uint64_t *nn);
uint32_t convert_ec17(uint64_t *nn);
uint32_t convert_naif(mpz_t n);
bool distinguished(mpz_t n);

void dlog_precompute();
