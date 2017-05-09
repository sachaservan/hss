#pragma once

#include <stdint.h>

#define strip_size 16
#define halfstrip_size ((strip_size)/2)

typedef __uint128_t uint128_t;


extern uint8_t lookup[256];
extern uint8_t offset[256];


uint32_t convert(uint64_t *nn);
void dlog_precompute();
