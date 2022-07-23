#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>

#define GF2_BITS 163 /** bits length of a field element */
#define GF2_WORDS 6
#define BIG_INT_BITS 192
#define BIG_INT_WORDS (BIG_INT_BITS / 32)

typedef uint8_t u8;
typedef uint32_t u32;

typedef struct
{
    u32 words[BIG_INT_WORDS];
} BigInt;

typedef struct
{
    u32 words[GF2_WORDS * 2]; /** for enough memory */
} GaloisField163;

void bigint_set(BigInt*,const u32*,u32);
void bigint_output(const BigInt*);
void bigint_rand(BigInt*,u32);
u8 bigint_getbit(const BigInt*,u32);
void bigint_setzero(BigInt*);

void set(GaloisField163*,const u32*, u32);
void output(const GaloisField163*);
bool cmp(const GaloisField163*,const GaloisField163*);
void gf_rand(GaloisField163*,u32);
u8 getbit(const GaloisField163*,u32);
void setzero(GaloisField163*);
void cpy(GaloisField163*,const GaloisField163*);
void add(GaloisField163*,GaloisField163*,GaloisField163*);
void shiftleft(GaloisField163*);
void reduce(GaloisField163*,GaloisField163*);
void multsa(GaloisField163*,GaloisField163*,GaloisField163*);
void multcomb(GaloisField163*,GaloisField163*,GaloisField163*);
void multwindow(GaloisField163*,GaloisField163,GaloisField163*);
void sqr(GaloisField163*,GaloisField163*);
void gf_exp(GaloisField163*,GaloisField163*,const BigInt*);
void inv(GaloisField163*,GaloisField163*);
