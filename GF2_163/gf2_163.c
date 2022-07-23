#include "gf2_163.h"

#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>



//static GaloisField163 IRR_POLY = {{0xC9,0x00,0x00,0x00,0x00,0x08}}; /** irreducible polynomial */
/**
static void err(const char* fnc)
{
    fprintf(stderr,"[-] function: %s failed. Error: %s\n",fnc,strerror(errno));
    exit(EXIT_FAILURE);
}
*/
void bigint_set(BigInt* element,const u32* dat,u32 len)
{
    for(u32 i = 0; (i < len)&&(i < BIG_INT_WORDS); i++)
        element->words[i] = dat[i];
    for(u32 i = len; i < BIG_INT_WORDS; i++)
        element->words[i] = 0;
}

void bigint_output(const BigInt* element)
{
    for(u32 i = 0; i < BIG_INT_WORDS; i++)
        printf("%08X ",element->words[i]);
    printf("\n");
}

void bigint_rand(BigInt* res, u32 bits)
{
    u32 i;
    for(i = 0; i*32 < bits; i++)
    {
        static int fd = -1;
        ssize_t ret;
        size_t len = 32;
        u32* out_ptr = &(res->words[i]);
        while(fd == -1)
        {
            fd = open("/dev/urandom",O_RDONLY);
            if(fd == -1 && errno == EINTR)
                continue;
            else if(fd == -1)
                abort();
        }
        while(len > 0)
        {
            ret = read(fd,out_ptr,len);
            if(ret == -1 && errno == EINTR)
                continue;
            else if(ret == -1)
                abort();
            out_ptr += ret;
            len -= ret;
        }
    }
    if(bits % 32 != 0)
        res->words[i-1]&=0xFFFFFFFF>>(32-bits%32);
}

inline u8 bigint_getbit(const BigInt* element,u32 bit)
{
    return (element->words[bit >> 0x5] >> (bit & 0x01F)) & 0x01;
}

void bigint_setzero(BigInt* element)
{
    for(u32 i = 0; i < BIG_INT_WORDS; i++)
        element->words[i] = 0;
}

void set(GaloisField163* element,const u32* dat,u32 length)
{
    for(u32 i = 0; i < length; i++)
        element->words[i] = dat[i];
    for(u32 i = length; i < GF2_WORDS*2; i++)
        element->words[i] = 0x00;
}

void output(const GaloisField163* element)
{
    for(u32 i = 0; i < GF2_WORDS; i++)
        printf("%08X ",element->words[i]);
    printf("\n");
}

bool cmp(const GaloisField163* first, const GaloisField163* second)
{
    for(u32 i = 0; i < GF2_WORDS*2; i++)
        if(first->words[i] != second->words[i])
            return false;
    return true;
}

void gf_rand(GaloisField163* res, u32 bits)
{
    u32 i;
    for(i = 0; i*32 < bits; i++)
    {
        static int fd = -1;
        ssize_t ret;
        size_t len = 32;
        u32* out_ptr = &(res->words[i]);
        while(fd == -1)
        {
            fd = open("/dev/urandom",O_RDONLY);
            if(fd == -1 && errno == EINTR)
                continue;
            else if(fd == -1)
                abort();
        }
        while(len > 0)
        {
            ret = read(fd,out_ptr,len);
            if(ret == -1 && errno == EINTR)
                continue;
            else if(ret == -1)
                abort();
            out_ptr += ret;
            len -= ret;
        }
    }
    if(bits % 32 != 0)
        res->words[i-1]&=0xFFFFFFFF>>(32-bits%32);
}

inline u8 getbit(const GaloisField163* element,u32 bit)
{
    return (element->words[bit>>5]>>(bit&0x01F))&0x01;
}

void setzero(GaloisField163* element)
{
    for(u32 i = 0; i < GF2_WORDS*2; i++)
        element->words[i] = 0;
}

void cpy(GaloisField163* dst,const GaloisField163* src)
{
    for(u32 i = 0; i < GF2_WORDS*2;i++)
        dst->words[i] = src->words[i];
}

void add(GaloisField163* res, GaloisField163* a, GaloisField163* b)
{
    for(u32 i = 0; i < GF2_WORDS*2; i++)
        res->words[i] = a->words[i]^b->words[i];
}

void shiftleft(GaloisField163* element)
{
    for(u32 i = GF2_WORDS*2 - 1; i > 0; i--)
        element->words[i] = (element->words[i] << 1) | (element->words[i-1] >> 31);
    element->words[0] <<= 1;
}

void reduce(GaloisField163* res,GaloisField163* element)
{
    for(u32 i = 10; i > 5; i--)
    {
        u32 t = element->words[i];
        element->words[i-6] = element->words[i-6]^(t << 29);
        element->words[i-5] = element->words[i-5]^(t<<4)^(t<<3)^t^(t>>3);
        element->words[i-4] = element->words[i-4]^(t>>28)^(t>>29);
    }
    u32 t = element->words[5] & 0xFFFFFFF8;
    element->words[0]=element->words[0]^(t<<4)^(t<<3)^t^(t>>3);
    element->words[1]=element->words[1]^(t>>28)^(t>>29);
    element->words[5]&=0x00000007;
    for(u32 i = 6; i < GF2_WORDS*2;i++)
        element->words[i] = 0x00;
    cpy(res,element);
}

void multsa(GaloisField163* res, GaloisField163* a, GaloisField163* b)
{
    if(getbit(b,0)==1)
        cpy(res,a);
    else
        setzero(res);

    for(u32 i = GF2_BITS-1; i > 0; i--)
    {
        shiftleft(a);
        reduce(a,a);
        if(getbit(b,i)==1)
            add(res,res,a);
    }
}

void multcomb(GaloisField163* res, GaloisField163* a, GaloisField163* b)
{
    setzero(res);

    for(u32 k = 0; k < 32; k++)
    {
        for(u32 j = 0; j < (u32)(GF2_BITS/32); j++)
        {
            if(getbit(a,32*j+k)==1)
                for(u32 i = 0; i <= (u32)(GF2_BITS/32); i++)
                    res->words[i+j] ^= b->words[i];
        }
        if(k != 31)
            shiftleft(b);
    }
}

void sqr(GaloisField163* res, GaloisField163* element)
{
    for(u32 i = 0; i < GF2_WORDS * 2; i++)
    {
        if(i%2 == 1)
            res->words[i] = 0x00;
        res->words[i] = element->words[i/2];
    }
    reduce(res,res);
}

void gf_exp(GaloisField163* res, GaloisField163* base,const BigInt* potent)
{
    /**
     * Using multiply-square algorithm to compute res = base^potent
     */
    setzero(res);
    res->words[5] = 0x01;

    for(u32 i = 0; i < BIG_INT_BITS; i++)
    {
        if(bigint_getbit(potent,i)==1)
            multcomb(res,res,base);
        if(bigint_getbit(potent,i)==0)
            sqr(base,base);
    }
}




