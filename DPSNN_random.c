// DPSNN_random.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "DPSNN_random.h"

/* getRandom_r(uint32_t *seed, uint32_t max1) is based on rand_r(unsigned int * seed)*/

/* Reentrant random function from POSIX.1c.
Copyright (C) 1996, 1999, 2009 Free Software Foundation, Inc.
This file is part of the GNU C Library.
Contributed by Ulrich Drepper <drepper@cygnus.com>, 1996.
The GNU C Library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
The GNU C Library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with the GNU C Library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA. */

/* This algorithm is mentioned in the ISO C standard, here extended
for 32 bits. */

uint32_t getRandom_r (uint32_t *seed, uint32_t max1)
{
uint32_t next = *seed;
uint32_t result;
next *= 1103515245;
next += 12345;
result = (uint32_t) (next / 65536) % 2048;
next *= 1103515245;
next += 12345;
result <<= 10;
result ^= (uint32_t) (next / 65536) % 1024;
next *= 1103515245;
next += 12345;
result <<= 10;
result ^= (uint32_t) (next / 65536) % 1024;
result = result % max1;
*seed = next;
 return result;
}


// WARNING THE KISS IS NEVER USED UP TO NOW IN THE CODE 

/*

 KISS - downloaded from www.helbreth.org 2012

 the idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
        x(n)=a*x(n-1)+1 mod 2^32
        y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
        z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
      2^32*(2^32-1)*(2^63+2^32-1) > 2^127

#define ulong unsigned long

static ulong kiss_x = 1;
static ulong kiss_y = 2;
static ulong kiss_z = 4;
static ulong kiss_w = 8;
static ulong kiss_carry = 0;
static ulong kiss_k;
static ulong kiss_m;
*/

// WARNING THE KISS FUNCTIONS ARE NEVER USED UP TO NOW IN THE DPSNN CODE
// the function defined by DPSNN_random.h is used 

static unsigned long kiss_x, kiss_y, kiss_z, kiss_w, kiss_carry, kiss_k, kiss_m;

void seed_rand_kiss(unsigned long seed) {
    kiss_x = seed | 1;
    kiss_y = seed | 2;
    kiss_z = seed | 4;
    kiss_w = seed | 8;
    kiss_carry = 0;
};

unsigned long rand_kiss() {
    kiss_x = kiss_x * 69069 + 1;
    kiss_y ^= kiss_y << 13;
    kiss_y ^= kiss_y >> 17;
    kiss_y ^= kiss_y << 5;
    kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
    kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
    kiss_z = kiss_w;
    kiss_w = kiss_m;
    kiss_carry = kiss_k >> 30;
    return(kiss_x + kiss_y + kiss_w);
};

void seed_rand(const uint32_t globNeuralId) {
  unsigned long shiftedId;
  shiftedId=globNeuralId;
  seed_rand_kiss(shiftedId<<4);
};

bool onceEveryThousandRand() {
  unsigned long temp1,temp2;
  temp1 = rand_kiss();
  temp2 = (temp1 >>6) & 0x000003FF;
  if (temp2 == 0x3FF) {return (1);} else {return(0);};  
};

