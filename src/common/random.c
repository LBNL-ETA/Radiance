#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Fast random number generator
 * for systems without drand48(), i.e., Windows
 *
 * Algorithm: xorshift64* (xorshift 64-bit Star)
 * Author: Sebastiano Vigna (vigna@acm.org)
 * Based on the original xorshift work by George Marsaglia.
 *
 * Reference:
 * S. Vigna, "An experimental exploration of Marsaglia's xorshift
 *	generators, scrambled",
 * ACM Transactions on Mathematical Software 42(4), 2016.
 *
 * Source: https://prng.di.unimi.it/
 */

#include <stdint.h>
#include "random.h"

static uint64_t		rand_state = 0x853c49e6748fea9bULL;

/* SplitMix64 for seeding */
static uint64_t
splitmix64(uint64_t x)
{
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

void
srandom(unsigned long s)
{
	rand_state = splitmix64((uint64_t)s);
	if (!rand_state)
		rand_state = 0x853c49e6748fea9bULL;
}

static uint64_t
xorshift64start(void)
{
	uint64_t x = rand_state;
	x ^= x >> 12;
	x ^= x << 25;
	x ^= x >> 27;
	rand_state = x;
	return x * 2685821657736338717ULL;
}

long
random(void)
{
	return (long)(xorshift64start() >> 33);
}

double
frandom(void)
{
	return (double)(xorshift64start() >> 11) * (1.0/9007199254740992.0);
}
