#include "mersenne.h"

// N.C. Cruz: This is a straightforward C adaptation of the C++ code provided in:
// "The Mersenne Twister Pseudo Random Number Generator"
// posted August 19, 2014 by Stephan Brumme, who is the real and original author of the code 
// https://create.stephan-brumme.com/mersenne-twister/

// Adapted from: https://github.com/stbrumme/mersenne-twister/blob/master/mersenne.cpp

static void twist(Mersenne_Twister* rng){ // Create new state (based on old one)
	const int M = 397;
	const int FirstHalf = SizeState - M;

	// first 624-397=227 words
	int i;
	for (i = 0; i < FirstHalf; i++){
    	uint32_t bits = (rng->state[i] & 0x80000000) | (rng->state[i + 1] & 0x7fffffff);
    	rng->state[i] = rng->state[i + M]         ^ (bits >> 1) ^ ((bits & 1) * 0x9908b0df);
  	}
  	// remaining words (except the very last one)
  	for ( ; i < SizeState - 1; i++){
    	uint32_t bits = (rng->state[i] & 0x80000000) | (rng->state[i + 1] & 0x7fffffff);
    	rng->state[i] = rng->state[i - FirstHalf] ^ (bits >> 1) ^ ((bits & 1) * 0x9908b0df);
  	}

  	// last word is computed pretty much the same way, but i + 1 must wrap around to 0
  	uint32_t bits = (rng->state[i] & 0x80000000) | (rng->state[0] & 0x7fffffff);
  	rng->state[i] = rng->state[M - 1] ^ (bits >> 1) ^ ((bits & 1) * 0x9908b0df);

	// word used for next random number
  	rng->next = 0;
}

void Create_MersenneTwister(Mersenne_Twister* rng, uint32_t seed){ //Expected default: 5489
	rng->next = 0;
	rng->state[0] = seed;
	for (int i = 1; i < SizeState; i++){
		rng->state[i] = 1812433253UL * (rng->state[i-1] ^ (rng->state[i-1] >> 30)) + i;
	}
	// let's twist'n'shout ...
  	twist(rng);
}

uint32_t getRandomNumber(Mersenne_Twister* rng){
	// compute new state ?
	if (rng->next >= SizeState){
    	twist(rng);
	}

  	// shuffle bits around
  	uint32_t x = rng->state[rng->next++];
  	x ^=  x >> 11;
  	x ^= (x <<  7) & 0x9d2c5680;
  	x ^= (x << 15) & 0xefc60000;
  	x ^=  x >> 18;
  	return x;
}
