#ifndef _MERSENNE_TWISTER
#define _MERSENNE_TWISTER

// N.C. Cruz: This is a straightforward C adaptation of the C++ code provided in:
// "The Mersenne Twister Pseudo Random Number Generator"
// posted August 19, 2014 by Stephan Brumme, who is the real and original author of the code 
// https://create.stephan-brumme.com/mersenne-twister/

// Adapted from: https://github.com/stbrumme/mersenne-twister/blob/master/mersenne.h

#include <stdint.h>

/// state size
enum   { SizeState = 624 };

typedef struct _mersennseTwister{
  	/// internal state
  	uint32_t state[SizeState];
  	/// offset of next state's word
  	int next;
} Mersenne_Twister;

void Create_MersenneTwister(Mersenne_Twister* rng, uint32_t seed); //Expected default: 5489

uint32_t getRandomNumber(Mersenne_Twister* rng); // Returns a random 32 bit number

#endif
