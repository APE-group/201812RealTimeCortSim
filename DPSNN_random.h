// DPSNN_random.h
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

#ifndef DPSNN_randomIncluded
#define DPSNN_randomIncluded

#define getRandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1

// Implementation based on the ISO C standard, here extended to 32 bits
uint32_t getRandom_r (uint32_t *seed, uint32_t max1);

#endif //multiple inclusion guard
