// DPSNN_stopWatch.c
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

#include "DPSNN_stopWatch.h"

void stopWatchClass::programSimTime_ms(const uint32_t maxSimTime_ms_initValue) {
  maxSimTime_ms = maxSimTime_ms_initValue;
  ms=0;
};

int32_t stopWatchClass::getSimTime_ms() {
  return(ms);
};

void stopWatchClass::advanceSimTime_ms() {
  ms++; 
};

stopWatchStatusEnum stopWatchClass::getStatus() {
  //increment the "millisecond clock hand"
  if(ms > maxSimTime_ms) {
      printf("ERROR programmed max sim time ms exceeded\n");
      exit(0);
      fflush(stdout);
  } else if (ms == maxSimTime_ms) {
     return(maxSimTime_ms_reached);
  } else if ((ms % 1000) == 999) {
    return(x999_ms_reached);
  } else if((ms % 1000) == 0) {
     return(new_s);
  } else {
    return(same_s);
  };
};
