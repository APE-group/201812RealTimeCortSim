// DPSNN_stopWatch.h
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
 
#ifndef DPSNN_stopWatchIncluded
#define DPSNN_stopWatchIncluded

enum stopWatchStatusEnum { 
  maxSimTime_ms_reached, x999_ms_reached, new_s, same_s};

class stopWatchClass 
{
 public:
  void programSimTime_ms(const uint32_t maxFrames_value);
  void advanceSimTime_ms();
  stopWatchStatusEnum getStatus();
  int32_t getSimTime_ms();
 private:
  int32_t maxSimTime_ms;;
  int32_t ms;
};

#endif //multiple inclusion guard
