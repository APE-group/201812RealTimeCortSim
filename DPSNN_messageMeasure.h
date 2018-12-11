// DPSNN_messageMeasure.h
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

#include "DPSNN_environmentSelection.h" 
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_spike.h"

#ifndef DPSNN_messageMeasureIncluded
#define DPSNN_messageMeasureIncluded
 
class msgMeasureClass { 
  private:
  //double startMsgMeasure; 
  //double stopMsgMeasure;
    uint32_t accumulatedMsgSize; 
    uint32_t numLaps;
    uint32_t minMsgSize;
    uint32_t maxMsgSize;
    double meanMsgSize;
    double varianceMsgSize;
    double sigmaMsgSize;
    double coeffOfVariationMsgSize;
    double M2;
 public:
    void clearMsgMeasure();
    //double startMsgMeasure();
    //double stopMsgMeasure();
    //void clearAndStartMsgMeasure();
    void accumulateMsgSize(uint32_t stepMsgSize);
    uint32_t getAccumulatedMsgSize();
    uint32_t getNumLaps();
    uint32_t getMinMsgSize();
    uint32_t getMaxMsgSize();
    double getMeanMsgSize();
    double getVarianceMsgSize();
    double getSigmaMsgSize();
    double getCoeffOfVariationMsgSize();
    void printMsgSize(const char *msgMeasureName,uint32_t loc_h,bool initFile);
};

#endif //multiple inclusion guard
