// DPSNN_chrono.h
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

#ifndef DPSNN_chronoIncluded
#define DPSNN_chronoIncluded
 
#ifdef DALonlyEnvironmentSelected
double MPI_Wtime();
#endif

class chronoClass { 
  private:
    double startChronoTime; 
    double stopChronoTime;
    double accumulatedChronoTime; 
    uint32_t numLaps;
    double minChronoTime;
    double maxChronoTime;
    double meanChronoTime;
    double varianceChronoTime;
    double sigmaChronoTime;
    double coeffOfVariationChronoTime;
    double M2;
    double *istogramData;
    uint32_t indexIstogramData;
    uint32_t startPartialChrono_ms;
    uint32_t stopPartialChrono_ms;
 public:
    void clearChrono();
    double startChrono();
    double stopChrono();
    double stopChronoIstogram();
    double stopChronoPartial(uint32_t thisSimTimeStep_ms);
    double getAccumulatedChrono();
    uint32_t getNumLaps();
    void clearAndStartChrono();
    double getMinChrono();
    double getMaxChrono();
    double getMeanChrono();
    double getVarianceChrono();
    double getSigmaChrono();
    double getCoeffOfVariationChrono();
    double *getIstogramData();
    uint32_t getstartPartialChrono_ms();
    uint32_t getstopPartialChrono_ms();
    void printChrono(const char *chronoName,uint32_t loc_h,bool initFile);
    void printIstogramData(const char *dataFile,uint32_t loc_h,uint32_t totalSimTime_ms);
};

#endif //multiple inclusion guard
