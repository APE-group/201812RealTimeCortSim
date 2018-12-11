// DPSNN_messageMeasure.c
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
#include "DPSNN_debug.h"
#include "DPSNN_environmentSelection.h"
#include "DPSNN_messageMeasure.h"


void msgMeasureClass::clearMsgMeasure() {
  //startMsgMeasure = 0.0;
  //stopMsgMeasure = 0.0;
  accumulatedMsgSize = 0;
  numLaps = 0;
  minMsgSize = 100000;
  maxMsgSize = 0;
  meanMsgSize = 0.0;
  varianceMsgSize = 0.0;
  sigmaMsgSize = 0.0;
  coeffOfVariationMsgSize = 0.0;
  M2 = 0.0;     
}; 

void msgMeasureClass::accumulateMsgSize(uint32_t stepMsgSize) {
  double deltaMsgSize;
#ifdef LIFCAneuron  
  stepMsgSize = stepMsgSize*3*sizeof(stepMsgSize);
#else
  stepMsgSize = stepMsgSize*2*sizeof(stepMsgSize);
#endif
  accumulatedMsgSize += stepMsgSize;
  numLaps ++;
  minMsgSize = (stepMsgSize < minMsgSize) ? stepMsgSize : minMsgSize;
  maxMsgSize = (stepMsgSize > maxMsgSize) ? stepMsgSize : maxMsgSize;
  deltaMsgSize = (double)stepMsgSize - meanMsgSize;
  meanMsgSize += (double)(deltaMsgSize / numLaps);
  M2 += (double)((double)deltaMsgSize * (double)((double)stepMsgSize - (double)meanMsgSize));
  varianceMsgSize = (double)(M2 / (numLaps -1));
  sigmaMsgSize = (double)(sqrt(varianceMsgSize));
  coeffOfVariationMsgSize = (double)(sigmaMsgSize / meanMsgSize);
};

uint32_t msgMeasureClass::getAccumulatedMsgSize() {
  return(accumulatedMsgSize);
}; 
uint32_t msgMeasureClass::getNumLaps() {
  return(numLaps);
};
uint32_t msgMeasureClass::getMinMsgSize() {
  return(minMsgSize);
}; 
uint32_t msgMeasureClass::getMaxMsgSize() {
  return(maxMsgSize);
}; 
double msgMeasureClass::getMeanMsgSize() {
  return((double)meanMsgSize);
}; 
double msgMeasureClass::getVarianceMsgSize() {
  return((double)varianceMsgSize);
}; 
double msgMeasureClass::getSigmaMsgSize() {
  return((double)sigmaMsgSize);
}; 
double msgMeasureClass::getCoeffOfVariationMsgSize() {
  return((double)coeffOfVariationMsgSize);
}; 
void msgMeasureClass::printMsgSize(const char *msgMeasureName,uint32_t loc_h,bool initFile){
  FILE *fp_output;
  char reportName[80];

  sprintf(reportName,"SummaryStatMeasures_h%d.stat",loc_h);
  if(initFile){
      fp_output=fopen(reportName,"w");
      fprintf(fp_output,"Summary of statistical measures for process H=%d\n",loc_h);
    }
  else
    fp_output=fopen(reportName,"a");

  fprintf(fp_output,"MESSAG-h=%.3d- %24.24s(byte) = %-11d min = %-11d max = %-11d mean = %-11.3f variance = %-11.3f sigma = %-11.3f coeffOfVariation = %-11.3f samples = %-11d\n",loc_h,msgMeasureName,accumulatedMsgSize,minMsgSize,maxMsgSize,meanMsgSize,varianceMsgSize,sigmaMsgSize,coeffOfVariationMsgSize,numLaps);
  fclose(fp_output);
}

