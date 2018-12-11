// DPSNN_chrono.c
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
#include "DPSNN_chrono.h"

#ifdef DALonlyEnvironmentSelected
  double MPI_Wtime() {
    return(1.0);
  };
#endif

void chronoClass::clearChrono() {
  
  
      uint32_t totalSimTime;
      uint32_t startChrono, chronoWindow;
      char *ps_startChrono, *ps_chronoWindow;
      totalSimTime = atoi(getenv("env_totalSimTime_ms"));
      ps_startChrono = getenv("env_startPartialChrono_ms");
      ps_chronoWindow = getenv("env_partialChronoWindow_ms");
      if(ps_startChrono!=NULL)
	startChrono = atoi(ps_startChrono);
      else
	startChrono = 0;
      if(ps_chronoWindow!=NULL)
	chronoWindow = atoi(ps_chronoWindow);
      else
	chronoWindow = 0;

      istogramData = new double [totalSimTime];
      if (istogramData == NULL) {
	printf("ERROR- memory for istogramData could not be allocated\n");
	fflush(stdout);
	exit(0);
      };
  
      startChronoTime = 0.0;
      stopChronoTime = 0.0;
      accumulatedChronoTime = 0.0;
      numLaps = 0;
      minChronoTime = 100000.0;
      maxChronoTime = 0.0;
      meanChronoTime = 0.0;
      varianceChronoTime = 0.0;
      sigmaChronoTime = 0.0;
      coeffOfVariationChronoTime = 0.0;
      M2 = 0.0;
      indexIstogramData = 0;
     
      startPartialChrono_ms = startChrono;
      stopPartialChrono_ms = startChrono + chronoWindow;
}; 

double chronoClass::startChrono() {
  startChronoTime = (double)MPI_Wtime();
  return((double)startChronoTime);
}; 

double chronoClass::stopChrono() {
  double stepChronoTime;
  double deltaChronoTime;

  stopChronoTime = MPI_Wtime();
  stepChronoTime = (double)stopChronoTime - (double)startChronoTime;
  accumulatedChronoTime = (double) stepChronoTime + 
                          (double) accumulatedChronoTime;
  numLaps ++;
  minChronoTime = (stepChronoTime < minChronoTime) ? stepChronoTime : minChronoTime;
  maxChronoTime = (stepChronoTime > maxChronoTime) ? stepChronoTime : maxChronoTime;;
  deltaChronoTime = stepChronoTime - meanChronoTime;
  meanChronoTime += (double)(deltaChronoTime / numLaps);
  M2 += (double)((double)deltaChronoTime * (double)((double)stepChronoTime - (double)meanChronoTime));
  varianceChronoTime = (double)(M2 / (numLaps -1));
  sigmaChronoTime = (double)(sqrt(varianceChronoTime));
  coeffOfVariationChronoTime = (double)(sigmaChronoTime / meanChronoTime);
  return((double)stopChronoTime);
};

double chronoClass::stopChronoIstogram() {
  double stepChronoTime;
  double deltaChronoTime;

  stopChronoTime = MPI_Wtime();
  stepChronoTime = (double)stopChronoTime - (double)startChronoTime;
  if ((indexIstogramData % 1000) == 999)
    stepChronoTime = istogramData[indexIstogramData-1];
  istogramData[indexIstogramData] = (double)stepChronoTime;
  indexIstogramData++;
  accumulatedChronoTime = (double) stepChronoTime + 
                          (double) accumulatedChronoTime;
  numLaps ++;
  minChronoTime = (stepChronoTime < minChronoTime) ? stepChronoTime : minChronoTime;
  maxChronoTime = (stepChronoTime > maxChronoTime) ? stepChronoTime : maxChronoTime;;
  deltaChronoTime = stepChronoTime - meanChronoTime;
  meanChronoTime += (double)(deltaChronoTime / numLaps);
  M2 += (double)((double)deltaChronoTime * (double)((double)stepChronoTime - (double)meanChronoTime));
  varianceChronoTime = (double)(M2 / (numLaps -1));
  sigmaChronoTime = (double)(sqrt(varianceChronoTime));
  coeffOfVariationChronoTime = (double)(sigmaChronoTime / meanChronoTime);
  return((double)stopChronoTime);
};

double chronoClass::stopChronoPartial(uint32_t thisSimTimeStep_ms) {
  double stepChronoTime;
  double deltaChronoTime;

  stopChronoTime = MPI_Wtime();
  if((thisSimTimeStep_ms >= startPartialChrono_ms)&&(thisSimTimeStep_ms < stopPartialChrono_ms))
    {
      stepChronoTime = (double)stopChronoTime - (double)startChronoTime;
      accumulatedChronoTime = (double) stepChronoTime + 
	(double) accumulatedChronoTime;
      numLaps ++;
      minChronoTime = (stepChronoTime < minChronoTime) ? stepChronoTime : minChronoTime;
      maxChronoTime = (stepChronoTime > maxChronoTime) ? stepChronoTime : maxChronoTime;;
      deltaChronoTime = stepChronoTime - meanChronoTime;
      meanChronoTime += (double)(deltaChronoTime / numLaps);
      M2 += (double)((double)deltaChronoTime * (double)((double)stepChronoTime - (double)meanChronoTime));
      varianceChronoTime = (double)(M2 / (numLaps -1));
      sigmaChronoTime = (double)(sqrt(varianceChronoTime));
      coeffOfVariationChronoTime = (double)(sigmaChronoTime / meanChronoTime);
    }
  return((double)stopChronoTime);
};

double chronoClass::getAccumulatedChrono() {
   return((double)accumulatedChronoTime);
}; 

uint32_t chronoClass::getNumLaps() {
   return(numLaps);
};

void chronoClass::clearAndStartChrono() {
  clearChrono(); startChrono();
};

double chronoClass::getMinChrono() {
   return((double)minChronoTime);
}; 
double chronoClass::getMaxChrono() {
   return((double)maxChronoTime);
}; 
double chronoClass::getMeanChrono() {
   return((double)meanChronoTime);
}; 
double chronoClass::getVarianceChrono() {
   return((double)varianceChronoTime);
}; 
double chronoClass::getSigmaChrono() {
   return((double)sigmaChronoTime);
}; 
double chronoClass::getCoeffOfVariationChrono() {
   return((double)coeffOfVariationChronoTime);
}; 
double *chronoClass::getIstogramData() {
   return(&istogramData[0]);
}
uint32_t chronoClass::getstartPartialChrono_ms() {
   return((uint32_t)startPartialChrono_ms);
};
uint32_t chronoClass::getstopPartialChrono_ms() {
   return((uint32_t)stopPartialChrono_ms);
};
void chronoClass::printChrono(const char *chronoName,uint32_t loc_h,bool initFile){
  FILE *fp_output;
  char reportName[80];

  sprintf(reportName,"SummaryStatMeasures_h%d.stat",loc_h);
  if(initFile){
      fp_output=fopen(reportName,"w");
      fprintf(fp_output,"Summary of statistical measures for process H=%d\n",loc_h);
    }
  else
    fp_output=fopen(reportName,"a");

  fprintf(fp_output,"CHRONO-h=%.3d- %30.30s = %-11.9f min = %-11.9f max = %-11.9f mean = %-11.9f variance = %-11.9f sigma = %-11.9f coeffOfVariation = %-11.9f samples = %-11d\n",loc_h,chronoName,accumulatedChronoTime,minChronoTime,maxChronoTime,meanChronoTime,varianceChronoTime,sigmaChronoTime,coeffOfVariationChronoTime,numLaps);
  fclose(fp_output);
}

void chronoClass::printIstogramData(const char *dataFile,uint32_t loc_h,uint32_t totalSimTime_ms){
  FILE *fp_output;
  char reportName[80];

  sprintf(reportName,"%s_h%d.stat",dataFile,loc_h);
  fp_output=fopen(reportName,"w");
  for(uint32_t i=0;i<totalSimTime_ms;i++)
    fprintf(fp_output,"%.9f\n",istogramData[i]);

}
