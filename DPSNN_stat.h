// DPSNN_stat.h
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

#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"

#ifndef DPSNN_stat
#define DPSNN_stat

#define DEF_stat_fileNameSize 80

class statBasicClass { 
  public:
    DPSNN_parameters lnp_par;
    char prefixName[DEF_stat_fileNameSize];
    char fName[DEF_stat_fileNameSize*2];
    FILE *fp;
    bool writeDisable;
    bool prepCalled;
    bool fopenStatus;
    void prep(const DPSNN_parameters lnp_par_initValue,
	      const char prefixName_value[]);
    void openFile(const uint32_t sec, const int moduloSec);
    void closeFile(const uint32_t sec, const int moduloSec);
    void test_prepCalled();
    void test_fopenStatus(const char methodName[], 
			  const bool stat_booleanValue, 
			  const int sec, const int moduloSec );
    void debugPrintStatStatus(const int sec, const int moduloSec);
};

class statParametersClass : public statBasicClass {
  public:
    void write(const uint32_t sec, const int moduloSec);
};

class statSpikingRatesClass : public statBasicClass {
  public:
  void write(const uint32_t sec, const int moduloSec,
	     const int N_firings,
	     const int N_thalamicInputs);
};


class statSpikingRatesPerPopClass : public statBasicClass {
  public:
    void write(const uint32_t sec,
	       const uint32_t milliSec, 
	       const uint32_t moduloSec, 
	       const uint32_t ratesSampling,
	       uint32_t *N_firingsPerPop,
	       const uint32_t *neuSubPopCount,
	       const uint32_t neuSubPopTotal);
};

class statSpikesClass : public statBasicClass  {
  public:
    void write(const uint32_t sec, const int moduloSec,
	       const int N_firings, const int thisMs,
	       const int data[],const double emissionTime[]);
};

class  statSpikeCountPerMsClass: public statBasicClass  {
  public:
    void write(const uint32_t sec, const int moduloSec, 
	       const int thisMs,
	       const int numberOfSpikesInMs);
};

class  statThalamicInputClass: public statBasicClass  {
  public:
    void write(const uint32_t sec, const int moduloSec, 
	       const int thisMs,
	       const int numberOfThalInputInMs);
};

class statFloatEvery_msClass : public statBasicClass {
  public:
    void write(const uint32_t sec, const int moduloSec, const uint32_t t_ms, 
               const int neuron_id, const float value);
};

class stat3FloatsEvery_msClass : public statBasicClass {
  public:
    void write(const uint32_t sec, const int moduloSec, const uint32_t t_ms, 
               const int neuron_id, 
	       const float I, const float v, const float u);
};

class statInputCurrentClass : public statFloatEvery_msClass { };
class stat_ivu_Class : public stat3FloatsEvery_msClass { };
class stat_LIFCAivu_Class : public statBasicClass { 
  public:
    void write(const uint32_t sec, const int moduloSec, const uint32_t t_ms, 
               const double currentTime, const int neuron_id, 
	       const double I, const double v, const double u);
};

class statSynClass : public statBasicClass {
  public:
  void write(const uint32_t sec, const int moduloSec,
	     const int this_ms, const int synIndex,
	     const int preSynNeu, const int postSynNeu,
	     const int synDelay,
	     const float synWeight, const float synDeriv);
};

class statSynPeriodicProbeClass : public statBasicClass {
 public:
 void write(const uint32_t this_s, const int moduloSec, const synapseClass synToBeReported);
};

class statSTDPeventClass : public statBasicClass {
 public:
 void write(const uint32_t this_s, const int moduloSec,
	    const uint32_t this_ms, const int preNeu,
	    const int postNeu,
	    const float contribute);
};

class statMessageTrafficClass : public statBasicClass {
 public:
 void write(const uint32_t sec, const int moduloSec,
	    const uint32_t t_ms, uint32_t *forwardAxonalSpikesCount);
};

class statClass {
 private:
  DPSNN_parameters lnp_par;
  uint32_t currentStat_sec;
 public:
  statParametersClass parameters;
  statSpikingRatesClass spikingRates;
  statSpikingRatesPerPopClass spikingRatesPerPop;
  statSpikeCountPerMsClass spikePerMs;
  statSpikesClass spikes;
  statThalamicInputClass thalamicInput;
  statInputCurrentClass inputCurrent;
  stat_ivu_Class ivu0,ivu1,ivu2,ivu999;
  stat_LIFCAivu_Class LIFCAivu0,LIFCAivu1,LIFCAivu2,LIFCAivu999;
  statSynClass syn;
  statSTDPeventClass statSTDPevent;
  statSynPeriodicProbeClass statSynPeriodicProbeA;
  statSynPeriodicProbeClass statSynPeriodicProbeB;
  statMessageTrafficClass messageTraffic;

  //statSynClass syn;

  void prep(DPSNN_parameters lnp_par_initValue);
  void openFiles_forThisSecond(const uint32_t sec, const int moduloSec);
  void closeFiles_forThisSecond(const uint32_t sec, const int moduloSec);
  void check_fileClosure(const int sec, const int moduloSec);
};

#endif //multiple inclusion guard
