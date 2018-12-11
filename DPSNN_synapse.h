// DPSNN_synapse.h
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

#ifndef DPSNN_synapseIncluded
#define DPSNN_synapseIncluded

#include "DPSNN_dataStructDims.h"
#include "DPSNN_chrono.h"
#include "DPSNN_parameters.h"

class forwardSynClass {
 public:
  //global identity of presynaptic neuron
  uint32_t pre_glob_n; 
  //global identity of post synaptic neuron
  uint32_t post_glob_n;
  //in ms 
  uint8_t delay; 
  //kind of presynaptic neuron
  uint8_t preSynNeuralKind; 
  //weight
  weightType weight;
};

class instrumentedForwardSyn: public forwardSynClass {
 public:
  void report(const forwardSynClass synapseToBeReported, FILE *fp_output, uint32_t locN) 
  {
    fprintf(fp_output,"pre_glob_n=%d, post_glob_n=%d, post_loc_h=%d, delay=%d, kind=%d \n",
	    synapseToBeReported.pre_glob_n, 
	    synapseToBeReported.post_glob_n,
	    synapseToBeReported.post_glob_n/locN,
	    synapseToBeReported.delay,
	    synapseToBeReported.preSynNeuralKind);
  };
};

class synapseClass {
 public:
  uint32_t pre_glob_n; //global identity of presynaptic neuron
  uint32_t post_glob_n; //global identity of post synaptic neuron
  uint8_t delay;//in ms
  uint8_t preSynNeuralKind; // kind of presynaptic neuron
  weightType weight;
#if defined(makeActiveLTD) || defined (makeActiveLTP)
  float timeDerivative;
  int32_t lastActivationTime;
#endif
};

class instrumentedSynapse: public synapseClass {
 public:
  void report(const synapseClass synapseToBeReported, FILE *fp_output,
	      uint32_t locN, float factorWeightType_2_Float) {
    float w;
#ifdef LIFCAneuron
    w=(float)(synapseToBeReported.weight * factorWeightType_2_Float);
    // fprintf(fp_output,"pre_glob_n=%d, post_glob_n=%d, post_loc_h=%d, delay=%d, kind=%d, weight=%f \n",
    fprintf(fp_output,"%d %d %d %d %d %f ",
	    synapseToBeReported.pre_glob_n, 
	    synapseToBeReported.post_glob_n,
            synapseToBeReported.post_glob_n/locN,
	    synapseToBeReported.delay,
	    synapseToBeReported.preSynNeuralKind,
	    w);
#if defined(makeActiveLTD) || defined (makeActiveLTP)
    fprintf(fp_output,"%f",
	    synapseToBeReported.timeDerivative);
#endif
    fprintf(fp_output,"\n");
#else
    fprintf(fp_output,"pre_glob_n=%d, post_glob_n=%d, post_loc_h=%d, delay=%d, kind=%d, weight=%f, timeDer=%f \n",
	    synapseToBeReported.pre_glob_n, 
	    synapseToBeReported.post_glob_n,
            synapseToBeReported.post_glob_n/locN,
	    synapseToBeReported.delay,
	    synapseToBeReported.preSynNeuralKind,
	    synapseToBeReported.weight,
	    synapseToBeReported.timeDerivative);
#endif
    };
};

class synapticDistributionClass {
 public:
  uint32_t synCount[DSD__maxGlobH];	
  uint32_t synOffset[DSD__maxGlobH];
};

class synapticDelayDistributionClass {
 public:
  uint32_t synCount[DSD__maxD][DSD__maxGlobH];	
  uint32_t synOffset[DSD__maxD][DSD__maxGlobH];
};

//critical assumption synapses are ordered by source neuron
//the search function must be reset at each timeStep
class synapticSearchClass { 
 public:
  void reset(const uint32_t synOffset, 
	     const uint32_t synCount, const int resetTime_ms);
  void setSynIndex(const uint32_t synOffsets);
  //get returns a synapse
  bool findAndGetNextSyn(const uint32_t sourceNeuron,			 
			 synapseClass *returnSynapse,
			 const synapseClass *backwardSynList,
			 uint32_t *backwardSynGlobalOffset,
			 const int thisTimeStep_ms);
 private:
  uint32_t synIndex;
  uint32_t maxSynIndex;
  uint32_t lastSearchedSourceNeuron;
  int32_t lastSearchResetTime_ms;

};

class synOffsetInSynList {
 public:
  uint32_t offsetByDelay[DSD__maxD];
};
#endif //inclusion guard
