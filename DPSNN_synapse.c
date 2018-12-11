// DPSNN_synapse.c
// DPSNN_*.* Distribution/Parallelization 
// performed by Pier Stanislao Paolucci (Roma, Italy project start date 2011),
// starting from the sequential code
// SPNET.cpp: Spiking neural network with axonal conduction delays and STDP
// Created by Eugene M. Izhikevich, May 17, 2004, San Diego, CA
  
// reference paper, named [IzhPol] in the following 
// Eugene M. Izhikevich "Polychronization: Computation with Spikes" 
// Neural Computation 18, 245-282 (2006)
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "DPSNN_debug.h"
#include "DPSNN_random.h"

#include "DPSNN_environmentSelection.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"
#include "DPSNN_localNet.h"

//critical assumption synapses are ordered by source neuron

void synapticSearchClass::reset(const uint32_t synOffset, 
const uint32_t synCount, const int resetTime_ms_initValue) 
{
  lastSearchResetTime_ms = resetTime_ms_initValue;
  synIndex = synOffset;
  maxSynIndex = synOffset + synCount;
  lastSearchedSourceNeuron = 0;
 
  DPSNNverboseStart(false,resetTime_ms_initValue,0);
    printf(
     "SynSearch Reset Start Syn Index %d Max Syn Index %d\n",
     synIndex, maxSynIndex);
    fflush(stdout);
  DPSNNverboseEnd();

};

void synapticSearchClass::setSynIndex(const uint32_t synOffset) 
{
  synIndex = synOffset;
};

bool synapticSearchClass::findAndGetNextSyn(
  const uint32_t sourceNeuron_id,
  synapseClass *returnSyn, 
  const synapseClass *synList,
  uint32_t *synGlobalOffset,
  const int thisTimeStep_ms) 
{
  uint32_t maxSourceNeuronFoundDuringThisSearch;

  if(thisTimeStep_ms != lastSearchResetTime_ms) {
    printf("ERROR synapticSearch has not been reset at time %d\n",
	   thisTimeStep_ms);fflush(stdout);exit(0);
  };

  if(sourceNeuron_id < lastSearchedSourceNeuron) {
    printf(
    "ERROR findAndGetNextSyn searches backward: for a sourceNeuron %d < lastSearchedNeuron %d\n",
	   sourceNeuron_id, lastSearchedSourceNeuron);
        fflush(stdout);
        exit(0);
  } else {
    //it is legal to search for the synapse 
    //in increasing order of source neuron
    lastSearchedSourceNeuron = sourceNeuron_id;
  };

  maxSourceNeuronFoundDuringThisSearch = 0;

  DPSNNverboseStart(false,thisTimeStep_ms,0);
    printf("findGetNextSyn: START source Neu=%d starting from synIndex=%d\n",sourceNeuron_id, synIndex);
  DPSNNverboseEnd();

  do 
  { 
    * returnSyn = synList[synIndex];
    * synGlobalOffset = synIndex;

    DPSNNverboseStart(false,thisTimeStep_ms,0);
          printf(
            "findGetNextSyn: search at synIndex=%d for source neu=%d \n",
            synIndex, sourceNeuron_id); 
          printf(
            "findGetNextSyn: synList[%d].pre_glob_n = %d \n",
            synIndex, synList[synIndex].pre_glob_n);
          printf(
            "findGetNextSyn: returnSyn->pre_glob_n = %d, sourceNeuron=%d\n",
            returnSyn->pre_glob_n, sourceNeuron_id);
    DPSNNverboseEnd();

    //this way we are sure that synList is ordered by source neuron
    if ((returnSyn->pre_glob_n) < maxSourceNeuronFoundDuringThisSearch) {
      printf("ERROR findAndGetNextSyn synList non ordered - source %d < max neu already found %d\n",
	     (returnSyn->pre_glob_n), maxSourceNeuronFoundDuringThisSearch);
        fflush(stdout);
        exit(0);
    } else {
      //this way we are doimg a check that synList is ordered by source neuron
      maxSourceNeuronFoundDuringThisSearch =  returnSyn->pre_glob_n;
    };

    if ((returnSyn->pre_glob_n) < sourceNeuron_id) {
      maxSourceNeuronFoundDuringThisSearch = (returnSyn->pre_glob_n);
      synIndex++;
    }else{
      if ((returnSyn->pre_glob_n) == sourceNeuron_id) {
         DPSNNverboseStart(false,thisTimeStep_ms,0);
         printf(
         "findAndGetNextSyn: found at synIndex=%d for source neu=%d \n",
            synIndex, sourceNeuron_id); 
         fflush(stdout);
        DPSNNverboseEnd();

        synIndex ++;
        return(true);
      } else {
        if ((returnSyn->pre_glob_n) > sourceNeuron_id) {
          maxSourceNeuronFoundDuringThisSearch = returnSyn->pre_glob_n;
          return(false);
        }else{
          printf("ERROR findAndGetNextSyn - A should never arrive here\n");
          fflush(stdout);
          exit(0);
          return(false);
	}
      }
    };
  } while(synIndex < maxSynIndex);
  return(false);
};

