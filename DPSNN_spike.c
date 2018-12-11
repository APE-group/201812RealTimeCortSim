// DPSNN_spike.c
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
#include "DPSNN_random.h"

#include "DPSNN_environmentSelection.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"
#include "DPSNN_neuron.h"
#include "DPSNN_spike.h"

void synapticSpikeSchedulerClass::reset(
  const uint32_t owner_glob_n_initValue,
  const uint32_t owner_loc_n_initValue,
  const uint32_t owner_loc_h_initValue) 
{
  owner_glob_n = owner_glob_n_initValue; 
  owner_loc_n  = owner_loc_n_initValue;
  owner_loc_h  = owner_loc_h_initValue;

  if(DSD__maxSimultaneousSpikesOnSameTarget < 4) {
    printf("ERROR DSD__maxSimultaneousSpikesOnSameTarget < 4\n");
    fflush(stdout);exit(0);        
  };

  lastReadTime=0;
  readCircTime=0;
  fillCountByCircTime=0;
  resetCalled=true;
}

void synapticSpikeSchedulerClass::insertSynSpike(
  const uint32_t target_owner_glob_n, 
  const synapticSpikeClass synSpikeToInsert) {
  int32_t writeTime, writeCircTime, fillCountAtWriteCircTime;

  if(resetCalled != true) {
  printf("ERROR miss synSpikeSched.reset() call on glob_n=%d loc_n=%d loc_h=%d\n",
	   owner_glob_n, owner_loc_n, owner_loc_h);
    fflush(stdout);exit(0);        
  };

  if(target_owner_glob_n != owner_glob_n) {
    printf("ERROR wrong target neuron, this is glob_n %d while target_owner_glob_n=%d \n",
	   owner_glob_n,target_owner_glob_n);
    fflush(stdout);exit(0);        
  } 

  writeTime = synSpikeToInsert.synDelayPlusOriginalEmissionTime;

  if((writeTime - lastReadTime) >= 2) {
    printf("ERROR: the event is too much in the future\n");
    fflush(stdout);exit(0);
  }

  if(synSpikeToInsert.synIndex >= (DSD__maxBackwardLocSyn)) {
    printf(
    "ERROR: insertSynSpike synSpikeToInsert.synIndex exceeds DSD__...\n");
    fflush(stdout);exit(0);
  }
  
  //writeCircTime = writeTime & DSD__circBufferBitwiseMaskD;
  writeCircTime = 0;
  fillCountAtWriteCircTime = fillCountByCircTime;

  DPSNNverboseStart(false,writeTime,0);
    printf("insertSynSpike-A futActT %d ms, synIdx=%d, n=%d done at t %d\n",
	   writeTime, synSpikeToInsert.synIndex, owner_glob_n, lastReadTime);
  DPSNNverboseEnd();

  DPSNNverboseStart(false,writeTime,0);
    printf(
    "insertSynSpike-B circT=%d, fillCount=%d\n",
    writeCircTime, fillCountAtWriteCircTime);
  DPSNNverboseEnd();

  circBuffer[fillCountAtWriteCircTime] = synSpikeToInsert.synIndex;

  fillCountByCircTime++;

  if(fillCountByCircTime >=
     DSD__maxSimultaneousSpikesOnSameTarget) {
    printf("ERROR/IMPROVEMENT NEEDED: too much %d simult spikes on same target neu=%d at time %d\n", 
	   fillCountByCircTime, owner_glob_n, writeTime);
    printf("ERROR/IMPROVEMENT NEEDED: synIndex=%d synTime %f\n", 
	   synSpikeToInsert.synIndex, synSpikeToInsert.synDelayPlusOriginalEmissionTime);
    fflush(stdout);exit(0);
  }
};

bool synapticSpikeSchedulerClass::extractSynSpike( 
       synapticSpikeClass * pReturnSynSpike) 
{
  synapticSpikeClass tempSynSpike;

 if(resetCalled != true) {
    printf("ERROR missing synapticSpikeSchedulet.reset() call\n");
    fflush(stdout);exit(0);        
  };

 DPSNNverboseStart(false,lastReadTime,0);
      printf("extractSynSpike-A %d ms started for neu=%d \n",
	     lastReadTime, owner_glob_n);
  DPSNNverboseEnd();

  //readCircTime = lastReadTime & DSD__circBufferBitwiseMaskD;
  readCircTime = 0;
  if(fillCountByCircTime == 0) {
    DPSNNverboseStart(false,lastReadTime,0);
      printf("extractSynSpike-BN t%d ms cT=%d neu=%d DO NOT return a spike\n",
	 lastReadTime, readCircTime, owner_glob_n);
    DPSNNverboseEnd();
    return(false);
  }else{ 
     fillCountByCircTime--;
     tempSynSpike.synDelayPlusOriginalEmissionTime = lastReadTime;
     tempSynSpike.synIndex = circBuffer[fillCountByCircTime];
     DPSNNverboseStart(false,lastReadTime,0);
     printf(
   "extrSynSpike-BY t%d cT%d will RETURN synSpike:actT%f synIdx=%d to neu=%d \n",
     lastReadTime, readCircTime, tempSynSpike.synDelayPlusOriginalEmissionTime, 
     tempSynSpike.synIndex, owner_glob_n);
     DPSNNverboseEnd();

     if(tempSynSpike.synDelayPlusOriginalEmissionTime !=
	lastReadTime) {
       printf(
	 "ERROR: at cT=%d extractSynSpike time in extracted synSpike %f != lastReadTime %d\n",
	 readCircTime, tempSynSpike.synDelayPlusOriginalEmissionTime,
	 lastReadTime);
       fflush(stdout);exit(0);
     }

     if(tempSynSpike.synIndex >= (DSD__maxBackwardLocSyn)) {
       printf(
     "ERROR: extractSynSpike synIndex exceeds DSD__...\n");
       fflush(stdout);exit(0);
     } 

     pReturnSynSpike->synDelayPlusOriginalEmissionTime = 
       tempSynSpike.synDelayPlusOriginalEmissionTime;
      pReturnSynSpike->synIndex = 
       tempSynSpike.synIndex;
    return(true);
  }
};

bool synapticSpikeSchedulerClass::mustBeEmpty(
  const int32_t readTime) 
{
  /* int32_t readCircTime; */
  if(readTime != lastReadTime) {
    printf("ERROR readTime %d != lastReadTime %d on spike queue of n%d\n",
	   readTime,lastReadTime,owner_glob_n);
    fflush(stdout);exit(0);
  };
  
  //readCircTime = lastReadTime & DSD__circBufferBitwiseMaskD;
  readCircTime = 0;
  if(fillCountByCircTime != 0) {
    printf("ERROR circular event buffer non empty at time %d on n%d\n",
	   lastReadTime, owner_glob_n);
    fflush(stdout);exit(0);
  } else {
    return(true);
  };
 return(false);
};

int32_t synapticSpikeSchedulerClass::getReadTime() 
{
 return(lastReadTime);
};

void synapticSpikeSchedulerClass::setReadTime(
  const int readTime_initValue) 
{
  lastReadTime = readTime_initValue;
};
						
// axonalSpikeSchedulerClass
void axonalSpikeSchedulerClass::reset()
{
  int32_t d;
  circBufferLen=DSD__maxD;
  readIndex=-1;
  nextFreeIndex=0;
  for(d=0;d<DSD__maxD;d++)
  {
    circBuffer[d].count=0;
    circBuffer[d].expectedCount=0;
    //EPA-PSP: added next line 2015-03-09
    circBuffer[d].list = (axonalSpikeDataOnlyClass *)NULL;
  };
  resetCalled=true;
}

void axonalSpikeSchedulerClass::insertAxSpike(const uint32_t source_h,
					      const backwardAxonalSpikesClass *tempAxSpike,
					      const uint32_t thisSimTimeStep_ms)
{
  uint32_t i;
  uint32_t axSpikeCount;

  if(resetCalled != true) {
    printf("ERROR miss axonalSpikeScheduler.reset() call \n");
    fflush(stdout);exit(0);        
  };
 
  DPSNNverboseStart(false,thisSimTimeStep_ms,0);
    printf("insertAxSpike PRE --- at ms=%d readIndex=%d nextFreeIndex=%d \n",
	   thisSimTimeStep_ms,readIndex,nextFreeIndex);
  DPSNNverboseEnd();

  axSpikeCount = tempAxSpike[source_h].count;
  circBuffer[nextFreeIndex].count = tempAxSpike[source_h].count;
  circBuffer[nextFreeIndex].expectedCount = tempAxSpike[source_h].expectedCount;
  for(i=0;i<axSpikeCount;i++)
    circBuffer[nextFreeIndex].list[i] = tempAxSpike[source_h].list[i];
    
  readIndex=nextFreeIndex;
  nextFreeIndex = (nextFreeIndex + 1 ) % circBufferLen;

  DPSNNverboseStart(false,thisSimTimeStep_ms,0);
    printf("insertAxSpike POST --- at ms=%d readIndex=%d nextFreeIndex=%d \n",
	   thisSimTimeStep_ms,readIndex,nextFreeIndex);
  DPSNNverboseEnd();
}

void axonalSpikeSchedulerClass::extractAxSpike(const uint32_t source_h,
					       backwardAxonalSpikesClass *returnAxSpike,
					       const uint32_t delay)
{
  uint32_t i;
  uint32_t axSpikeCount;
  uint32_t index;

  DPSNNverboseStart(false,1,0);    
  printf("extractAxSpike: ENTER source_h=%03d delay=%d \n", source_h,delay);
  DPSNNverboseEnd();


  if(resetCalled != true) {
    printf("ERROR miss axonalSpikeScheduler.reset() call \n");
    fflush(stdout);exit(0);        
  };

  if(readIndex == -1) {
    printf("ERROR axonalSpikeScheduler empty \n");
    fflush(stdout);exit(0);        
  };
  
  index = (readIndex-delay+circBufferLen) % circBufferLen;
  //EPA-PSP 2015-03-09
  //returnAxSpike = &circBuffer[index];
  axSpikeCount = circBuffer[index].count;
  returnAxSpike[source_h].count = circBuffer[index].count;
  returnAxSpike[source_h].expectedCount = circBuffer[index].expectedCount;
  for(i=0;i<axSpikeCount;i++)
    returnAxSpike[source_h].list[i] = circBuffer[index].list[i];

  DPSNNverboseStart(false,1,0);    
  printf("extractAxSpike: EXIT source_h=%03d delay=%d \n", source_h,delay);
  DPSNNverboseEnd();

}



