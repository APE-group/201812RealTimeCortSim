// DPSNN_spike.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#ifndef DPSNN_spikeIncluded
#define DPSNN_spikeIncluded

#include "DPSNN_dataStructDims.h"

class axonalSpikeDataOnlyClass {
 public:
  double originalEmissionTime;//original emission time

  // BEWARE: the .pre_glob_n field of axSpike_forwardBuffer[.]  was
  // hijacked (by a union) with the total number of spikes expected
  union {
    uint32_t pre_glob_n;//global identity of presynaptic neuron
    uint32_t spikes_in_packet;// signature of last (fake) spike in
			      // packet bearing how many spikes are
			      // expected in the transfer
  };
  // ATTENTION: the behavoiur of MPI seems uncorrect if the uint32 is before the double
  // Just a padding problem??
};

class axonalSpikeClass : public axonalSpikeDataOnlyClass {
 public:
  void axonalSpikeToSynapticSpikes();
};

class synapticSpikeClass {
 public:
#ifdef LIFCAneuron
  double synDelayPlusOriginalEmissionTime;
#else
  int32_t synDelayPlusOriginalEmissionTime; 
#endif
  //index of target synapse in backwardSynList   
  uint32_t synIndex;
};

class synapticSpikeSchedulerClass {
public:
  void reset(  
     const uint32_t owner_glob_n_initValue,
     const uint32_t owner_loc_n_initValue,
     const uint32_t owner_loc_h_initValue);

  void insertSynSpike(const uint32_t target_owner_glob_n,
    const synapticSpikeClass tempSynSpike);
  bool extractSynSpike(synapticSpikeClass *returnSynSpike);
  bool mustBeEmpty(const int32_t time);
  int32_t getReadTime();
  void setReadTime(const int32_t time);
private:
  uint32_t circBuffer[ DSD__maxSimultaneousSpikesOnSameTarget ];
  uint32_t fillCountByCircTime;
  int32_t lastReadTime;
  int32_t readCircTime;
  uint32_t owner_glob_n, owner_loc_n, owner_loc_h;
  bool resetCalled;
};

class axonalSpikesClass {
 public:
  //EPA-PSP 2015-03-09
  //axonalSpikeDataOnlyClass list[DSD__maxAxonalSpike]; 
  axonalSpikeDataOnlyClass *list;
  //axonalSpikeClass 
  uint32_t count; // number of outbound spikes?
  uint32_t expectedCount; // number of inbound spikes?
};

class forwardAxonalSpikesClass : public axonalSpikesClass{ 
};
class backwardAxonalSpikesClass : public axonalSpikesClass{ 
};

class axonalSpikeSchedulerClass {
public:
  void reset();

  void insertAxSpike(const uint32_t source_h,
		     const backwardAxonalSpikesClass *tempAxSpike,
		     const uint32_t thisSinTimeStep_ms);
  void extractAxSpike(const uint32_t source_h,
		      backwardAxonalSpikesClass *returnAxSpike,
		      const uint32_t delay);

  backwardAxonalSpikesClass 
    circBuffer[ DSD__maxD ];
private:
  int32_t readIndex;
  int32_t nextFreeIndex;
  int32_t circBufferLen;
  bool resetCalled;
};



#endif //inclusion guard end
