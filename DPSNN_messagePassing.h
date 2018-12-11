// DPSNN_messagePassing.h
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
#include <stddef.h>

#ifndef DPSNN_messagesIncluded
#define DPSNN_messagesIncluded

#include "DPSNN_random.h"

#include "DPSNN_environmentSelection.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"
#include "DPSNN_spike.h"
#include "DPSNN_messageMeasure.h"
#ifdef APENET_RDMA
#include "DPSNN_apenet.h"
#endif

//if included the next line generates the following errors in MPI environment
//#include "localNetProcess.h"
//error: field ‘localNeuralNetObject’ has incomplete type
//error: ‘messagePassingClass’ does not name a type

class messagePassingClass {
  DPSNN_parameters lnp_par;
  DALProcess *p;
  int axSpike_forwardPrep[DSD__maxGlobH];
  int axSpike_backwardPrep[DSD__maxGlobH];

  int axSpike_forwardCount[DSD__maxGlobH];
  int axSpike_forwardOffset[DSD__maxGlobH];
  axonalSpikeDataOnlyClass axSpike_forwardBuffer[DSD__maxAxonalSpike * DSD__maxGlobH];

  int axSpike_backwardCount[DSD__maxGlobH];
  int axSpike_backwardOffset[DSD__maxGlobH];
  axonalSpikeDataOnlyClass axSpike_backwardBuffer[DSD__maxAxonalSpike * DSD__maxGlobH];
 
public:
#ifdef APENET_RDMA
  APE_RDMA_desc_t* ape_descriptor;

  char * apeForwardBuffer[DSD__maxGlobH];
  volatile char * apeBackwardBuffer[DSD__maxGlobH];
  char * apeRemoteBackwardBuffer[DSD__maxGlobH];
#define APE_BUFFER_SIZE AM_MAX_MSG_SIZE*32

#endif

  /* here starts a part needed for fixed packet length: */
  // define an MPI_Type for the spike to send a vector of spikes

  // BEWARE: remember that 'pktLength' cannot be lesser
  // than 2, so the packet bears at least 1 spike!
  const int blocklengths[2] = {1, 1};
  const MPI_Aint offsets[2] = {
    offsetof(axonalSpikeDataOnlyClass, originalEmissionTime),
    offsetof(axonalSpikeDataOnlyClass, spikes_in_packet)
    };
  const MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT}; // field types
						       // in axonalSpikeDataOnlyClass
  MPI_Datatype MPI_spike;
  MPI_Aint lb_MPI_spike = 0;

  MPI_Aint MPI_fwdOffsetBaseAddress;
  MPI_Aint MPI_bwdOffsetBaseAddress;

  void init_OffsetBaseAddr(axonalSpikeDataOnlyClass*, axonalSpikeDataOnlyClass*);

  void set_fwdOffsets(axonalSpikeDataOnlyClass* , int);
  void set_bwdOffsets(axonalSpikeDataOnlyClass* , int);
  /* here ends a part needed for fixed packet length */

  void init(const DPSNN_parameters lnp_par_initValue,
	    DALProcess *p_initValue, const int thisTime_ms);

  void barrier();
  void sendForwardSynListDimToRemoteHosts(
	synapticDistributionClass *synTargetHostDistribution, 
	synapticDistributionClass *synSourceHostDistribution,
	const uint32_t thisTime_ms);

  void sendForwardSynListToRemoteHosts(
	    synapseClass *forwardSynList, 
	    synapticDistributionClass *synTargetHostDistribution, 
	    synapseClass *intermSynList, 
	    synapticDistributionClass *synSourceHostDistribution,
	    const uint32_t thisTime_ms);

  void exchangeAxonalSpikesDim(
    const forwardAxonalSpikesClass *pointForwardAxonalSpikes,
    const synapticDistributionClass *synTargetHostDistribution, 
    backwardAxonalSpikesClass *pointBackwardAxonalSpikes,
    const synapticDistributionClass *synSourceHostDistribution,
    const uint32_t thisTime_ms);

  void exchangeAxonalSpikes(
    forwardAxonalSpikesClass *pointForwardAxonalSpikes,
    backwardAxonalSpikesClass *pointBackwardAxonalSpikes,
    const uint32_t thisTime_ms);

  msgMeasureClass spikeDimSize;
  msgMeasureClass spikePayloadSize;

};

#endif
