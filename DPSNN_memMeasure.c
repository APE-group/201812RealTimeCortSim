// DPSNN_memMeasure.c
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
#include "DPSNN_memMeasure.h"
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_random.h"
#include "DPSNN_neuron.h"
#include "DPSNN_messagePassing.h"
#include "DPSNN_stopWatch.h"
#include "DPSNN_spike.h"
#include "DPSNN_stat.h"
#include "DPSNN_chrono.h"
#include "DPSNN_localNet.h"

#ifdef DALonlyEnvironmentSelected
#endif

void memMeasureClass::clear() {
  totalInitMM=0;
  totalSimMM=0;
}; 

void memMeasureClass::measure(const struct DPSNN_parameters lnp_par) {
  memoryPoolMM = sizeof(synapseClass) * lnp_par.locN*lnp_par.M*2;//DSD__maxBackwardLocSyn;  
  localNetMM=sizeof(localNetClass);
  forwardSynMM = sizeof(forwardSynClass);
  messagePassingMM = sizeof(messagePassingClass);
  synTargetHostDistributionMM = sizeof(synapticDistributionClass);
  synapseMM = sizeof(synapseClass);
  synSourceHostDistributionMM = sizeof(synapticDistributionClass);
  synDelaySourceHostDistributionMM = sizeof(synapticDelayDistributionClass);
  synapticSearchMM = sizeof(synapticSearchClass);
  backwardSynapticSearchMM = sizeof(synapticSearchClass) * lnp_par.globH;
  neuronMM=sizeof(neuronClass);
  synapticSpikeSchedulerMM = sizeof(synapticSpikeSchedulerClass);
  spikingStatusInThisTimeStepMM = sizeof(uint32_t) * lnp_par.locN;
  spikingNeuronIdsInThisTimeStepMM = sizeof(uint32_t) * lnp_par.locN;
  forwardAxonalSpikesMM = sizeof(forwardAxonalSpikesClass) * lnp_par.globH;
  backwardAxonalSpikesMM = sizeof(backwardAxonalSpikesClass) * lnp_par.globH;
  axonalSpikeSchedulerMM = sizeof(axonalSpikeSchedulerClass) * lnp_par.globH;
  delayedBackwardAxonalSpikesMM = sizeof(backwardAxonalSpikesClass) * lnp_par.globH;

  forwardSynListMM = forwardSynMM * lnp_par.locN*lnp_par.M;
  radixSortForwardSynMM = forwardSynMM * lnp_par.locN*lnp_par.M;
  intermSynListMM = forwardSynMM * lnp_par.locN*lnp_par.M*2;
  backwardSynListMM = synapseMM * lnp_par.locN*lnp_par.M*2;

  totalInitMM+=localNetMM;
  totalInitMM+= memoryPoolMM;
  totalInitMM += intermSynListMM;
  totalInitMM += messagePassingMM;
  totalSimMM += localNetMM;
  totalSimMM += memoryPoolMM;
  totalSimMM += messagePassingMM;

};

void memMeasureClass::report(const struct DPSNN_parameters lnp_par) {
  uint64_t totalSyn;
  totalSyn = (uint64_t)DSD__maxLocN*(uint64_t)DSD__maxGlobH*(uint64_t)DSD__maxM;
  printf("memMeasure-001 maxSynapses=maxLocN*maxGlobH*maxM = %llu\n",
	 (long long unsigned int)totalSyn);fflush(stdout);
  totalSyn = (uint64_t)lnp_par.locN*(uint64_t)lnp_par.globH*(uint64_t)lnp_par.M;
  printf("memMeasure-002 actualSynapses=locN*globH*M = %llu\n",
	 (long long unsigned int)totalSyn);fflush(stdout);
  printf("memMeasure-005 totalInitMM=%llu\n",(long long unsigned int)totalInitMM);fflush(stdout);
  printf("memMeasure-010 totalSimMM=%llu\n",(long long unsigned int)totalSimMM);fflush(stdout);
  printf("memMeasure-015 localNet=%llu\n",(long long unsigned int)localNetMM);fflush(stdout);
  printf("memMeasure-020 memoyPoolMM=%llu\n",(long long unsigned int)memoryPoolMM);fflush(stdout);
  printf("memMeasure-025 messagePassingMM = %llu\n",
	 (long long unsigned int)messagePassingMM);fflush(stdout);
  printf("memMeasure-030 forwardSynList=%llu\n",(long long unsigned int)forwardSynListMM);fflush(stdout); 
  printf("memMeasure-035 sort of forwardSyn=%llu\n",(long long unsigned int)radixSortForwardSynMM);fflush(stdout);
  printf("memMeasure-040 intermSynList=%llu\n",(long long unsigned int)intermSynListMM);fflush(stdout);
  printf("memMeasure-045 backwardSynList=%llu\n",(long long unsigned int)backwardSynListMM);fflush(stdout);
  printf("memMeasure-050 each forwardSynapse=%llu\n",(long long unsigned int)forwardSynMM);fflush(stdout);
  printf("memMeasure-055 each synapse=%llu\n",(long long unsigned int)synapseMM);fflush(stdout);
  printf("memMeasure-060 synTargetHostDistribution=%llu\n",
	 (long long unsigned int)synTargetHostDistributionMM);fflush(stdout);
  printf("memMeasure-065 synSourceHostDistribution=%llu\n",
	 (long long unsigned int)synSourceHostDistributionMM);fflush(stdout);
    printf("memMeasure-070 synDelaySourceHostDistribution=%llu\n",
	 (long long unsigned int)synDelaySourceHostDistributionMM);fflush(stdout);
  printf("memMeasure-075 backwardSynapticSearch=%llu\n",
	 (long long unsigned int)backwardSynapticSearchMM);fflush(stdout);
  printf("memMeasure-080 each neuron=%llu\n",
	 (long long unsigned int)neuronMM);fflush(stdout); fflush(stdout); 
  printf("memMeasure-090 n[locN=%d]=%llu\n",
	 lnp_par.locN,(long long unsigned int)(neuronMM * lnp_par.locN));fflush(stdout);
  printf("memMeasure-095 spikingStatusInThisTimeStepMM=%llu\n",
	 (long long unsigned int)spikingStatusInThisTimeStepMM);fflush(stdout);
  printf("memMeasure-100 spikingNeuronIdsInThisTimeStepMM=%llu\n",
	 (long long unsigned int)spikingNeuronIdsInThisTimeStepMM);fflush(stdout);
  printf("memMeasure-105 forwardAxonalSpikes[H=%d]=%llu\n",
	 lnp_par.globH,(long long unsigned int)forwardAxonalSpikesMM);fflush(stdout);
  printf("memMeasure-110 backwardAxonalSpikes[H=%d]=%llu\n",
	 lnp_par.globH,(long long unsigned int)backwardAxonalSpikesMM);fflush(stdout);
  printf("memMeasure-115 axonalSpikeScheduler=%llu\n",
	 (long long unsigned int)axonalSpikeSchedulerMM);fflush(stdout);
  printf("memMeasure-120 delayedBackwardAxonalSpikes[H=%d]=%llu\n",
	 lnp_par.globH,(long long unsigned int)delayedBackwardAxonalSpikesMM);fflush(stdout);
  printf("memMeasure-125 DSD__maxForwardLocSyn=%llu\n",
	 (long long unsigned int)DSD__maxForwardLocSyn);fflush(stdout);
  printf("memMeasure-130 DSD__maxBackwardLocSyn=%llu\n",
	 (long long unsigned int)DSD__maxBackwardLocSyn);fflush(stdout);
  printf("memMeasure-135 DSD__maxLocN=%d\n",DSD__maxLocN);fflush(stdout);
  printf("memMeasure-140 actual locN=%d\n",lnp_par.locN);fflush(stdout);
  printf("memMeasure-145 DSD__maxM=%d\n",DSD__maxM);fflush(stdout);
  printf("memMeasure-150 actual M=%d\n",lnp_par.M);fflush(stdout);
  printf("memMeasure-155 DSD__maxM_pre=%d\n",DSD__maxM_pre);fflush(stdout);
  printf("memMeasure-160 DSD__maxGlobH=%d\n",DSD__maxGlobH);fflush(stdout);
  printf("memMeasure-165 actual globH=%d\n",lnp_par.globH);fflush(stdout);
  printf("memMeasure-170 DSD__maxGlobN=%d\n",DSD__maxGlobN);fflush(stdout);
  printf("memMeasure-175 actual globN=%d\n",lnp_par.globN);fflush(stdout);
  printf("memMeasure-180 DSD__maxAxonalSpike=%llu\n",
	 (long long unsigned int)DSD__maxAxonalSpike);fflush(stdout);
  printf("memMeasure-185 DSD__maxMemoryPool=%llu\n",
	 (long long unsigned int)DSD__maxSizeOfSyn * (long long unsigned int)DSD__maxBackwardLocSyn);fflush(stdout);
	 //	 (long long unsigned int)DSD__maxMemoryPool);fflush(stdout);
};  

