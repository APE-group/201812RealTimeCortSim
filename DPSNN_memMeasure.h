// DPSNN_memMeasure.h
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

#ifndef DPSNN_memMeasureIncluded
#define DPSNN_memMeasureIncluded
 
#ifdef DALonlyEnvironmentSelected
#endif

class memMeasureClass { 
  private:
    uint64_t totalInitMM;
    uint64_t totalSimMM;
    uint64_t forwardSynMM;
    uint64_t localNetMM;
    uint64_t messagePassingMM;
    uint64_t synapseMM;
    uint64_t forwardSynListMM;
    uint64_t radixSortForwardSynMM;
    uint64_t intermSynListMM;
    uint64_t synTargetHostDistributionMM;
    uint64_t backwardSynListMM;
    uint64_t synSourceHostDistributionMM;
    uint64_t synDelaySourceHostDistributionMM;
    uint64_t synapticSearchMM;
    uint64_t backwardSynapticSearchMM;
    uint64_t neuronMM;
    uint64_t spikingStatusInThisTimeStepMM;
    uint64_t spikingNeuronIdsInThisTimeStepMM;
    uint64_t synapticSpikeSchedulerMM;
    uint64_t forwardAxonalSpikesMM;
    uint64_t backwardAxonalSpikesMM;
    uint64_t axonalSpikeSchedulerMM;
    uint64_t delayedBackwardAxonalSpikesMM;
    uint64_t memoryPoolMM;
  public:
    void clear();
    void measure(const struct DPSNN_parameters lnp_par);
    void report(const struct DPSNN_parameters lnp_par);
};

#endif //multiple inclusion guard
