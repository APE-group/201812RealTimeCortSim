// DPSNN_localNet.h 
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#ifndef DPSNN_localNetIncluded
#define DPSNN_localNetIncluded

#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_random.h"
#include "DPSNN_neuron.h"
#include "DPSNN_messagePassing.h"
#include "DPSNN_stopWatch.h"
#include "DPSNN_spike.h"
#include "DPSNN_stat.h"
#include "DPSNN_chrono.h"
#include "DPSNN_memMeasure.h"
#ifndef LIFCAneuron
#include "DPSNN_connectome.h"
#else
#include "DPSNN_LIFCAconnectome.h"
#include "erflib.h"
#include "randdev.h"
#endif

/*number of outgoing synapses per neuron*/
/*max number of incoming synapses per neuron ORROR 3*M already there in Inzhikevich code */
const char forwardString[]="FORWARD";
const char backwardString[]="BACKWARD";

class localNetClass {

private:

  //the parameters of the simulation
  struct DPSNN_parameters lnp_par;
  //the distributed message passing system
  messagePassingClass *pMessagePassing;
  //the simulation stopWatch
  stopWatchClass *pStopWatch;
  //the measurement and reporting system
  statClass *pStat;

  //lnp_par.loc_h; local host on which the local cluster of neurons resides
  //lnp_par.globH; total number of hosts in the system
  //lnp_par.locN, lnp_par.log2_locN  actual number of neurons in local network
  //lnp_par.globN; actual number of neurons in the system
  //lnp_par.M number of output synapses per neuron 
  //lnp_par.D maximum possible axo-synaptic delay in ms
  //lnp_par.locNe number of excitatory neurons in this cluster
  //lnp_par.locNi number of inhibitory neurons in this cluster

  char name[80];

  uint32_t thisSimTimeStep_ms;
  uint32_t startPartialChrono_ms;
  uint32_t stopPartialChrono_ms;
  
  uint8_t *memPoolA;
  uint8_t *memPoolB;
 
  synapseClass *localSynList;
  uint32_t localSynCount;

  synapseClass *forwardSynList;
  synapseClass *q1;//used in radixSort
  synapseClass *intermSynList;

  //old
  //synapseClass forwardSynList[DSD__maxBackwardLocSyn];

  synapticDistributionClass  synTargetHostDistribution;
  uint32_t forwardSynCount;
  uint32_t projectedSynCount; //local + forward

  //        Implementation of the
  //  backwardSynapses:backwardSynapsesClass
  //        box in the UML class diagram
  //synapseClass backwardSynList[DSD__maxBackwardLocSyn];
  synapseClass *backwardSynList;
  synapticDistributionClass  synSourceHostDistribution;
  synapticDelayDistributionClass synDelaySourceHostDistribution; 
  synapticSearchClass backwardSynapticSearch[DSD__maxGlobH];
  synOffsetInSynList backwardSynOffsetInSynList[DSD__maxGlobN];

  uint32_t backwardSynCount;

  // set of neurons composing this localNet
  // distribution of synaptic host targets and source hosts
  neuronClass n[DSD__maxLocN];
 
  uint32_t spikingStatusInThisTimeStep[DSD__maxLocN];
  uint32_t spikingNeuronIdsInThisTimeStep[DSD__maxAxonalSpike];
  uint32_t spikeCountInThisTimeStep;
  uint32_t compressedSpikingOffset;

  forwardAxonalSpikesClass forwardAxonalSpikes[DSD__maxGlobH];
  backwardAxonalSpikesClass backwardAxonalSpikes[DSD__maxGlobH];
  axonalSpikeSchedulerClass axonalSpikeScheduler[DSD__maxGlobH];
  backwardAxonalSpikesClass delayedBackwardAxonalSpikes[DSD__maxGlobH];
  void prepareAxonalSpikeBuffers();

  uint32_t N_firingsTotInFrame; //number of spikes in this frame (millisec)
  uint32_t N_thalamicInputs; //number of thalamic inputs in this frame
  //uint32_t N_firingsPerPop[lnp_par.globCFT * neuSubPopTotal]; //number of spikes for each population in this frame (millisec)

  // Firings in ChronoWindow and some statistic on it
  uint32_t N_firingsInChronoWindow;
  uint32_t minFiringsInChronoWindow;
  uint32_t maxFiringsInChronoWindow;
  double meanFiringsInChronoWindow;
  double varianceFiringsInChronoWindow;
  double sigmaFiringsInChronoWindow;
  double coeffOfVariationFiringsInChronoWindow;
  double M2FiringsInChronoWindow;
  uint32_t numLapsInChronoWindow;

  uint32_t reportCount;

  enum synListKind {forward, backward, unknown};

  float inputDistrib[DSD__maxLocN];
  
public:
  //initial construction of synapses and neurons
  void init(const struct DPSNN_parameters lnp_par_initValue, 
	    messagePassingClass * pMessagePassing_initValue,
	    stopWatchClass * pStopWatch_initValue,
	    statClass * pStat);

  void sendReceiveSynListWithOtherLocalNets();
  void inflateIntermToBackwardSynList(
    const synapseClass *intermSynList, 
    synapseClass *backwardSynList, const uint32_t backwardSynCount);

  chronoClass chronoSendReceiveSynList;

  void reportLocalNetAfterInit();

  void printStatFiringsInChronoWindow(uint32_t loc_h);

  //simulation steps
  void newFrameReset_s();
  void completeTimeStep_ms();
  void clearAllChronometers();
  void printAllInitChronoResults();
  void printAllSimulChronoResults();
  void printAllStatChronoResults();
  			
private:
  void prepareForwardSynapses();  

  memMeasureClass memMeasure;
  chronoClass chronoTimeStep;
  chronoClass chronoLPTAndAfterSpikeCalc;
  chronoClass chronoAfterSendRecSpikes;
  chronoClass chronoDynamicOfNeu;
  chronoClass chronoThalamicInput;
  chronoClass chronoAddSynCurrAndLTD;
  chronoClass chronoNeuralDynamic;
  chronoClass chronoRastergram;
  chronoClass chronoPlasticity;
  chronoClass chronoComputeDetail;
  chronoClass chronoComputeDetail2;
  chronoClass chronoBarrier1;
  chronoClass chronoBarrier2;
  chronoClass chronoBarrier3;
  chronoClass chronoPartialTimeStep;
  chronoClass chronoStatFunctions;

  //private functions to implement the simulation steps
  void longTermPlasticity();

  void addThalamicInput(const int time_ms, uint32_t *pThalInputCounter);

  void clearSpikingRecordForThisTimeStep();

  void compressNeuralSpikesInThisTimeStep();
  chronoClass chronoCompressNeuralSpikes;

  void logCompressedSpikesOfThisTimeStep();

  void neuralSpikesToAxonalSpikes();
  chronoClass chronoNeuralSpikesToAxonalSpikes;

  void sendReceiveAxonalSpikes();
  chronoClass chronoSendReceiveAxonalSpikes;
  chronoClass chronoExchangeAxonalSpikesDim;
  chronoClass chronoExchangeAxonalSpikes;

  void backwardAxonalSpikesToSynapticSpikeSchedulers();
  chronoClass chronoBackwardAxonalSpikersToSchedulers;

  void backwardAxonalSpikesToAxonalSpikeSchedulers();
  void axonalSpikeSchedulersToSynapticSpikeSchedulers();

  void sortTargetHostInForwardSynList(
	    synapseClass *synList,
	    synapseClass *ausiliarySynList,
	    uint32_t synCount);
  chronoClass chronoSortTargetHostInForwardSynList;

  void radixSortBinForwardSynList(synapseClass *synList, 
       synapseClass *ausiliarySynList, uint32_t synCount);
  chronoClass chronoRadixSortBinForwardSynList;

  //debug, check  reporting
  void printBackwardAxonalSpikes();

  void forwardSynListReport( 
     synapseClass *synList, 
     uint32_t synCount, uint32_t repNum);

  void synapticListReport( 
     synapseClass *synList, 
     uint32_t synCount, uint32_t repNum);
	
  void check_lnp_par_initValues();	
  void checkTotalSynCountInSynapticDistribution(
	    synapticDistributionClass *pSynapticDistribution, 
	    uint32_t expectedSynCount);
  void checkSortedSynListAndCreateSynapticDistribution(
	    synapseClass synList[], 
            synapticDistributionClass *pSynHostDistribution, 
            uint32_t synCount);
	
  void checkForwardSynListInitValues( 
            synapseClass *synList, 
            uint32_t synCount);

  void checkSynListInitValues(
            synapseClass *synList, 
            uint32_t synCount);

#if defined(makeActiveLTD) || defined (makeActiveLTP)
  void checkBackwardSynCountAfterInit();
#endif

  void timeStepSimulation();

  // Methods to manage ChronoWindow and Firing in ChronoWindow
  void initChronoWindow();
  bool isChronoWindow(uint32_t thisSimTimeStep_ms);
  void clearFiringsInChronoWindow();
  void firingsInChronoWindow(uint32_t spikeCountInThisTimeStep);

  //section added to implement Cortical Modules
  //numerosity of different kind of neurons
private: 
  simpleCM_connectomeClass simpleCM_connectome;
private:
  //methods to describe Cortical Modules
  void simpleCM_prepareForwardSynapses();
  void simpleCM_invokeIndivSynGenerators();

uint32_t intraModule_targetNeu_globId(
  uint32_t sourceModular_nId,
  uint32_t mId, uint32_t cxId, uint32_t cyId,
  neuralKindEnum sourceNeuKind,
  synGenEnum synGenKind);

uint32_t probableInterModule_targetNeu_globId(
  uint32_t sourceNeuIdInCM,
  uint32_t mId, uint32_t source_cxId, uint32_t source_cyId,
  uint32_t target_cxId, uint32_t target_cyId,
  neuralKindEnum sourceNeuKind,
  synGenEnum synGenKind);

uint32_t intraMod_nonRand_targetNeu_globId(
  uint32_t sourceNeu_idInCM,
  uint32_t mId, uint32_t cxId, uint32_t cyId,
  neuralKindEnum sourceNeuKind);

uint32_t probableInterMod_nonRand_targetNeu_globId(
  uint32_t sourceNeu_idInCM,
  uint32_t mId, uint32_t source_cxId, uint32_t source_cyId,
  uint32_t target_cxId, uint32_t target_cyId,
  neuralKindEnum sourceNeuKind);  

#ifdef LIFCAneuron
  erflibClass *pErfLib;
  randdevClass localNetRandDev;
  neuSubPopParamStruct neuSubPopParam[DSD__subPopTot];
  double spikingNeuronEmissionTime[DSD__maxLocN];
  uint32_t N_firingsPerPop[DSD__maxGlobCFT*DSD__subPopTot];
  uint32_t neuSubPopCount[DSD__subPopTot];

  double tableLUT_synMatrix[DSD__subPopTot][DSD__subPopTot][ANALOG_DEPTH];
  double tableLUT_synBathMatrix[DSD__subPopTot][ANALOG_DEPTH];

  double tableLUT_synBath_Exc[ANALOG_DEPTH];
  double tableLUT_synBath_Inh[ANALOG_DEPTH];

  uint32_t getRandomDelay_EXP(neuSubPopParamStruct *neuSubPopParam,neuSubPopEnum sourceSubPop,neuSubPopEnum targetSubPop);
  uint32_t getRandomDelay_UNI(neuSubPopParamStruct *neuSubPopParam,neuSubPopEnum sourceSubPop,neuSubPopEnum targetSubPop);
  void initLUT();
  void writeIniFiles();
#endif
  uint32_t forwardAxonalSpikesCount[DSD__maxGlobH];
  
};

#endif //multiple inclusion guard



