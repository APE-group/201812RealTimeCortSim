// DPSNN_localNet_chrono.cpp 
// DPSNN-STDP project
// DPSNN_*.*  Distribution/Parallelization 
// of Polychronous Spiking Neural Networks 
// with synaptic Spiking Time Dependent Plasticity
// Pier Stanislao Paolucci 
// (Roma, Italy, project start date 2011)

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "DPSNN_debug.h"
#include "DPSNN_chrono.h"
#include "DPSNN_localNet.h"
#include "DPSNN_messagePassing.h"

void localNetClass::clearAllChronometers(){
  uint32_t totalSimTime;
  totalSimTime = atoi(getenv("env_totalSimTime_ms"));

  chronoSendReceiveSynList.clearChrono();
  chronoSortTargetHostInForwardSynList.clearChrono();

  chronoSendReceiveAxonalSpikes.clearChrono();
  chronoExchangeAxonalSpikesDim.clearChrono();
  chronoExchangeAxonalSpikes.clearChrono();

  chronoTimeStep.clearChrono();
  chronoLPTAndAfterSpikeCalc.clearChrono();
  chronoAfterSendRecSpikes.clearChrono();
  chronoDynamicOfNeu.clearChrono();
  chronoThalamicInput.clearChrono();
  chronoAddSynCurrAndLTD.clearChrono();
  chronoNeuralDynamic.clearChrono();
  chronoRastergram.clearChrono();
  chronoPlasticity.clearChrono();
  chronoComputeDetail.clearChrono();
  chronoComputeDetail2.clearChrono();
  chronoBarrier1.clearChrono();
  chronoBarrier2.clearChrono();
  chronoBarrier3.clearChrono();
  chronoPartialTimeStep.clearChrono();
  chronoStatFunctions.clearChrono();

  chronoCompressNeuralSpikes.clearChrono();
  chronoNeuralSpikesToAxonalSpikes.clearChrono();
  chronoBackwardAxonalSpikersToSchedulers.clearChrono();

}

void localNetClass::printAllInitChronoResults(){
  DPSNNverboseStart(true,1,0); 
    if(lnp_par.loc_h <= 1 ||  lnp_par.loc_h>=(lnp_par.globH-2)) {
      printf("CHRONO:sortTargetHostInForwardSynList h=%d = %f sec \n",
        lnp_par.loc_h, 
	chronoSortTargetHostInForwardSynList.getAccumulatedChrono());
      printf("CHRONO:sendReceiveSynList h=%d = %f sec \n",
        lnp_par.loc_h, chronoSendReceiveSynList.getAccumulatedChrono());
      fflush(stdout);};
  DPSNNverboseEnd();
};

void localNetClass::printAllSimulChronoResults(){
  double totalTime;
  uint32_t numSamples;
  double meanTime;
  uint32_t startPartialChrono_ms, stopPartialChrono_ms;

  totalTime = chronoLPTAndAfterSpikeCalc.getAccumulatedChrono();
  numSamples = chronoLPTAndAfterSpikeCalc.getNumLaps();
  meanTime =0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(true,1,0);
    printf("CHRONO-h=%d-A-BEFORE: each = %f sec sampled on %d iters, total time = %f sec\n",lnp_par.loc_h, meanTime, numSamples, totalTime);
    fflush(stdout);
  DPSNNverboseEnd();

  totalTime = chronoBarrier1.getAccumulatedChrono();
  numSamples = chronoBarrier1.getNumLaps();
  meanTime =0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(true,1,0);
    printf("CHRONO-h=%d-B-BARRIER: each = %f sec sampled on %d iters, total time = %f sec\n",lnp_par.loc_h, meanTime, numSamples, totalTime);
    fflush(stdout);
  DPSNNverboseEnd();

  totalTime =  chronoExchangeAxonalSpikesDim.getAccumulatedChrono();
  numSamples = chronoExchangeAxonalSpikesDim.getNumLaps();
  meanTime =0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(true,1,0);
    printf("CHRONO-h=%d-C-DIM: each = %f sec sampled on %d iters, total time = %f sec\n",lnp_par.loc_h, meanTime, numSamples, totalTime);
    fflush(stdout);
  DPSNNverboseEnd();

  totalTime =  chronoExchangeAxonalSpikes.getAccumulatedChrono();
  numSamples = chronoExchangeAxonalSpikes.getNumLaps();
  meanTime =0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(true,1,0);
    printf("CHRONO-h=%d-D-PAYLOAD: each = %f sec sampled on %d iters, total time = %f sec\n",lnp_par.loc_h, meanTime, numSamples, totalTime);
  DPSNNverboseEnd();

  totalTime =  chronoAfterSendRecSpikes.getAccumulatedChrono();
  numSamples = chronoAfterSendRecSpikes.getNumLaps();
  meanTime =0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(true,1,0);
  printf("CHRONO-h=%d-E-AFTER: each = %f sec sampled on %d iters, total time = %f sec\n",
	 lnp_par.loc_h, meanTime, numSamples, totalTime);
  fflush(stdout);
  DPSNNverboseEnd();

  totalTime = chronoTimeStep.getAccumulatedChrono();
  numSamples = chronoTimeStep.getNumLaps();
  meanTime = 0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(true,1,0);
  printf("CHRONO-h=%d-F-TOTAL TotalTimeStep: each timeStep = %f sec sampled on %d iters, total time = %f sec\n",lnp_par.loc_h, meanTime, numSamples, totalTime );
  fflush(stdout);
  DPSNNverboseEnd();

  totalTime =  chronoSendReceiveAxonalSpikes.getAccumulatedChrono();
  numSamples = chronoSendReceiveAxonalSpikes.getNumLaps();
  meanTime =0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(false,1,0);
    printf("CHRONO-h=%d-G-TOTAL SendReceiveSpikes: each sendRec = %f sec sampled on %d iters, total time = %f sec\n",lnp_par.loc_h, meanTime, numSamples, totalTime);
    fflush(stdout);
  DPSNNverboseEnd();
  
  totalTime = chronoPartialTimeStep.getAccumulatedChrono();
  numSamples = chronoPartialTimeStep.getNumLaps();
  startPartialChrono_ms = chronoPartialTimeStep.getstartPartialChrono_ms();
  stopPartialChrono_ms = chronoPartialTimeStep.getstopPartialChrono_ms();
  meanTime = 0.0; 
  if(numSamples>0) meanTime = (double)(totalTime/(double)numSamples);
  DPSNNverboseStart(false,1,0);
  printf("CHRONO-h=%d-H-***PARTIAL*** PartialTimeStep: simulation from ms %d to ms %d: partial time = %f sec\n",lnp_par.loc_h,startPartialChrono_ms,stopPartialChrono_ms,totalTime );
  fflush(stdout);
  DPSNNverboseEnd();
  
};

void localNetClass::printAllStatChronoResults(){
  
  chronoLPTAndAfterSpikeCalc.printChrono("BEFORE",lnp_par.loc_h,1);
  chronoBarrier1.printChrono("BARRIER",lnp_par.loc_h,0);
  chronoExchangeAxonalSpikesDim.printChrono("DIM",lnp_par.loc_h,0);
  chronoExchangeAxonalSpikes.printChrono("PAYLOAD",lnp_par.loc_h,0);
  chronoAfterSendRecSpikes.printChrono("AFTER",lnp_par.loc_h,0);
  chronoTimeStep.printChrono("TOTAL",lnp_par.loc_h,0);

  //chronoBarrier2.printChrono("AfterSendRec_Barrier",lnp_par.loc_h,0);
  chronoThalamicInput.printChrono("ThalamicInput",lnp_par.loc_h,0);
  chronoAddSynCurrAndLTD.printChrono("AddSynapticCurrents+LTD",lnp_par.loc_h,0);
  chronoNeuralDynamic.printChrono("NeuronDynamic",lnp_par.loc_h,0);
  chronoRastergram.printChrono("Rastergram",lnp_par.loc_h,0);
  chronoPlasticity.printChrono("Plasticity",lnp_par.loc_h,0);  
  chronoStatFunctions.printChrono("StatisticalFunctions",lnp_par.loc_h,0);
  chronoPartialTimeStep.printChrono("PartialTimeStep",lnp_par.loc_h,0);
  pMessagePassing->spikeDimSize.printMsgSize("MessageDimSize",lnp_par.loc_h,0);
  pMessagePassing->spikePayloadSize.printMsgSize("MessagePayloadSize",lnp_par.loc_h,0);

  DPSNNverboseStart(false,1,0);
    chronoDynamicOfNeu.printIstogramData("IstogramData1",lnp_par.loc_h,lnp_par.totalSimTime_ms);
    chronoThalamicInput.printIstogramData("IstogramData2",lnp_par.loc_h,lnp_par.totalSimTime_ms);
    chronoAddSynCurrAndLTD.printIstogramData("IstogramData3",lnp_par.loc_h,lnp_par.totalSimTime_ms);
    chronoNeuralDynamic.printIstogramData("IstogramData4",lnp_par.loc_h,lnp_par.totalSimTime_ms);
  DPSNNverboseEnd();
}
