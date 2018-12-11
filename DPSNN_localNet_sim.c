// DPSNN_localNet_sim.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include "DPSNN_debug.h"
#include "DPSNN_random.h"
#include "DPSNN_neuron.h"
#include "DPSNN_messagePassing.h"
#include "DPSNN_spike.h"
#include "DPSNN_localNet.h"

void localNetClass::completeTimeStep_ms() {

  stopWatchStatusEnum stopWatchStatus;
  uint32_t thisSim_sec;
  uint32_t i; 

  //the time is advanced by the localNetProcess fire and init 
  //here inside we read infos
  thisSimTimeStep_ms = pStopWatch->getSimTime_ms();

  DPSNNverboseStart(true,1,0);
  if(lnp_par.loc_h==0  && ((thisSimTimeStep_ms % 1000)==0)) {
    printf("START of completeTimeStep_ms at %d ms\n",
	   thisSimTimeStep_ms);
    fflush(stdout);
  }
  DPSNNverboseEnd();

 #ifdef MPIandDALenvironmentSelected
  {
    char hostNode[256];
    double gigabyte = 1024 * 1024 * 1024;
    struct sysinfo si;
  
    DPSNNverboseStart(false,1,0);
    if(lnp_par.loc_h==0 && ((thisSimTimeStep_ms % 100)==0)) {
      sysinfo (&si);
      gethostname(hostNode,sizeof(hostNode));
      printf ("During run at %d ms free memory on %s is %5.1f GB\n",
	      thisSimTimeStep_ms,hostNode?hostNode:"NULL",si.freeram / gigabyte);
      fflush(stdout);
    }
    DPSNNverboseEnd();
  }
  #endif
  
  chronoTimeStep.startChrono();
  chronoPartialTimeStep.startChrono();
  chronoLPTAndAfterSpikeCalc.startChrono();

  thisSim_sec = thisSimTimeStep_ms/1000;

  stopWatchStatus = pStopWatch->getStatus();
  
  if(stopWatchStatus == new_s) {
    N_firingsTotInFrame=0;//to avoid the location 0 on the firings[][] buffer      
    N_thalamicInputs=0;

    //reset of spike statistic over the frame
    pStat->openFiles_forThisSecond(thisSim_sec, lnp_par.moduloSec);

    //print synapses
    DPSNNverboseStart(false,thisSimTimeStep_ms,0);
      for (i=0;i<backwardSynCount;i++)  {
        pStat->syn.write(thisSim_sec, lnp_par.moduloSec, 
	  thisSimTimeStep_ms,i, 
	  backwardSynList[i].pre_glob_n,
	  backwardSynList[i].post_glob_n,
	  backwardSynList[i].delay,
	  (float)backwardSynList[i].weight,
	  (float)backwardSynList[i].timeDerivative);
      };
    DPSNNverboseEnd();
  }

  //reset of ms statistic
  clearSpikingRecordForThisTimeStep();

  //only for those neurons which fired:
  //- axonal spikes are added to the list of signals to be sent
  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
      printf("just before test of spiking at %d ms\n",
            thisSimTimeStep_ms); 
  DPSNNverboseEnd();

  //check of spiking (in previous time step, actually)
  //and consequent actions
  //this order (check before dyn) 
  //is momentarily kept for compatibility
  //with the sequential version 

  for (i=0;i<lnp_par.locN;i++) {
    if(n[i].didItFire) {
      n[i].didItFire = false;
      spikingStatusInThisTimeStep[i] = 1;
      spikeCountInThisTimeStep++;

      n[i].setLastEmittedSpikeTime_ms(thisSimTimeStep_ms);
#ifdef makeActiveLTP
      if((thisSimTimeStep_ms>=lnp_par.startPlasticity_ms)&&
	 (thisSimTimeStep_ms<=lnp_par.stopPlasticity_ms)) {
        n[i].causalSTDP_ms_ofBackwardSynapses(thisSimTimeStep_ms);
      }; 
#endif
      DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
        printf("YES on h=%d spiking glob_n=%d dat %d ms\n",
	       lnp_par.loc_h, n[i].getGlob_n(), thisSimTimeStep_ms);
      DPSNNverboseEnd();

      //n[i].afterSpikeNeuralDynamic(thisSimTimeStep_ms);
    } else {
      DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
        printf("NON spiking loc neu %d at %d ms\n",
             i,thisSimTimeStep_ms);
      DPSNNverboseEnd();
      spikingStatusInThisTimeStep[i]=0;
    };
  };

  //a list of spikes if created
  compressNeuralSpikesInThisTimeStep();

  //neural spikes are added to the list 
  //of axonal spikes to be transmitted
  neuralSpikesToAxonalSpikes();

  //each cluster of neuron sends the axonal spikes 
  //generated in this time step
  chronoLPTAndAfterSpikeCalc.stopChronoPartial(thisSimTimeStep_ms);

#ifdef makeInsertBarrierBeforeSendRec
  chronoBarrier1.startChrono();
  pMessagePassing->barrier();
  chronoBarrier1.stopChronoPartial(thisSimTimeStep_ms);
#endif

  sendReceiveAxonalSpikes();

  /*
  printf("********       at %d ms    source_h=%d",thisSimTimeStep_ms,lnp_par.loc_h);
  for(uint32_t target_h=0;target_h<lnp_par.globH;target_h++) 
    printf("    target_h=%d   count=%d",target_h,forwardAxonalSpikes[target_h].count);
  printf("  \n");
  */
  for(uint32_t target_h=0;target_h<lnp_par.globH;target_h++) 
    forwardAxonalSpikesCount[target_h] = forwardAxonalSpikes[target_h].count;
  pStat->messageTraffic.write(thisSim_sec,lnp_par.moduloSec,thisSimTimeStep_ms,
			      forwardAxonalSpikesCount);

#ifdef makeInsertBarrierAfterSendRec
  chronoBarrier2.startChrono();
  pMessagePassing->barrier();
  chronoBarrier2.stopChronoPartial(thisSimTimeStep_ms);
#endif
  
  //debugging
  //#define debugAxonSpikes
  #ifdef debugAxonSpikes
    printBackwardAxonalSpikes();
  #endif
  #undef debugAxonSpikes
    
  chronoAfterSendRecSpikes.startChrono();
  chronoAddSynCurrAndLTD.startChrono();

  for (i=0;i<lnp_par.locN;i++) {
    //first action for dynamic: clear input currents
    n[i].clearInputCurrent();
  };

  backwardAxonalSpikesToAxonalSpikeSchedulers();
  axonalSpikeSchedulersToSynapticSpikeSchedulers();

  chronoAddSynCurrAndLTD.stopChronoPartial(thisSimTimeStep_ms);

#ifndef LIFCAneuron
  chronoThalamicInput.startChrono();
  //computes thalamic input
  addThalamicInput(thisSimTimeStep_ms, &N_thalamicInputs);
  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  printf("N_thalamicInputs=%d\n",N_thalamicInputs);
  DPSNNverboseEnd();
  chronoThalamicInput.stopChronoPartial(thisSimTimeStep_ms);
#endif

  //---------here starts the dynamic------------

  chronoNeuralDynamic.startChrono();

  //regular update of each neuron state variables
  //using diffential equation
  for (i=0;i<lnp_par.locN;i++) {
    n[i].regularDynamicTimeStep(thisSimTimeStep_ms, inputDistrib, &N_thalamicInputs);
  };

  DPSNNverboseStart(false,thisSimTimeStep_ms,1999);
    printf("NEAR END: Just before spikes.write on loc_h=%d at ms=%d\n",
	 lnp_par.loc_h,thisSimTimeStep_ms);
  fflush(stdout);
  DPSNNverboseEnd();

  N_firingsTotInFrame += spikeCountInThisTimeStep;
  if(isChronoWindow(thisSimTimeStep_ms))
    firingsInChronoWindow(spikeCountInThisTimeStep);

  chronoNeuralDynamic.stopChronoPartial(thisSimTimeStep_ms);
  chronoAfterSendRecSpikes.stopChronoPartial(thisSimTimeStep_ms);

  chronoRastergram.startChrono();  

  pStat->spikes.write(thisSim_sec,  lnp_par.moduloSec, 
		      spikeCountInThisTimeStep, 
		      thisSimTimeStep_ms,
		      (const int *) spikingNeuronIdsInThisTimeStep,
		      (const double *) spikingNeuronEmissionTime);
  {
      uint32_t N_column;
      uint32_t N_globSubPop;
      uint32_t N_locSubPop;
      uint32_t N_firstLocSubPop;
      uint32_t spikingNeuInItsColumn;
      uint16_t subPop;

      N_firstLocSubPop = (lnp_par.first_locCFX + lnp_par.first_locCFY * lnp_par.globCFX) * lnp_par.subPopNumber;
      for(i=0;i<spikeCountInThisTimeStep;i++){
        N_column = (uint32_t)floor((double)(spikingNeuronIdsInThisTimeStep[i]/lnp_par.neuronsPerCM));
        spikingNeuInItsColumn = spikingNeuronIdsInThisTimeStep[i] - N_column * lnp_par.neuronsPerCM;
        subPop = 0;
        while(subPop < lnp_par.subPopNumber){
 	  if(spikingNeuInItsColumn < neuSubPopParam[subPop].count) break;
	  spikingNeuInItsColumn-=neuSubPopParam[subPop].count;
	  subPop++;
        }

        N_globSubPop = N_column * lnp_par.subPopNumber + subPop;
        N_locSubPop = N_globSubPop - N_firstLocSubPop;
        N_firingsPerPop[N_locSubPop]++;
      }

      pStat->spikingRatesPerPop.write(thisSim_sec,thisSimTimeStep_ms,lnp_par.moduloSec,
	       lnp_par.ratesSampling,N_firingsPerPop,neuSubPopCount,lnp_par.subPopNumber);
  }
 
   
  pStat->spikePerMs.write(thisSim_sec,  lnp_par.moduloSec, 
			  thisSimTimeStep_ms, spikeCountInThisTimeStep);  

  DPSNNverboseStart(false,thisSimTimeStep_ms,1999);
    printf("NEAR END: Just after spikes.write on loc_h=%03d at ms=%d\n",
	 lnp_par.loc_h,thisSimTimeStep_ms);
  fflush(stdout);
  DPSNNverboseEnd();
  
  chronoRastergram.stopChronoPartial(thisSimTimeStep_ms);

  if(stopWatchStatus == x999_ms_reached) {
 
    DPSNNverboseStart(true,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
    if(lnp_par.loc_h==0) {
	  printf("closing %d ms simul\n", thisSimTimeStep_ms);
	  fflush(stdout);}
    DPSNNverboseEnd();
    
#if defined(makeActiveLTD) || defined (makeActiveLTP)
    
    DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.startPlasticity_ms);
      chronoPlasticity.startChrono();
    DPSNNverboseEnd();
    
    if((thisSimTimeStep_ms>=lnp_par.startPlasticity_ms)&&
       (thisSimTimeStep_ms<=lnp_par.stopPlasticity_ms)) {
      longTermPlasticity();
    };

    DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.startPlasticity_ms);
      chronoPlasticity.stopChronoPartial(thisSimTimeStep_ms);  
    DPSNNverboseEnd();

    DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
    synapticListReport(backwardSynList, backwardSynCount, thisSim_sec);
    DPSNNverboseEnd();
#endif

    chronoStatFunctions.startChrono();

    //pStat->spikingRates.write(thisSim_sec,  
    //	      lnp_par.moduloSec, N_firingsTotInFrame, N_thalamicInputs);
 
    pStat->closeFiles_forThisSecond(thisSim_sec, lnp_par.moduloSec);

    chronoStatFunctions.stopChronoPartial(thisSimTimeStep_ms);   
  };
  
  chronoPartialTimeStep.stopChronoPartial(thisSimTimeStep_ms);
  chronoTimeStep.stopChronoPartial(thisSimTimeStep_ms); 
};

#if defined(makeActiveLTD) || defined (makeActiveLTP)
#ifdef LIFCAneuron

void localNetClass::longTermPlasticity() 
{
  //sequential code // modify only exc connections
  //for (i=0;i<cNe;i++) for (j=0;j<cM;j++) {
  //   s[i][j]+=0.01 + sd[i][j]; sd[i][j]*=0.9;			
  //   if (s[i][j] > cSm) s[i][j]= cSm; if (s[i][j] < 0) s[i][j]=0.0;}
  uint32_t i;
  uint32_t targetNeu_glob_n;
  uint32_t targetNeu_loc_n;
  uint32_t sourceNeu_glob_n;
  neuSubPopEnum sourcePop, targetPop;
  float maxSynValue;
  float reverseMaxSynValue;
  
  float weightF;
  float timeDerivativeF;
  float clampCorrectionF;
  /* instrumentedSynapse instrumentedSynapseDummy; */
  synapseClass tempSyn;
  
  /* FILE *fp_outputA, *fp_outputB; */
  /* char reportName[80]; */

  DPSNNverboseStart(true,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if(lnp_par.loc_h==0) {
    printf("longTermPlasticity Start: backwardSynCount=%d on h=%d at %d ms\n",backwardSynCount,lnp_par.loc_h,thisSimTimeStep_ms);
  }
  DPSNNverboseEnd();
  
  if(lnp_par.loc_h==0) {
    if(lnp_par.SongSTDP_param.tauDerivativeDecay_s < 1.0) {
      printf("ERROR in longTermPlasticity, tauDerivativeDecay_s < 1.0\n");
      fflush(stdout);exit(0);} 
  };
 
  for (i=0; i< backwardSynCount; i++)
  {
    sourcePop = (neuSubPopEnum)backwardSynList[i].preSynNeuralKind;
    sourceNeu_glob_n = backwardSynList[i].pre_glob_n;
    targetNeu_glob_n = backwardSynList[i].post_glob_n;
    targetNeu_loc_n = targetNeu_glob_n % lnp_par.locN;
    targetPop = n[targetNeu_loc_n].getSubPop_n();
    // 2.0 * mean value, to keep a range simmetrical bewteen 0.0 and 2.0 * meanValue
    maxSynValue = 2.0*neuSubPopParam[sourcePop].J[targetPop];

    //reverseMaxSynValue: if>0 the targetPop is excitatory 
    reverseMaxSynValue = 2.0*neuSubPopParam[targetPop].J[sourcePop];
    clampCorrectionF=0.0;

    if(sourceNeu_glob_n < lnp_par.locN) {
      tempSyn = backwardSynList[i];
      pStat->statSynPeriodicProbeA.write(thisSimTimeStep_ms/1000,lnp_par.moduloSec,tempSyn);
    }
    
    if (maxSynValue > 0.0 &&      //i.e. excitatory synapse
	reverseMaxSynValue > 0.0) //toward an excitatory pop
                                  //(dirty trick to be removed)
    {
      weightF = (float)(backwardSynList[i].weight * lnp_par.factorWeightType_2_Float);

      timeDerivativeF = backwardSynList[i].timeDerivative;
 
      //add the possible constant drift (it will be subtracted before the decay)
      timeDerivativeF += lnp_par.SongSTDP_param.derivativeDrift;      
      weightF = weightF + timeDerivativeF;
	    
      if(weightF > maxSynValue) {
	clampCorrectionF = (weightF - maxSynValue);
        timeDerivativeF = timeDerivativeF - clampCorrectionF;
	weightF = maxSynValue; } 
        
      if(weightF < 0.0) {
	clampCorrectionF = weightF;
        timeDerivativeF = timeDerivativeF - clampCorrectionF;
	weightF = 0.0 ;
      }
 
      backwardSynList[i].weight = (int16_t)(weightF * lnp_par.factorFloat_2_weightType);
      backwardSynList[i].timeDerivative = timeDerivativeF;
    };


    if(sourceNeu_glob_n < lnp_par.locN) {
      tempSyn = backwardSynList[i];
      pStat->statSynPeriodicProbeB.write(thisSimTimeStep_ms/1000,lnp_par.moduloSec,tempSyn);
    }
    
    //decay of the time derivative.
    //If tau is set to 1 sec, each second is independent
    //if tau is set to 10

    if(lnp_par.SongSTDP_param.tauDerivativeDecay_s == 1.0) {
      backwardSynList[i].timeDerivative = 0.0;
    }else{
      backwardSynList[i].timeDerivative =
	  (backwardSynList[i].timeDerivative -
	   lnp_par.SongSTDP_param.derivativeDrift +
	   clampCorrectionF) *
	  (1 - 1/lnp_par.SongSTDP_param.tauDerivativeDecay_s);
    };
  };

};
#else
void localNetClass::longTermPlasticity() 
{
  //sequential code // modify only exc connections
  //for (i=0;i<cNe;i++) for (j=0;j<cM;j++) {
  //   s[i][j]+=0.01 + sd[i][j]; sd[i][j]*=0.9;			
  //   if (s[i][j] > cSm) s[i][j]= cSm; if (s[i][j] < 0) s[i][j]=0.0;}
  uint32_t i;
  float synapseMaxValue;
  //DSD_cSm is the value to be set for M = 100
  //synapseMaxValue = DSD__cSm;
  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  printf("backwardSynCount=%d on h=%d at %d ms\n",backwardSynCount,lnp_par.loc_h,thisSimTimeStep_ms);
  DPSNNverboseEnd();

  synapseMaxValue = (float) makeInitExcitWeight * (10.0/6.14) * (100.0 / (float)lnp_par.M) ;

  for (i=0; i< backwardSynCount; i++)
  {
    if ((backwardSynList[i].preSynNeuralKind == excitatoryRS) ||
	(backwardSynList[i].preSynNeuralKind == excitatoryLbExc) ||
	(backwardSynList[i].preSynNeuralKind == excitatoryLaExc))
    {

      backwardSynList[i].weight = 
	(weightType)((float)backwardSynList[i].weight + 0.01 + (float)backwardSynList[i].timeDerivative);
      backwardSynList[i].timeDerivative = (weightType)((float)backwardSynList[i].timeDerivative * 0.9);
      if((float)backwardSynList[i].weight > synapseMaxValue) 
        backwardSynList[i].weight = (weightType)synapseMaxValue;
      if((float)backwardSynList[i].weight < 0.0) 
        backwardSynList[i].weight = (weightType)0.0;
    } else {
      if ((backwardSynList[i].preSynNeuralKind != inhibitoryFS) &&
	  (backwardSynList[i].preSynNeuralKind != inhibitoryLaInh)) {
	printf("ERROR preSynNeuralKind=%d wrong in longTermPlasticity for neu=%d \n",
	       backwardSynList[i].preSynNeuralKind,backwardSynList[i].pre_glob_n);
	printf("It should be inhibitoryFS(%d) or inhibitoryLaInh(%d) \n",
	       inhibitoryFS,inhibitoryLaInh);
	fflush(stdout);exit(0);
      };
    };
  };
};
#endif
#endif

void localNetClass::clearSpikingRecordForThisTimeStep()
{
    uint32_t i;
    spikeCountInThisTimeStep = 0;
    compressedSpikingOffset = 0;

    for (i=0;i<lnp_par.locN;i++) {
      spikingStatusInThisTimeStep[i] = 0;
      spikingNeuronIdsInThisTimeStep[i] = lnp_par.globN;//absurd value
      #ifdef LIFCAneuron
	spikingNeuronEmissionTime[i] = 0.0;
      #endif
    }
};

void localNetClass::addThalamicInput(
  const int thisSimTimeStep_ms, 
  uint32_t * pThalInputCounter)
{

  float thalamicInputCurrent;
  //default value
  //
  //WARNIG: Thal.InCurrent changed from 20 to 60 for DEBUGGING purp.
  //
  thalamicInputCurrent=20.0;

  if(lnp_par.fastDebugDyn == fastDebugDyn_1) {
      //false value to speed up simulation
	thalamicInputCurrent = 10000.0;
  }

  switch (lnp_par.thalamicInput) {
  case no_thalamicInput_0:
      break;
  case default_random_thalamicInput_1:
      {
	uint32_t i_loc_cfx, i_loc_cfy;
	uint32_t iThal, max_iThal;
	uint32_t i;
	uint32_t this_CF,loc_CF;
	uint32_t inputDistribLen;

	inputDistribLen = lnp_par.locCFT >= 1 ? lnp_par.locN : lnp_par.neuronsPerCM;
	max_iThal = (lnp_par.neuronsPerCM / 1024) * lnp_par.thalamicInputFreq;

	for (i=0; i<inputDistribLen;i++)
	  inputDistrib[i] = 0.0;

	for(i_loc_cfx  = lnp_par.first_locCFX; 
	    i_loc_cfx <= lnp_par.last_locCFX; i_loc_cfx++) {
	  for(i_loc_cfy  = lnp_par.first_locCFY; 
	      i_loc_cfy <= lnp_par.last_locCFY; i_loc_cfy++) {

	    this_CF = i_loc_cfx + i_loc_cfy * lnp_par.globCFX;
	    loc_CF = lnp_par.locCFT >= 1 ? this_CF % (uint32_t) lnp_par.locCFT : 0;
	    
	    srand((this_CF + 1) * thisSimTimeStep_ms);
	    for(iThal=0; iThal<max_iThal; iThal++)
	      inputDistrib[getRandom(lnp_par.RSPerCM)+loc_CF*lnp_par.neuronsPerCM] = 1.0;
	  }
	}
      };
      break;
  case linear_neuronsPerCM_vs_time_thalamicInput_3:
      { 
	uint32_t target_n_inCM, target_glob_n, target_loc_n;
	uint32_t i_loc_cfx, i_loc_cfy;

	for(i_loc_cfx  = lnp_par.first_locCFX; 
	    i_loc_cfx <= lnp_par.last_locCFX; i_loc_cfx++) {
	  for(i_loc_cfy  = lnp_par.first_locCFY; 
	      i_loc_cfy <= lnp_par.last_locCFY; i_loc_cfy++) {
	    //the target moves linearly with time (ms and sec) 
	    //and cm coordinates
	    
	    target_n_inCM =
	      ((thisSimTimeStep_ms%1000 + 
	      i_loc_cfx + i_loc_cfy * lnp_par.globCFX) )
	      % lnp_par.neuronsPerCM;
	    target_glob_n = target_n_inCM + 
	      i_loc_cfx * lnp_par.neuronsPerCM + 
	      i_loc_cfy * lnp_par.globCFX * lnp_par.neuronsPerCM;
	        	    
	    if((target_glob_n >= lnp_par.first_glob_n) &&
		(target_glob_n <= lnp_par.last_glob_n))
	    {
	      target_loc_n = target_glob_n % lnp_par.locN;
	      n[target_loc_n].addInputCurrent(thalamicInputCurrent);
              (*pThalInputCounter) ++;
              pStat->thalamicInput.write(
	      thisSimTimeStep_ms/1000, lnp_par.moduloSec, 
	      thisSimTimeStep_ms,  
              target_glob_n);
	    };
	  };
	};	
      };
      break;
  case random_perCM_perMs_thalamicInput_4:
      { 
	uint32_t target_n_inCM, target_glob_n, target_loc_n;
	uint32_t i_loc_cfx, i_loc_cfy;
	uint32_t iThal, max_iThal;
	uint32_t this_CF;
	uint32_t seedForThalInput;
	//--//long int seedForThalInput;
	//--//struct drand48_data *randBufferThalInput;
	//--//long int randomValue;


	max_iThal = (lnp_par.neuronsPerCM / 1024) * lnp_par.thalamicInputFreq;
	for(i_loc_cfx  = lnp_par.first_locCFX; 
	    i_loc_cfx <= lnp_par.last_locCFX; i_loc_cfx++) {
	  for(i_loc_cfy  = lnp_par.first_locCFY; 
	      i_loc_cfy <= lnp_par.last_locCFY; i_loc_cfy++) {

	      //srand((this_CF + 1) * thisSimTimeStep_ms);
	      this_CF = i_loc_cfx + i_loc_cfy * lnp_par.globCFX;
	      seedForThalInput = (this_CF + 1) * thisSimTimeStep_ms;  

	    for(iThal=0; iThal<max_iThal; iThal++) {
	      //target_n_inCM = getRandom(lnp_par.RSPerCM);
	      target_n_inCM = getRandom_r(&seedForThalInput,lnp_par.RSPerCM);  
	
	      //--//seedForThalInput = (this_CF + 1) * thisSimTimeStep_ms;
	      //--//srand48_r(seedForThalInput,randBufferThalInput);
	      //--//lrand48_r(randBufferThalInput,&randomValue);
	      //--//target_n_inCM = randomValue %(long int)lnp_par.RSPerCM;

	      target_glob_n = target_n_inCM + 
		  i_loc_cfx * lnp_par.neuronsPerCM + 
		  i_loc_cfy * lnp_par.globCFX * lnp_par.neuronsPerCM;
		
	      if((target_glob_n >= lnp_par.first_glob_n) &&
		   (target_glob_n <= lnp_par.last_glob_n))
	      {
		target_loc_n = target_glob_n % lnp_par.locN;
		n[target_loc_n].addInputCurrent(thalamicInputCurrent);
		(*pThalInputCounter) ++;
		pStat->thalamicInput.write(
		   thisSimTimeStep_ms/1000, lnp_par.moduloSec, 
		   thisSimTimeStep_ms,  
		   target_glob_n);
	      };
	    };
	  };
	};
      };
      break;

  default:
      {printf(
         "ERROR unrecognized thalamic input set by the environment\n");
      fflush(stdout);exit(0);
      };
  };
};

void localNetClass::compressNeuralSpikesInThisTimeStep() 
{
  //this could be parallelized with the following tricks
  //- divide the spikeId in segments
  //- divide the neurons in segments
  //- use lists for each segment
  //  at the end attach the lists of each segment
  uint32_t i;
  
  for (i=0; i < lnp_par.locN; i++) {
    if(spikingStatusInThisTimeStep[i]==1) 
    {
      if(compressedSpikingOffset>=DSD__maxAxonalSpike){
	printf("ERROR in compressNeuralSpikesInThisTimeStep: number of spikes greater then DSD__maxAxonalSpike\n");
	fflush(stdout); exit(0);
      }
      spikingNeuronIdsInThisTimeStep[compressedSpikingOffset]
               = i + lnp_par.loc_h*lnp_par.locN;
      #ifdef LIFCAneuron
      spikingNeuronEmissionTime[compressedSpikingOffset] = n[i].getLastEmissionTime();
      #endif
      DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
      printf("compressed spike of glob_n %d at %d ms\n",
	       spikingNeuronIdsInThisTimeStep[compressedSpikingOffset],
               thisSimTimeStep_ms);
      DPSNNverboseEnd();

      if(compressedSpikingOffset + 1 > spikeCountInThisTimeStep) {
	printf("ERROR while compressing spikes\n"); fflush(stdout); exit(0);
      };
      compressedSpikingOffset ++;
    };
  };

    if(compressedSpikingOffset != spikeCountInThisTimeStep) {
      printf("ERROR in compression of spikes %d != %d\n",
	   compressedSpikingOffset, spikeCountInThisTimeStep);
      fflush(stdout); exit(0);};
  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if(spikeCountInThisTimeStep!=0) {
    printf("compressed a total of %d Neural Spikes at %d ms \n",
	 spikeCountInThisTimeStep, thisSimTimeStep_ms);
  }
  DPSNNverboseEnd();
};

void localNetClass::neuralSpikesToAxonalSpikes() {
  //deubug printf range 652-698
  uint32_t spike_i;
  uint32_t spiking_glob_n,spiking_loc_n;
  double spiking_time;
  uint32_t target_h;  

  uint32_t forwardCount;
  uint32_t totalConnectionsForThisSpike;

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if(spikeCountInThisTimeStep!=0)
     printf(
     "653- neuSpikesToAxSpikes: h=%03d START there are %d spikes at %d ms\n",
       lnp_par.loc_h, spikeCountInThisTimeStep, thisSimTimeStep_ms);
  DPSNNverboseEnd();

  for(target_h=0;target_h<lnp_par.globH;target_h++) {
    forwardAxonalSpikes[target_h].count = 0;
    forwardAxonalSpikes[target_h].expectedCount = 0;
  };

  if(spikeCountInThisTimeStep>=DSD__maxAxonalSpike){
    printf("ERROR in neuralSpikesToAxonalSpikes: number of spikes greater then DSD__maxAxonalSpike\n");
    fflush(stdout); exit(0);
  }

  for (spike_i = 0;
       spike_i < spikeCountInThisTimeStep;
       spike_i++) 
  {
    totalConnectionsForThisSpike=0;
    spiking_glob_n = spikingNeuronIdsInThisTimeStep[spike_i];
    spiking_loc_n  = spiking_glob_n % lnp_par.locN;
    spiking_time = spikingNeuronEmissionTime[spike_i];

    DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
      printf(
         "654- neuSpikToAxSpikes: neuSpikes = %05d EXISTS on h=%03d at %d ms\n",
         spiking_glob_n, lnp_par.loc_h, thisSimTimeStep_ms);
    DPSNNverboseEnd();
    
    for(target_h=0;target_h<lnp_par.globH;target_h++) {
      forwardCount=0;
      if(n[spiking_loc_n].getForwardTargetDistrInHost(target_h) !=0 )
      {
        totalConnectionsForThisSpike++;
	DPSNNverboseStart(false, thisSimTimeStep_ms, lnp_par.debugPrintEnable_ms)
          printf(
         "655- neuSpikToAxSpikes: spike neu=%05d NEEDS s_h%03d->t_h%03d\n at %d ms\n",
          spiking_glob_n, lnp_par.loc_h, target_h, thisSimTimeStep_ms);
	DPSNNverboseEnd();
       forwardCount = forwardAxonalSpikes[target_h].count;
       DPSNNverboseStart(false, thisSimTimeStep_ms, lnp_par.debugPrintEnable_ms)
         printf("656: from h=%d to h=%d neuron=%d forwardCountInList=%d\n",
		lnp_par.loc_h,target_h,spiking_glob_n,forwardCount);
       DPSNNverboseEnd();
       forwardAxonalSpikes[target_h].list[forwardCount].pre_glob_n =
	     spiking_glob_n;

       forwardAxonalSpikes[target_h].list[forwardCount].originalEmissionTime = spiking_time;
       forwardAxonalSpikes[target_h].count++;
      }

    }

    if(totalConnectionsForThisSpike == 0) {
      printf(
      "ERROR neuSpikToAxSp: no forw host targ for neu = %d spiking on h = %03d\n",
      spiking_glob_n, lnp_par.loc_h);
      fflush(stdout);exit(0);
    };

    DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
      printf(
      "656- neuSpikToAxSpikes: neu=%05d TOTAL %d host->host conn required from h=%03d at %d ms\n",
      spiking_glob_n, totalConnectionsForThisSpike, lnp_par.loc_h, thisSimTimeStep_ms);
    DPSNNverboseEnd();
  }
};


void localNetClass::sendReceiveAxonalSpikes() {

chronoSendReceiveAxonalSpikes.startChrono();

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
    printf("!sendReceiveAxonalSpikes: %d ms h=%d START with %d spikes to send\n", 
	    thisSimTimeStep_ms, lnp_par.loc_h, spikeCountInThisTimeStep);
    fflush(stdout);
  DPSNNverboseEnd();

  chronoExchangeAxonalSpikesDim.startChrono();

  pMessagePassing->exchangeAxonalSpikesDim(
    forwardAxonalSpikes, 
    &synTargetHostDistribution,
    backwardAxonalSpikes, 
    &synSourceHostDistribution,
    thisSimTimeStep_ms);

  chronoExchangeAxonalSpikesDim.stopChronoPartial(thisSimTimeStep_ms);

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
    printf("sendReceiveAxonalSpikes: %d ms h=%d AFTER MPI DIM EXCHANGE\n", 
	   thisSimTimeStep_ms, lnp_par.loc_h);
    fflush(stdout);
  DPSNNverboseEnd();

  chronoExchangeAxonalSpikes.startChrono();

  pMessagePassing->exchangeAxonalSpikes(
    forwardAxonalSpikes, 
    backwardAxonalSpikes, 
    thisSimTimeStep_ms);
  chronoExchangeAxonalSpikes.stopChronoPartial(thisSimTimeStep_ms);

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  {  uint32_t s,h;
    axonalSpikeDataOnlyClass axSp;
    for(h=0;h<lnp_par.globH;h++) {
      for(s=0;s<backwardAxonalSpikes[h].count;s++) {
	axSp=backwardAxonalSpikes[h].list[s];
	printf(
	       "sendReceiveAxonalSpikes: %d ms h=%d REC AxSp from glob_n=%d,orig %.7f ms\n",
	       thisSimTimeStep_ms, lnp_par.loc_h, 
	       axSp.pre_glob_n, axSp.originalEmissionTime);
      };
    };
    fflush(stdout);
  }
  DPSNNverboseEnd();

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
    printf("sendReceiveAxonalSpikes: %d ms h=%d END\n", 
	   thisSimTimeStep_ms, lnp_par.loc_h );
    fflush(stdout);
  DPSNNverboseEnd();

  chronoSendReceiveAxonalSpikes.stopChronoPartial(thisSimTimeStep_ms);
};

void localNetClass::printBackwardAxonalSpikes() {
  uint32_t h,c,s;
  for(h=0;h<lnp_par.globH;h++) {
    c = backwardAxonalSpikes[h].count;
    printf("backward Axonal Spikes received at %d ms from h=%d -  total=%d\n",
	   thisSimTimeStep_ms, h, c);
    for(s=0;s<c;s++) {
      printf("ax Spike: source_neu=%d, originalEmissionTime=%f\n",
	     backwardAxonalSpikes[h].list[s].pre_glob_n,
     	     backwardAxonalSpikes[h].list[s].originalEmissionTime);
    }
  }
}

void localNetClass::backwardAxonalSpikesToAxonalSpikeSchedulers() {
  uint32_t source_h;

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
    printf("AxSpikeToAxSpikeSched: START on h=%03d at %d ms\n", 
	   lnp_par.loc_h, thisSimTimeStep_ms);
  DPSNNverboseEnd();

  //fill axonalSpikeScheduler delay line
  //axonalSpikeScheduler.insertAxSpike(lnp_par.globH,backwardAxonalSpikes,thisSimTimeStep_ms);
  for(source_h=0;source_h<lnp_par.globH;source_h++)
    axonalSpikeScheduler[source_h].insertAxSpike(source_h,backwardAxonalSpikes,thisSimTimeStep_ms);

};

void localNetClass::axonalSpikeSchedulersToSynapticSpikeSchedulers() {
  uint32_t source_h;
  uint32_t axonalSpike_index;
  uint32_t lastSearchedSourceNeuron_glob_n;
  uint32_t atLeastOneSynFound;
  int32_t delay;
  uint32_t synIndexToBeSet;

  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
  printf("AxSpikeSchedToSynSpikeSched: START on h=%03d at %d ms\n", 
	 lnp_par.loc_h, thisSimTimeStep_ms);
  DPSNNverboseEnd();

  //critical assuntion axonal spikes ordered first by delay
  //then by source spiking neuron

  atLeastOneSynFound=0;

  // for(delay=0;delay<lnp_par.D;delay++){
  // inverted loop index in order to obtain the same "identical" references
  // when compiled with and without compactSpikeScheduler compiler option
  for(delay=lnp_par.D-1;delay>=0;delay--) {

    //axonalSpikeScheduler.extractAxSpike(lnp_par.globH,delayedBackwardAxonalSpikes,delay);
    for(source_h=0;source_h<lnp_par.globH;source_h++){
      axonalSpikeScheduler[source_h].extractAxSpike(source_h,delayedBackwardAxonalSpikes,delay);
    }

    for(source_h=0;source_h<lnp_par.globH;source_h++) {	  
      DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);    
      if(delayedBackwardAxonalSpikes[source_h].count !=0 ) 
	printf("AxSpikToSynSpike: FROM h=%03d TO h=%03d ARRIVING %d AX spikes at %d ms\n", 
	       source_h, lnp_par.loc_h, 
	       delayedBackwardAxonalSpikes[source_h].count, thisSimTimeStep_ms);
      DPSNNverboseEnd();

      if(synDelaySourceHostDistribution.synCount[delay][source_h]!=0) {
	backwardSynapticSearch[source_h].reset(synDelaySourceHostDistribution.synOffset[delay][source_h],
					       synDelaySourceHostDistribution.synCount[delay][source_h],
					       thisSimTimeStep_ms);
	  
	lastSearchedSourceNeuron_glob_n = 0;
	  
	for (axonalSpike_index = 0;
	     axonalSpike_index < delayedBackwardAxonalSpikes[source_h].count;
	     axonalSpike_index ++) { 
	  axonalSpikeDataOnlyClass tempAxonalSpike;
	  uint32_t sourceNeuron_glob_n;
	  synapseClass tempReturnBackwardSyn;
	  uint32_t tempBackwardSynGlobalOffset;
	  synapticSpikeClass tempSynSpike;
	  uint32_t targetNeu_glob_n, targetNeu_loc_n;
	  bool foundBackwardSyn;
	  uint32_t countFoundBackwardSyn;
	    
	  foundBackwardSyn      = false;
	  countFoundBackwardSyn = 0; 

	  tempAxonalSpike = delayedBackwardAxonalSpikes[source_h].list[axonalSpike_index];
	  sourceNeuron_glob_n = tempAxonalSpike.pre_glob_n;

	  //critical assumption      
	  if(sourceNeuron_glob_n < lastSearchedSourceNeuron_glob_n) {
	    printf("ERROR backwAxSpLst non ordered h=%d: neu=%d < last searc=h%d\n",
		   lnp_par.loc_h, sourceNeuron_glob_n, 
		   lastSearchedSourceNeuron_glob_n);
	    fflush(stdout); exit(0);
	  } else {
	    lastSearchedSourceNeuron_glob_n = sourceNeuron_glob_n;
	  };

	  DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
	  printf("AxSpikToSynSpike: FROM h=%d TO h=%03d axSpIdx=%d origEmissTime=%f sourceNeu=%d on delay=%d BEFORE synSearch\n",
		 source_h, lnp_par.loc_h, axonalSpike_index,
		 tempAxonalSpike.originalEmissionTime,sourceNeuron_glob_n,delay);
	  DPSNNverboseEnd();

	  synIndexToBeSet = backwardSynOffsetInSynList[sourceNeuron_glob_n].offsetByDelay[delay];
	  if(synIndexToBeSet!=0xFFFFFFFF){
	    backwardSynapticSearch[source_h].setSynIndex(synIndexToBeSet);

	    foundBackwardSyn = backwardSynapticSearch[source_h].findAndGetNextSyn(sourceNeuron_glob_n, 
										  &tempReturnBackwardSyn, 
										  backwardSynList, 
										  &tempBackwardSynGlobalOffset,
										  thisSimTimeStep_ms);

	    if(foundBackwardSyn==true){
	      atLeastOneSynFound=1;
	      do {
		countFoundBackwardSyn++;
		tempSynSpike.synIndex = tempBackwardSynGlobalOffset;
		tempSynSpike.synDelayPlusOriginalEmissionTime =
		  (double)tempReturnBackwardSyn.delay +
		  tempAxonalSpike.originalEmissionTime;

		targetNeu_glob_n = tempReturnBackwardSyn.post_glob_n;
		targetNeu_loc_n = targetNeu_glob_n % lnp_par.locN;

		DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
		printf("AxSpikToSynSpike: FROM h=%d TO h=%03d targN=%d synDel=%d futureActivT=%f synIdx=%d BEFORE insertSynSpike\n",
		       source_h, lnp_par.loc_h,
		       targetNeu_glob_n,
		       tempReturnBackwardSyn.delay,
		       tempSynSpike.synDelayPlusOriginalEmissionTime,
		       tempSynSpike.synIndex);
		DPSNNverboseEnd();
#ifdef LIFCAneuron
	        //n[targetNeu_loc_n].addInputSpike(backwardSynList[tempBackwardSynGlobalOffset].weight,
		//				 tempAxonalSpike.originalEmissionTime);
			
		{
		  float weightF;
		  //conversion factor from int16 to float
		  weightF = ((float)(backwardSynList[tempBackwardSynGlobalOffset].weight)) 
		    * lnp_par.factorWeightType_2_Float;
		  n[targetNeu_loc_n].addInputSpike(weightF,tempAxonalSpike.originalEmissionTime);		
		}
		//LTD
                #ifdef makeActiveLTD
		if((thisSimTimeStep_ms>=lnp_par.startPlasticity_ms) &&
		   (thisSimTimeStep_ms<=lnp_par.stopPlasticity_ms)) {
		  n[targetNeu_loc_n].synapseSetLastActiveTime_antiCausalSTDP(&backwardSynList[tempBackwardSynGlobalOffset],thisSimTimeStep_ms);
		}
                #endif
		
#else
		//Current injection
		n[targetNeu_loc_n].addInputCurrent(backwardSynList[tempBackwardSynGlobalOffset].weight);
		//LTD
            #ifdef makeActiveLTD
		if((thisSimTimeStep_ms>=lnp_par.startPlasticity_ms) &&
		   (thisSimTimeStep_ms<=lnp_par.stopPlasticity_ms)){
		  n[targetNeu_loc_n].synapseSetLastActiveTime_LTD(&backwardSynList[tempBackwardSynGlobalOffset],thisSimTimeStep_ms);
		}
            #endif
#endif
	      } while( true == backwardSynapticSearch[source_h].findAndGetNextSyn(
										  sourceNeuron_glob_n, 
										  &tempReturnBackwardSyn,
										  backwardSynList, 
										  &tempBackwardSynGlobalOffset,
										  thisSimTimeStep_ms));

	      DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
	      printf("AxSpikToSynSpike: FROM h=%d TO h=%03d TOTAL %d syn found from neu %d\n",
		     source_h, lnp_par.loc_h, 
		     countFoundBackwardSyn, 
		     sourceNeuron_glob_n); fflush(stdout);
	      DPSNNverboseEnd();
	    }

	    if(atLeastOneSynFound==0 && delay==0 && source_h==lnp_par.globH){
	      printf("ERROR AxSpikToSynSpike: on h=%03d, %d ms not found any syn from neu %d\n",
		     lnp_par.loc_h, thisSimTimeStep_ms, sourceNeuron_glob_n);
	      fflush(stdout);exit(0);
	    }
	  }
	}// end of loop on "axonalSpike_index"
      }// end of if on synDelaySourceHostDistribution.synCount[][]
    }// end of loop on "source_h"
  }// end of loop on "delay"
}

bool localNetClass::isChronoWindow(uint32_t thisSimTimeStep_ms){
  if((thisSimTimeStep_ms >= startPartialChrono_ms) && (thisSimTimeStep_ms < stopPartialChrono_ms))
    return true;
  else
    return false;
}

void localNetClass::firingsInChronoWindow(uint32_t spikeCountInThisTimeStep) {
  double delta;

  N_firingsInChronoWindow += spikeCountInThisTimeStep;
  numLapsInChronoWindow ++;
  minFiringsInChronoWindow = (spikeCountInThisTimeStep < minFiringsInChronoWindow) ? spikeCountInThisTimeStep : minFiringsInChronoWindow;
  maxFiringsInChronoWindow = (spikeCountInThisTimeStep > maxFiringsInChronoWindow) ? spikeCountInThisTimeStep : maxFiringsInChronoWindow;
  delta = spikeCountInThisTimeStep - meanFiringsInChronoWindow;
  meanFiringsInChronoWindow += (double)(delta / numLapsInChronoWindow);
  M2FiringsInChronoWindow += (double)((double)delta * (double)((double)spikeCountInThisTimeStep - (double)meanFiringsInChronoWindow));
  varianceFiringsInChronoWindow = (double)(M2FiringsInChronoWindow / (numLapsInChronoWindow -1));
  sigmaFiringsInChronoWindow = (double)(sqrt(varianceFiringsInChronoWindow));
  coeffOfVariationFiringsInChronoWindow = (double)(sigmaFiringsInChronoWindow / meanFiringsInChronoWindow);
};

void localNetClass::printStatFiringsInChronoWindow(uint32_t loc_h){
  FILE *fp_output;
  char reportName[80];
  double firingRate;

  sprintf(reportName,"FiringsInChronoWindow_h%d.stat",loc_h);
  fp_output=fopen(reportName,"w");

  firingRate = double((double)N_firingsInChronoWindow * 1000) / double((double)numLapsInChronoWindow * (double)lnp_par.locN);

  fprintf(fp_output,"FIRING-h=%.3d     Firings in ms %6d-%6d = %-11d min = %-11d max = %-11d mean = %-11f variance = %-11f sigma = %-11f coeffOfVariation = %-11f samples = %d\n",loc_h,startPartialChrono_ms,stopPartialChrono_ms,N_firingsInChronoWindow,minFiringsInChronoWindow,maxFiringsInChronoWindow,meanFiringsInChronoWindow,varianceFiringsInChronoWindow,sigmaFiringsInChronoWindow,coeffOfVariationFiringsInChronoWindow,numLapsInChronoWindow);
fprintf(fp_output,"FIRING-h=%.3d  FiringRate in ms %6d-%6d = %-11f\n",loc_h,startPartialChrono_ms,stopPartialChrono_ms,firingRate);
  fclose(fp_output);
}
