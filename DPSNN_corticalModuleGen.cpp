// DPSNN_corticalModuleGen.cpp
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_debug.h"
#include "DPSNN_neuron.h"
#include "DPSNN_localNet.h"
#include "DPSNN_connectome.h"

void localNetClass::simpleCM_prepareForwardSynapses() {
  chronoClass chronoDescribeConnectome;  
  chronoClass chronoInvokeIndivSynGen;
  DPSNNverboseStart(true,1,0);

    chronoDescribeConnectome.clearAndStartChrono();
  DPSNNverboseEnd();

  simpleCM_connectome.describeConnectome(&lnp_par);
  
  DPSNNverboseStart(true,1,0);
    chronoDescribeConnectome.stopChrono(); 
    if(lnp_par.loc_h <= 1 ||  lnp_par.loc_h>=(lnp_par.globH-2)) {
      printf("CHRONO:describeConnectome h=%d = %f sec \n",
        lnp_par.loc_h, chronoDescribeConnectome.getAccumulatedChrono());
      fflush(stdout);};
  DPSNNverboseEnd();

  DPSNNverboseStart(true,1,0);
    chronoInvokeIndivSynGen.clearAndStartChrono();
  DPSNNverboseEnd();

  simpleCM_invokeIndivSynGenerators();

  DPSNNverboseStart(true,1,0);
    chronoInvokeIndivSynGen.stopChrono(); 
    if(lnp_par.loc_h <= 1 ||  lnp_par.loc_h>=(lnp_par.globH-2)) {
      printf("CHRONO:invokeIndivSynGen h=%d = %f sec \n",
        lnp_par.loc_h, chronoInvokeIndivSynGen.getAccumulatedChrono());
      fflush(stdout);};
  DPSNNverboseEnd();
};


void localNetClass::simpleCM_invokeIndivSynGenerators()
{
  uint32_t i;
  uint32_t jSynIdInNeu;
  uint32_t iTot,jTot;
  simpleCM_neuCoordinatesStruct sourceNeuIds;
  simpleCM_neuCoordinatesStruct targetNeuIds;

  iTot=0;jTot=0;

  for(i=0;i < lnp_par.locN; i++) {
    sourceNeuIds = 
      simpleCM_connectome.convert_loc_n_h_to_neuCMCoordinates(
      i,lnp_par.loc_h);
    srand(sourceNeuIds.glob_n);
    simpleCM_connectome.seedForSynapses = sourceNeuIds.glob_n;
    n[i].initGlobalAwareness(lnp_par);
    n[i].initStat(pStat);
    n[i].set_glob_n(sourceNeuIds.glob_n);
    n[i].set_loc_n(i);
    n[i].set_loc_h(lnp_par.loc_h);
    n[i].initM(lnp_par.M);
    n[i].initD(lnp_par.D);
    n[i].clearForwardConnections();
    n[i].initNeuralKind(sourceNeuIds.neuralKind); 
 
    simpleCM_connectome.countRandSynGen = 0;
    for(jSynIdInNeu=0;jSynIdInNeu < lnp_par.M;jSynIdInNeu++) 
      simpleCM_connectome.synListOfThisNeu [jSynIdInNeu] = lnp_par.globN; //initialized to an absurd value

    for(jSynIdInNeu=0;jSynIdInNeu < lnp_par.M;jSynIdInNeu++) 
    {
      targetNeuIds = simpleCM_connectome.generateTargetNeu(sourceNeuIds,jSynIdInNeu);

      //final section: synapse filling 
      forwardSynList[i*lnp_par.M+jSynIdInNeu].pre_glob_n = 
	      sourceNeuIds.glob_n;
      forwardSynList[i*lnp_par.M+jSynIdInNeu].post_glob_n = 
	      targetNeuIds.glob_n;
      n[i].forwardNeuralTargetDistrInHost[
	      targetNeuIds.loc_h] ++;

      switch (sourceNeuIds.neuralKind) {
      case excitatoryRS:
        forwardSynList[i*lnp_par.M+jSynIdInNeu].preSynNeuralKind=excitatoryRS; 
        break;
      case inhibitoryFS:
        forwardSynList[i*lnp_par.M+jSynIdInNeu].preSynNeuralKind=inhibitoryFS; 
        break;
      default:
        printf(
        "Error in initOutputWeightsAndDerivatives() - unknown neuron kind\n");
        fflush(stdout);exit(0);
        break;
      };
    };

    //generation of delays
    if( (lnp_par.M / lnp_par.D) * lnp_par.D != lnp_par.M) {
        printf("ERROR lnp_par.M=%d is not a multiple of lnp_par.D=%d\n",
	     lnp_par.M, lnp_par.D);fflush(stdout);exit(0);
    };
    {uint32_t m,d;
    //generation of delays
      switch (sourceNeuIds.neuralKind) {
      case excitatoryRS:
        //uniform distribution of exc. synaptic delays
        for(d=0; d < lnp_par.D; d++) {
          for (m =    d  * lnp_par.M/lnp_par.D; 
             m < (d+1) * lnp_par.M/lnp_par.D ; 
             m ++) {
            forwardSynList[i*lnp_par.M+m].delay = d;
          };
        };
        break;
      case inhibitoryFS:
#ifdef LIFCAneuron
        // uniform distribution of exc. synaptic delays in LIFCA
	// NOTE: this trial didn't give good results,
	// a different distribution of the synapses delays of inhib. neurons mus be tested
	// At the moment it's used delay=0 as in Izhikevich model
	
        for(d=0; d < 10; d++) {
          for (m =    d  * lnp_par.M/10; 
             m < (d+1) * lnp_par.M/10 ; 
             m ++) {
            forwardSynList[i*lnp_par.M+m].delay = d;
          };
        };

        // all inhibitory delays are 1 ms
        //for (m=0; m<lnp_par.M; m++) 
	//  forwardSynList[i*lnp_par.M+m].delay = 0;
#else
        // all inhibitory delays are 1 ms
        for (m=0; m<lnp_par.M; m++) 
	  forwardSynList[i*lnp_par.M+m].delay = 0;
#endif
        break;
      default:
        printf(
	"Error in neuron::generateForwardDelays() - unknown neu kind\n");
        fflush(stdout);exit(0);
        break;
      };
    };	    
  };

  //END loop on all source neurons 

  DPSNNverboseStart(false,1,0);
    printf(
    "loc_h=%d - iTot=%d  expected %d , jTot=%d expected %d \n", 
    lnp_par.loc_h, iTot, lnp_par.globN, jTot, lnp_par.globN * lnp_par.M);
    printf(
    "loc_h=%d - CFT=%d, CFX=%d,CFY=%d,neuPerCM=%d\n", 
    lnp_par.loc_h, lnp_par.globCFT, lnp_par.globCFX, 
    lnp_par.globCFY, lnp_par.neuronsPerCM);
  DPSNNverboseEnd();
};

