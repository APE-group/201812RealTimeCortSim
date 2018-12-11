// DPSNN_neuron_init.c
// DPSNN_*.* Distribution/Parallelization 
// performed by Pier Stanislao Paolucci (Roma, Italy project start date 2011),
// starting from the sequential code
// SPNET.cpp: Spiking neural network with axonal conduction delays and STDP
// Created by Eugene M. Izhikevich, May 17, 2004, San Diego, CA
 
// reference paper, named [IzhPol] in the following 1// Eugene M. Izhikevich "Polychronization: Computation with Spikes" 
// Neural Computation 18, 245-282 (2006)
 
// Saves spiking data each second in file spikes.dat
// To plot spikes, use MATLAB code: load spikes.dat;plot(spikes(:,1),spikes(:,2),'.');
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "DPSNN_debug.h"
#include "DPSNN_random.h"

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"
#include "DPSNN_neuron.h"
#include "DPSNN_stat.h"

uint32_t neuronClass::getForwardTargetDistrInHost(const uint32_t h) {
  return (forwardNeuralTargetDistrInHost[h]);
}

uint32_t neuronClass::set_glob_n(const uint32_t n) {
  glob_n=n;
  return 0;
};

uint32_t neuronClass::set_loc_h(const uint32_t h) {
  loc_h=h;
  return 0;
};

uint32_t neuronClass::set_loc_n(const uint32_t n) {
  loc_n=n;
  return 0;
};

uint32_t neuronClass::getGlob_n() {
  return (glob_n);
};

neuSubPopEnum neuronClass::getSubPop_n() {
  return (subPop);
};

void neuronClass::initGlobalAwareness(
  const struct DPSNN_parameters lnp_par_initValue) {
  lnp_par=lnp_par_initValue;
};

void neuronClass::initStat(
  statClass *pStat_initValue) {
  pStat = pStat_initValue;
};

/* EPA 13-Nov-2015: Currently not used */
/*
void neuronClass::initNeuralId(const uint32_t loc_n_initValue) {
  loc_n = loc_n_initValue;
  switch (kind) {
  case excitatoryRS:
  case excitatoryLbExc:
  case excitatoryLaExc:
    glob_n = loc_n + lnp_par.loc_h * lnp_par.locN;
    break;
  case inhibitoryFS:
  case inhibitoryLaInh:
    glob_n = loc_n + lnp_par.loc_h * lnp_par.locN;
    break;
  default:
    kind=unknown;
    printf("Error::initNeuralKind() - unknown neuronKind\n");
    fflush(stdout);exit(0);
    break;
  };

  DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
    printf("initNeuralId glob_n=%d loc_n=%d loc_h%d completed\n",
	 glob_n,loc_n,lnp_par.loc_h);fflush(stdout);
  DPSNNverboseEnd();

  if (glob_n>=lnp_par.globN){
    printf("ERROR: neuron::initId() glob_n out of range\n");
    fflush(stdout);exit(0);
  }      
};
*/

#ifdef LIFCAneuron
void neuronClass::initNeuralKind(const enum neuralKindEnum initKind,
				 const enum neuSubPopEnum subPopulation,
				 const neuSubPopParamStruct *neuSubPopParam) {
  kind=initKind;
  subPop=subPopulation;
  Tau=neuSubPopParam[subPop].Tau;
  Theta=neuSubPopParam[subPop].Theta;
  H=neuSubPopParam[subPop].H;
  Tarp=neuSubPopParam[subPop].Tarp;
  AlphaC=neuSubPopParam[subPop].AlphaC;
  TauC=neuSubPopParam[subPop].TauC;
  gC=neuSubPopParam[subPop].gC;
  v=H; c=0.0; 
  Te=-1000;
  assignedExternalSpikeTime=0.0;//getPoissonLapse();
  
};
void neuronClass::setTableLUT_synBath(
  double *tableLUT_synBath_Exc,double *tableLUT_synBath_Inh) {
  switch (kind) {
  case excitatoryLbExc:
  case excitatoryLaExc:
    pTableLUT_synBath=tableLUT_synBath_Exc;
    break;
  case inhibitoryLaInh:
    pTableLUT_synBath=tableLUT_synBath_Inh;
    break;
  default:
    kind=unknown;
    printf("Error::initNeuralKind() - unknown neuronKind\n");
    fflush(stdout);exit(0);
    break;
  };
};
void neuronClass::setTableLUT_synBathMatrix(
  uint64_t address_of_tableLUT_synBathMatrix) {
  pTableLUT_synBath=(double *)address_of_tableLUT_synBathMatrix;
};

void neuronClass::initRandDevPointer(randdevClass *pRandDevPoint) {
  pLocalNetRandDev = pRandDevPoint; 
};


#else
void neuronClass::initNeuralKind(const enum neuralKindEnum initKind) {
  switch (initKind) {
  case excitatoryRS:
    kind=initKind;
    a=0.02; d=8.0; v=-65.0; u=0.2*v; 
    break;
  case inhibitoryFS:
    kind=initKind;
    a=0.1; d=2.0; v=-65.0; u=0.2*v; 
    break;
  default:
    kind=unknown;
    printf("Error::initNeuralKind() - unknown neuronKind\n");
    fflush(stdout);exit(0);
    break;
  };
};

#endif

#if defined(makeActiveLTD) || defined (makeActiveLTP)
void neuronClass::initBackwardSynListPointer(synapseClass *synListPoint) {
  pointBackwardSynList = synListPoint; 
};

void neuronClass::addBackwardSyn(const uint32_t synOffset) {
  backwardSynOffset[N_pre]=synOffset; 
  N_pre++;//incrementing the counter of total incoming synapses to this neuron 
  checkSizeOfBackwardSynList();//checks if space exists for further synapses
};

void neuronClass::clearBackwardSynList() {
  //number of incoming synapses
  N_pre=0;
};

void neuronClass::checkSizeOfBackwardSynList() {
  if(N_pre>=DSD__maxM_pre) {
    printf(
    "ERROR too many syn incoming on neuron glob_n=%d (on loc_h=%d)\n", 
	 glob_n, lnp_par.loc_h);
    fflush(stdout);exit(0);
  };
};
#endif

void neuronClass::initM(uint32_t M_initValue) {
  if (M_initValue<=DSD__maxM) { M=M_initValue;
  }else{ printf("ERROR: neuron::initM() M_initValue out of range\n");fflush(stdout);exit(0);} 
};

void neuronClass::initD(uint32_t D_initValue) {
  if (D_initValue<=DSD__maxD) { D=D_initValue;
  }else{ printf("ERROR: neuron::initD() D_initValue out of range\n");fflush(stdout);exit(0);}; 
};

void neuronClass::setLastEmittedSpikeTime_ms(const int32_t initValue) {
  lastNeuralEmittedSpikeTime_ms=initValue;
};

void neuronClass::clearForwardConnections() {
  uint32_t h;
  //at the beginning the neurons has no target host
  for (h=0;h<lnp_par.globH;h++) {
   forwardNeuralTargetDistrInHost[h]=0;
  };
};

void neuronClass::set_hashId(){
  uint32_t target_glob_n, target_loc_n;
  uint32_t directNeuIndex,reverseNeuIndex,haschedNeuIndex;
  uint32_t xorNeuIndex,maskNeuIndex;
  uint32_t order;
  uint32_t v;
  int32_t r;

  if(lnp_par.globN <= lnp_par.neuronsPerCM){
     hashId = glob_n;
     return;
  }

  target_glob_n = glob_n;
  target_loc_n = target_glob_n % lnp_par.locN;
  directNeuIndex = target_glob_n;

  if (lnp_par.globN > lnp_par.neuronsPerCM) {
    order = 0;
    v = directNeuIndex;
    while (v >>= 1)
      order++;
  }

  reverseNeuIndex = 0;
  maskNeuIndex = 0;
  for(v=directNeuIndex;v;v>>=1) {
    reverseNeuIndex <<= 1;
    reverseNeuIndex |= v & 1;
  }
  xorNeuIndex = (directNeuIndex ^ reverseNeuIndex);

  r = order + 1 - 10;
  if(r >= 0 )
    maskNeuIndex = (1 << r) - 1;
	     
  haschedNeuIndex = directNeuIndex ^ ((directNeuIndex ^ xorNeuIndex) & maskNeuIndex);
  haschedNeuIndex &= (lnp_par.neuronsPerCM - 1);

  DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
  if (directNeuIndex==116543)
    printf("addThalamicInput: random generation - H=%d directNeuIndex=%d reverseNeuIndex=%d order=%d xorNeuIndex=%d maskNeuIndex=%d haschedNeuIndex=%d \n",lnp_par.loc_h,directNeuIndex,reverseNeuIndex,order,xorNeuIndex,maskNeuIndex,haschedNeuIndex);
  DPSNNverboseEnd();

  hashId = haschedNeuIndex;

}

uint32_t neuronClass::get_hashId(){
  return hashId;
}
