// DPSNN_neuron.h
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

#ifndef DPSNN_neuronIncluded
#define DPSNN_neuronIncluded

#include "DPSNN_random.h"
#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"
#include "DPSNN_spike.h"
#include "DPSNN_stat.h"
#ifdef LIFCAneuron
//#include "DPSNN_LIFCAconnectome.h"
#include "erflib.h"
#include "randdev.h"
#endif

class neuronClass { // individual neuron database

public:
    float InputCurrent; // (total) current incoming from backward synapses
    uint32_t getGlob_n();
    neuSubPopEnum getSubPop_n();
    uint32_t get_hashId();
    uint32_t set_glob_n(uint32_t glob_n_initValue);
    uint32_t set_loc_n(uint32_t loc_n_initValue);
    uint32_t set_loc_h(uint32_t loc_h_initValue);
    bool didItFire;
    void set_hashId();
    uint32_t forwardNeuralTargetDistrInHost[DSD__maxGlobH];
    //initialization phase
    void initGlobalAwareness(const DPSNN_parameters lnp_par_initValue);
    //void initNeuralId(const uint32_t loc_n_initValue);
    void initNeuralKind(const enum neuralKindEnum initKind,
			const enum neuSubPopEnum subPop,
			const neuSubPopParamStruct *neuSubPopParam);
    void initM(uint32_t M_initValue);//with check for range
    void initD(uint32_t D_initValue);//with check for range
    void clearForwardConnections(); //sets an impossible value
    void initStat(statClass *pStat);

    void initLongTermPlasticity();

    uint32_t getForwardTargetDistrInHost(const uint32_t h);
    void initRandDevPointer(randdevClass *pLocalNetRandDev_initValue);
  
    void clearInputCurrent();
    void addInputCurrent(const float inputValue);
    void setLastEmittedSpikeTime_ms(const int32_t inputValue);
    void spikeLogging(const uint32_t spikeTime, const uint32_t neuralId);

    //simulation of dynamic steps
    void synapseSetLastActiveTime_antiCausalSTDP(synapseClass *pointerActiveSynapse, const uint32_t thisTimeStep_ms);
    void regularDynamicTimeStep(const uint32_t thisTimeStep_ms, const float *pInputDistrib, uint32_t * pThalInputCounter); 
    //bool didItFire();
    //void afterSpikeNeuralDynamic(const int32_t thisTimeStep_ms);
    float causalSTDP_ms(const int32_t localNeuralSpikeTime_ms,
	    const int32_t previousSynapticActivationTime_ms);
    void causalSTDP_ms_ofBackwardSynapses(const uint32_t thisNeuronLastSpike_ms);
    float antiCausalSTDP_ms(const int32_t synapticActivationTime_ms);
public:
    //initializations
    void initForwardWeightsAndDerivatives();//synapses are set according to the neural kind
#ifdef LIFCAneuron   
  void setTableLUT_synBath(double *tableLUT_synBath_Exc,double *tableLUT_synBath_Inh);
  void setTableLUT_synBathMatrix(uint64_t address_of_tableLUT_synBathMatrix);
#endif
    //initial generation of synaptic details
    //e.g.  target neurons, delays and neurons
private:
    void checkSizeOfBackwardSynList();
private:	
  //PSP There are 4 parameters in the Izhikevich equations...
  //...for each neuron: a, b, c, d. 
  //PSP Setting those parameters you obtain a different
  //PSP kind of neuron (there are 20 types of actual neurons)
  //PSP c is set to -65 mv for all neurons, directly in the time update routine
  //PSP b is set to 0.02 for all neurons, directly in the time update routine
  //PSP there are at least two types of neurons in this model:
  //PSP cortical pyramidal neurons exhibiting regular spiking (RS) 
  //PSP firing pattern
  //PSP The Izhikevic equation used in the code is described ...
  //... at page 276 of [IzhPol]
#ifdef LIFCAneuron
  double Tau;                    //Decay time constant for the membrane potential
  double Theta;                  //The emission threshold
  double H;                      //The reset potential
  double Tarp;                   //The absolute refractory period
  double AlphaC;                 //Increase of Ca concentration after the emission of a spike
  double TauC;                   //Characteristic time of Ca channel inactivation
  double gC;                     //K-current due to a unitary Ca concentration
  //uint16_t InitType;

  double Tr;                   //Arriving time of the last pre-synaptic spike
  double Te;                   //Emission time of the last spike
  double assignedExternalSpikeTime;
  double *pTableLUT_synBath;

  double v, c;			// activity variables

  void setLIFCAParams(int NumRealParams, float *RealParams, 
	       int NumStringParams, char **StringParams);
  void LIFCADynamic(const uint32_t thisTimeStep_ms, const double currentTime, const double currentValue);
  double getPoissonLapse();
  double getUniformLapse();
  randdevClass *pLocalNetRandDev;
  void sortCurrentsInPlace(inputCurrentArrayStruct *array, uint32_t elemN);
  void mergeCurrentsArrays(inputCurrentArrayStruct *arrayOut, inputCurrentArrayStruct *array1, inputCurrentArrayStruct *array2, uint32_t count1, uint32_t count2);
  void mergeFunction(inputCurrentArrayStruct *A, uint32_t p, uint32_t q, uint32_t r);
  void mergeSort(inputCurrentArrayStruct *A, uint32_t p, uint32_t r);

 public:  
  double invNuExt;
  uint32_t maxExternCurrentCount;
  double getLastEmissionTime();
  float inputCurrentPlus;
  float inputCurrentMinus;
  uint32_t inputCurrentPlusCount;
  uint32_t inputCurrentMinusCount;
  //EPA - 2015-03-10
  //inputCurrentClass inputCurrents[DSD__maxSimultaneousSpikesOnSameTarget];
  inputCurrentClass *inputCurrents;
  uint32_t inputCurrentsCount;
  void addInputSpike(const float current,const double time);

 private:
#else
  float	a, d;			// neuronal dynamics parameters 
  float	v, u;			// activity variables
#endif
  enum neuralKindEnum kind;	// neuron kind (excitatoryRS, inhibitoryFS, ...,
  enum neuSubPopEnum subPop;    // in LIFCA is the subPopulation

  int32_t lastNeuralEmittedSpikeTime_ms;
  //we keep a local copy of the synapses, 
  // This could be release after creation and dump
  //into the general list of localNet

  synapseClass *pointBackwardSynList;
 
  //knowledge of parameters of global system and local cluster of neurons

  DPSNN_parameters lnp_par;
  statClass *pStat;
  //lnp_par.globH knowledge of total number of hosts in the system	
  //lnp_par.globN knowledge of total number of neurons in the system
  //lnp_par.locN and lnp_par.log2_locN knowledge of total number of... 
  //... neurons in the local cluster
  //locNe knowledge of number of excitatory neurons in the local cluster
  //lnp_par.loc_h  local host on which the neuron resides

  uint32_t loc_n;  // local neuron identifier <lnp_par.locN
  uint32_t loc_h;  // identifier of the host hosting the neuron <lnp_par.globH
  uint32_t glob_n; // global neuron identifier <lnp_par.globN
  uint32_t hashId;  // hashed neuron index
  uint32_t M;	   // number of synapses in output
  uint32_t D;	   // max possible value for synaptic delay in ms
 
  // PSP Excitatory synapses are plastic and evolve 
  // according to STDP (Synaptic Time Dependent Plasticity)
  // PSP Inhibitory synapses are not plastic

  // PSP the simulation proceeds in chunks of 1000 ms steps. ...
  // ... Sometimes you need informations from the past...
  // PSP due to the delays introduced by axons/synapses. 
  // PSP Therefore certain vectors need D additional locations
  // PSP STDP (Synaptic Time Dependent Plasticity) functions 

  // PSP N_pre: how many synapses arriving to the neuron
  // PSP N_pre: can ne larger than M, 
  // because the outgoing synapses are assigned random ...
  // ... and therefore more than M could converge on a target
  // PSP arbitrary assumption: 
  // 3*M the max number of incoming synapses. 
  uint32_t N_pre;//how many synapses arriving to the local neuron
public: 
#if defined(makeActiveLTD) || defined (makeActiveLTP)
  uint32_t backwardSynOffset[DSD__maxM_pre]; 
    //initial generation of synapses
    void addBackwardSyn(const uint32_t synOffset);
    void clearBackwardSynList();//clears the number of incoming synapses
    uint32_t getSizeBackwardSynList() {return(N_pre);};
    void initBackwardSynListPointer(synapseClass *pointBackwardSynList);
#endif
  enum neuralKindEnum getNeuralKind(){return(kind);}
};  //end of declaration of class representing individual neuron //

#endif //multiple inclusion guard

