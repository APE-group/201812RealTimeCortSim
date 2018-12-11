// DPSNN_connectome.h
// DPSNN_*.*  Distribution/Parallelization 
//of polychronous Spiking Neural Networks
// performed by Pier Stanislao Paolucci 
// (Roma, Italy, project start date 2011)
#ifndef DPSNN_connectomeIncluded
#define DPSNN_connectomeIncluded

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_neuron.h" 

enum neuSubPopEnum {RS,FS,neuSubPopTotal};
enum synFiberEnum {
  intraCM_RS_RS,intraCM_RS_FS,intraCM_FS_RS,
  interCM_xp,interCM_xm,interCM_yp,interCM_ym,
  interCM_xp_yp,interCM_xp_ym,interCM_xm_yp,interCM_xm_ym,
  interCM_2xp,interCM_2xm,interCM_2yp,interCM_2ym,
  synFiberTotal};

struct deltaCMStruct {
  int32_t cfx;
  int32_t cfy;
};



struct simpleCM_neuCoordinatesStruct {
  uint32_t inCM_n;
  uint32_t cfx_n;
  uint32_t cfy_n;
  uint32_t glob_n;
  uint32_t loc_n;
  uint32_t loc_h;
  neuSubPopEnum subPop;
  neuralKindEnum neuralKind;
};

class simpleCM_connectomeClass {
 private:
  //numerosity (perMil) of synapses projected by RS neurons
  uint32_t intraCM_RS_FS_probPerMil;
  uint32_t intraCM_RS_RS_probPerMil;
  uint32_t interCM_RS_1stNeighb_probPerMil;
  uint32_t interCM_RS_2ndNeighb_probPerMil;
  uint32_t interCM_RS_3rdNeighb_probPerMil;

  //numerosity of synapses projected by RS neurons
  uint32_t intraCM_FS_RS_probPerMil;

  void initProbPerMil();

  uint32_t *pRandTable;

public:
  uint32_t countRandSynGen;
  uint32_t synListOfThisNeu[DSD__maxM];
  simpleCM_neuCoordinatesStruct convert_loc_n_h_to_neuCMCoordinates(
			        uint32_t loc_n, uint32_t loc_h);

  simpleCM_neuCoordinatesStruct generateTargetNeu(
				const simpleCM_neuCoordinatesStruct sourceNeuIds,
				uint32_t jSynIdInNeu);

  void describeConnectome(struct DPSNN_parameters *p_lnp_par);

  void initRandTable(uint32_t *pRandTable_initValue);

  uint32_t neuSubPopInCM_count[neuSubPopTotal];
  uint32_t neuSubPopInCM_offset[neuSubPopTotal];

  uint32_t seedForSynapses;

private:
  struct DPSNN_parameters *p_lnp_par;
  //indexes
  //YES intraCM / interCM
  //NOT YET NECESSARY sourceCMid: could be necessary in the future
  //source neural sub pop name
  //target neural sub pop name
  uint32_t neuSubPopInCM_probPerMil[neuSubPopTotal];

  uint32_t synFiber_probPerMil[neuSubPopTotal][synFiberTotal];
  uint32_t synFiber_count[neuSubPopTotal][synFiberTotal];
  uint32_t synFiber_offset[neuSubPopTotal][synFiberTotal]; 

deltaCMStruct deltaCM[synFiberTotal];
  
synFiberEnum getSynFiberType(
  uint32_t synIdInSourceNeu, neuSubPopEnum sourceNeuSubPop);

bool isInterCMSynFiber(synFiberEnum synFiberKind);

bool isIntraCMSynFiber(synFiberEnum synFiberKind);

bool isThisSynFiberKind(
  uint32_t jSynIdInNeu, neuSubPopEnum neuSubPopInCM, 
  synFiberEnum synFiberKind);

neuSubPopEnum getNeuralSubPop(uint32_t iNeuInCM);

neuralKindEnum getNeuralKind(uint32_t iNeuInCM);

uint32_t conv_neuIdInCM_to_glob_n(
  uint32_t iNeuIdInCM, uint32_t icfx, uint32_t icfy); 

simpleCM_neuCoordinatesStruct
  intraModule_targetNeu(
    uint32_t mId, 
    synGenEnum synGenKind, 
    simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct
  interModule_targetNeu(
    uint32_t mId,
    synFiberEnum synFiberKind,
    synGenEnum synGenKind,
    simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct 
  intraModule_random_targetNeu(
    uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct 
  interModule_random_targetNeu(
    uint32_t mId,  
    synFiberEnum synFiberKind, 
    simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct 
  intraModule_nonRand_targetNeu(
    uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct 
  interModule_nonRand_targetNeu(
    uint32_t mId,  
    synFiberEnum synFiberKind, 
    simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct 
  intraModule_randTable_targetNeu(
    uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu);

simpleCM_neuCoordinatesStruct 
  interModule_randTable_targetNeu(
    uint32_t mId,  
    synFiberEnum synFiberKind, 
    simpleCM_neuCoordinatesStruct sourceNeu);

 void report();
};

#endif
