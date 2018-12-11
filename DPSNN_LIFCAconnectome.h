// DPSNN_LIFCAconnectome.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
#ifndef DPSNN_LIFCAconnectomeIncluded
#define DPSNN_LIFCAconnectomeIncluded

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "randdev.h"

class simpleCM_connectomeClass {
public:
  uint32_t countRandSynGen;
  uint32_t synListOfThisNeu[DSD__maxM];
  simpleCM_targetNeuListStruct targetNeuList[DSD__maxM];
  uint32_t neuSubPopInCM_count[DSD__subPopTot];
  uint32_t neuSubPopInCM_offset[DSD__subPopTot];
  uint32_t seedForSynapses;

  void describeConnectome(struct DPSNN_parameters *p_lnp_par, struct neuSubPopParamStruct *neuSubPopParam);
  simpleCM_neuCoordinatesStruct convert_loc_n_h_to_neuCMCoordinates(
    uint32_t loc_n, uint32_t loc_h,struct DPSNN_parameters *lnp_par);
  uint32_t generateTargetNeuList(const simpleCM_neuCoordinatesStruct sourceNeuIds,
				 simpleCM_targetNeuListStruct *targetNeuList,
				 struct DPSNN_parameters *lnp_par);
  float generateSourceNeuCExt(uint32_t cfx, uint32_t cfy, struct DPSNN_parameters *lnp_par);
  char* getSubPopName(enum neuSubPopEnum neuSubPop);
  neuSubPopEnum getSubPopEnum(uint32_t neuSubPop);
  void initLocalNetRandDevPointer(randdevClass *pRandDevPoint);
  float getConnectProb(uint32_t sourceSubPop,uint32_t targetSubPop,uint32_t stencilX,uint32_t stencilY);

private:
  //struct DPSNN_parameters *p_lnp_par;
  randdevClass  *pLocalNetRandDev;
  
  float specConnectivityMatrix[DSD__subPopTot][DSD__subPopTot][DSD__stencilX_Max][DSD__stencilY_Max];
  
  float probConnectMatrix[DSD__subPopTot][DSD__subPopTot][DSD__stencilX_Max][DSD__stencilY_Max];
  uint32_t connectMatrix[DSD__subPopTot][DSD__subPopTot][DSD__stencilX_Max][DSD__stencilY_Max];
  uint32_t bathEfficacyTemplate[DSD__bathEfficacyTemplateX_Max][DSD__bathEfficacyTemplateY_Max];

  void initConnectivityParam(neuSubPopParamStruct *neuSubPopParam, struct DPSNN_parameters *lnp_par);
  void initProbConnectMatrix(struct DPSNN_parameters *lnp_par);
  void initConnectMatrix(struct DPSNN_parameters *lnp_par);
  void initBathEfficacyTemplate(struct DPSNN_parameters *lnp_par);
  void initSubPopParamStruct(neuSubPopParamStruct *neuSubPopParam, struct DPSNN_parameters *lnp_par);
  uint32_t generateTargetNeu(uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu, uint32_t subPopId, struct DPSNN_parameters *lnp_par);  
  neuSubPopEnum getNeuralSubPop(uint32_t iNeuInCM,struct DPSNN_parameters *lnp_par );
  neuralKindEnum getNeuralKind(uint32_t iNeuInCM, struct DPSNN_parameters *lnp_par);
  uint32_t conv_neuIdInCM_to_glob_n(uint32_t iNeuIdInCM, uint32_t icfx,
     uint32_t icfy, struct DPSNN_parameters *lnp_par);
  void report();
  uint32_t getEmptyRandomSynapses (uint32_t sourceSubPop, uint32_t targetSubPop, uint32_t X, uint32_t Y,uint32_t maxN);

};

#endif
