// DPSNN_connectome.cpp
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include "DPSNN_connectome.h"
#include "DPSNN_neuron.h"
#include "DPSNN_debug.h"

void simpleCM_connectomeClass::initProbPerMil() {

  //relative numerosity of neural sub populations
#ifdef LIFCAneuron
  neuSubPopInCM_probPerMil[RS]=800;
  neuSubPopInCM_probPerMil[FS]=200;
#else
  neuSubPopInCM_probPerMil[RS]=813;//800;//813 in the flat conf
  neuSubPopInCM_probPerMil[FS]=187;//200;//187
#endif

  //numerosity (perMil) of synapses projected by RS neurons

  //excitatoryRS to inhibitoryFS
  intraCM_RS_FS_probPerMil         = 200;
  //excitatoryRS to excitatoryRS (local or remote)
  intraCM_RS_RS_probPerMil         = 560;
  interCM_RS_1stNeighb_probPerMil  = 120/4;
  interCM_RS_2ndNeighb_probPerMil  = 80/4;
  interCM_RS_3rdNeighb_probPerMil  = 40/4;

  //numerosity of synapses projected by FS neurons
  //inhibitory only to excitatory
  intraCM_FS_RS_probPerMil         = 1000;

  if(((  intraCM_RS_FS_probPerMil + 
         intraCM_RS_RS_probPerMil +
     4 * interCM_RS_1stNeighb_probPerMil + 
     4 * interCM_RS_2ndNeighb_probPerMil + 
     4 * interCM_RS_3rdNeighb_probPerMil
                             ) != 1000)   ||
     (intraCM_FS_RS_probPerMil != 1000)) {
       printf("ERROR sum of probPerMil do not sum to 1000\n");
	   fflush(stdout);exit(0);}; 
};

void simpleCM_connectomeClass::initRandTable(
  uint32_t *pRandTable_initValue) {
  pRandTable = pRandTable_initValue;
};

void simpleCM_connectomeClass::describeConnectome(
struct DPSNN_parameters *p_lnp_par_initValue) {
  uint32_t iNeuSubPop, jSynFiber;

  p_lnp_par = p_lnp_par_initValue;

  //before assigning actual connectivity, clear the matrixes
  for(iNeuSubPop = 0; iNeuSubPop < neuSubPopTotal; iNeuSubPop++) {
    neuSubPopInCM_count[iNeuSubPop]=0;
    neuSubPopInCM_offset[iNeuSubPop]=0;
    for(jSynFiber = 0; jSynFiber < synFiberTotal;jSynFiber++) {
      synFiber_probPerMil[iNeuSubPop][jSynFiber]=0;
      synFiber_offset[iNeuSubPop][jSynFiber]=0;
      synFiber_count[iNeuSubPop][jSynFiber]=0;
    }; 
  };

  //KEY POINT change connectome constants only inside the following 
  //function  
  initProbPerMil();

  {//setting offset and counts for neural subpopulations
      neuSubPopInCM_offset[RS] = 0;
      neuSubPopInCM_count[RS] = 
        (uint32_t) 
        ((neuSubPopInCM_probPerMil[RS] * 
	  p_lnp_par->neuronsPerCM) / 1000);
      p_lnp_par->RSPerCM = neuSubPopInCM_count[RS];

      neuSubPopInCM_offset[FS] = neuSubPopInCM_count[RS];
      neuSubPopInCM_count[FS] = p_lnp_par->neuronsPerCM - 
      neuSubPopInCM_count[RS];
      #warning "should be generalized for more subpopulations"
  }

  {//reporting neural subPop probPerMil, Offset and Count 
    DPSNNverboseStart(false,1,0);
    for(iNeuSubPop = 0; iNeuSubPop < neuSubPopTotal; iNeuSubPop++) {
      printf("\n NEW NEURON SUBPOP: iNeuSubPop=%d \n",
	     iNeuSubPop); fflush(stdout); 
      printf("neuSubPopInCM_probPerMil[iNeuSubPop] = %d\n ",
	     neuSubPopInCM_probPerMil[iNeuSubPop]); 
             fflush(stdout); 
      printf("neuSubPopInCM_count[iNeuSubPop] = %d\n ",
	     neuSubPopInCM_count[iNeuSubPop]); 
             fflush(stdout); 
      printf("neuSubPopInCM_offset[iNeuSubPop] = %d\n ",
	     neuSubPopInCM_offset[iNeuSubPop]);
             fflush(stdout);
    };
    DPSNNverboseEnd();
  }
  { //defining neural count AT GLOBAL SCALE
    p_lnp_par->globNe = neuSubPopInCM_count[RS] * p_lnp_par->globCFT;
    p_lnp_par->globNi = neuSubPopInCM_count[FS] * p_lnp_par->globCFT;
    if((p_lnp_par->globNe + p_lnp_par->globNi) != p_lnp_par->globN) {
      printf(
	"ERROR simpleCM: globNe =%d +globNi=%d != globN=%d \n",
	p_lnp_par->globNe, p_lnp_par->globNi, p_lnp_par->globN);
      fflush(stdout);exit(0);
    }; 
  }; 
 
  //relative numerosity of synaptic fibers != 0
  //typical delta in cortical module connectivity

  //source neuron sub population RS
  synFiber_probPerMil[RS][intraCM_RS_FS]
                        = intraCM_RS_FS_probPerMil;
                  deltaCM[intraCM_RS_FS].cfx = 0; 
                  deltaCM[intraCM_RS_FS].cfy = 0;

  synFiber_probPerMil[RS][intraCM_RS_RS]
                        = intraCM_RS_RS_probPerMil;
                  deltaCM[intraCM_RS_RS].cfx = 0; 
                  deltaCM[intraCM_RS_RS].cfy = 0;

  {//RS -> first neighbours 
    synFiber_probPerMil[RS][interCM_xp]    
                          = interCM_RS_1stNeighb_probPerMil;
                    deltaCM[interCM_xp].cfx = 1; 
                    deltaCM[interCM_xp].cfy = 0;

    synFiber_probPerMil[RS][interCM_xm]
                          = interCM_RS_1stNeighb_probPerMil;
                    deltaCM[interCM_xm].cfx =-1; 
                    deltaCM[interCM_xm].cfy = 0;

    synFiber_probPerMil[RS][interCM_yp]
                          = interCM_RS_1stNeighb_probPerMil;
                    deltaCM[interCM_yp].cfx = 0; 
                    deltaCM[interCM_yp].cfy = 1;

    synFiber_probPerMil[RS][interCM_ym]
                          = interCM_RS_1stNeighb_probPerMil;
                    deltaCM[interCM_ym].cfx = 0; 
                    deltaCM[interCM_ym].cfy =-1;
  };
  {//RS -> second neighbours  
    synFiber_probPerMil[RS][interCM_xp_yp]
                          = interCM_RS_2ndNeighb_probPerMil;
                    deltaCM[interCM_xp_yp].cfx = 1; 
                    deltaCM[interCM_xp_yp].cfy = 1;

    synFiber_probPerMil[RS][interCM_xm_yp]
                          = interCM_RS_2ndNeighb_probPerMil;
                    deltaCM[interCM_xm_yp].cfx =-1; 
                    deltaCM[interCM_xm_yp].cfy = 1;

    synFiber_probPerMil[RS][interCM_xp_ym]
                          = interCM_RS_2ndNeighb_probPerMil;
                    deltaCM[interCM_xp_ym].cfx = 1; 
                    deltaCM[interCM_xp_ym].cfy =-1;

    synFiber_probPerMil[RS][interCM_xm_ym]
                          = interCM_RS_2ndNeighb_probPerMil;
                    deltaCM[interCM_xm_ym].cfx =-1; 
                    deltaCM[interCM_xm_ym].cfy =-1;
  };

  {//RS -> third neighbours   
    synFiber_probPerMil[RS][interCM_2xp]
                          = interCM_RS_3rdNeighb_probPerMil;
                    deltaCM[interCM_2xp].cfx = 2; 
                    deltaCM[interCM_2xp].cfy = 0;

    synFiber_probPerMil[RS][interCM_2xm]
                          = interCM_RS_3rdNeighb_probPerMil;
                    deltaCM[interCM_2xm].cfx =-2; 
                    deltaCM[interCM_2xm].cfy = 0;

    synFiber_probPerMil[RS][interCM_2yp]
                          = interCM_RS_3rdNeighb_probPerMil;
                    deltaCM[interCM_2yp].cfx = 0; 
                    deltaCM[interCM_2yp].cfy = 2;

    synFiber_probPerMil[RS][interCM_2ym]
                          = interCM_RS_3rdNeighb_probPerMil;
                    deltaCM[interCM_2ym].cfx = 0; 
                    deltaCM[interCM_2ym].cfy = -2;
  };

  //source neuron FS
  synFiber_probPerMil[FS][intraCM_FS_RS] 
                        = intraCM_FS_RS_probPerMil;
                  deltaCM[intraCM_FS_RS].cfx = 0; 
                  deltaCM[intraCM_FS_RS].cfy = 0;


  {//reporting probPerMil of synFibers
    DPSNNverboseStart(false,1,0);
      printf("report of synFiber_probPerMil starts:\n"); 
        fflush(stdout);
      printf("intraCM_FS_RS_probPerMil =%d\n",
	intraCM_FS_RS_probPerMil);fflush(stdout);
      for(iNeuSubPop = 0; iNeuSubPop < neuSubPopTotal; iNeuSubPop++) {
        for(jSynFiber = 0; jSynFiber < synFiberTotal; jSynFiber++) {
	  if(synFiber_probPerMil[iNeuSubPop][jSynFiber] != 0) {
	    printf(
	      "synFiber_probPerMil[iNeuSubPop=%d][jSynFiber=%d]=%d \n",
		 iNeuSubPop, jSynFiber, 
		 synFiber_probPerMil[iNeuSubPop][jSynFiber]); 
	      fflush(stdout);
	};
      };
    };
    DPSNNverboseEnd();
  };   
  

  { //generate synaptic counts and offsets and
    //check if the probabilities and total synapses are normalized
    
    uint32_t synFiber_totalPerMil[neuSubPopTotal];
    uint32_t synFiber_totalCount[neuSubPopTotal];

    for(iNeuSubPop = 0; iNeuSubPop < neuSubPopTotal; iNeuSubPop++) {
      synFiber_totalPerMil[iNeuSubPop]=0;
      synFiber_totalCount[iNeuSubPop]=0;
    };


    for(iNeuSubPop = 0; iNeuSubPop < neuSubPopTotal; iNeuSubPop++) {
      for(jSynFiber = 0; jSynFiber < synFiberTotal;jSynFiber++) {   
	synFiber_count[iNeuSubPop][jSynFiber] =
	  (uint32_t) (p_lnp_par->M *  
	   synFiber_probPerMil[iNeuSubPop][jSynFiber])/1000;
	if(jSynFiber > 0) {
	  synFiber_offset[iNeuSubPop][jSynFiber]=
	    synFiber_offset[iNeuSubPop][jSynFiber-1]+
	    synFiber_count[iNeuSubPop][jSynFiber-1];	    
	};
	synFiber_totalCount[iNeuSubPop] +=
	  synFiber_count[iNeuSubPop][jSynFiber];
	synFiber_totalPerMil[iNeuSubPop] += 
	  synFiber_probPerMil[iNeuSubPop][jSynFiber];
      };
      if(synFiber_totalPerMil[iNeuSubPop] != 1000) {
	printf(
        "ERROR synFiber_probPerMil[iNeuSubPop=%d]=%d should be 1000\n",
          iNeuSubPop,synFiber_totalPerMil[iNeuSubPop]); 
	  fflush(stdout);exit(0);
      };
      if(synFiber_totalCount[iNeuSubPop] != p_lnp_par->M) {
	printf("ERROR synFiber_totalCount not normalized to M \n"); 
	fflush(stdout);exit(0);
      };
    };
  };

  

  {//report counts and offsets of synFibers
    DPSNNverboseStart(false,1,0);
    for(iNeuSubPop = 0; iNeuSubPop < neuSubPopTotal; iNeuSubPop++) {
      for(jSynFiber = 0; jSynFiber < synFiberTotal; jSynFiber++) {
	if(synFiber_count[iNeuSubPop][jSynFiber] != 0) {
	  printf("\n NEW synFiber KIND: iNeuSubPop=%d jSynFiber=%d \n",
	       iNeuSubPop, jSynFiber); fflush(stdout);
	  printf("synFiber_count[iNeuSubPop][jSynFiber] = %d \n",
	       synFiber_count[iNeuSubPop][jSynFiber]); 
	       fflush(stdout);
	  printf("synFiber_offset[iNeuSubPop][jSynFiber] = %d \n",
	       synFiber_offset[iNeuSubPop][jSynFiber]);
	       fflush(stdout); 
	};
      };
    }; 
    DPSNNverboseEnd();
  };   
};

bool simpleCM_connectomeClass::isInterCMSynFiber(
  synFiberEnum synFiberKind) {
  if((synFiberKind >= interCM_xp) &&
     (synFiberKind <= interCM_2ym))
       return(true);
  else
       return(false);
};

bool simpleCM_connectomeClass::isIntraCMSynFiber(
  synFiberEnum synFiberKind){
  if((synFiberKind >= intraCM_RS_RS) &&
     (synFiberKind <= intraCM_FS_RS)) 
       return(true);
  else
       return(false);
};

simpleCM_neuCoordinatesStruct 
simpleCM_connectomeClass::convert_loc_n_h_to_neuCMCoordinates(
  uint32_t loc_n, uint32_t loc_h)
{ //compute coordinates for this neuron i in
      //other coordinate systems: global system and cortical modules system
  simpleCM_neuCoordinatesStruct neuCMCoordinates;
      neuCMCoordinates.loc_n = loc_n;
      neuCMCoordinates.loc_h = loc_h;
      neuCMCoordinates.glob_n = 
	neuCMCoordinates.loc_n + p_lnp_par->locN * loc_h;
      neuCMCoordinates.cfy_n = 
        (neuCMCoordinates.glob_n / p_lnp_par->neuronsPerCM) / 
	p_lnp_par->globCFX;
      neuCMCoordinates.cfx_n = 
        (neuCMCoordinates.glob_n / p_lnp_par->neuronsPerCM) % 
	p_lnp_par->globCFX;
      neuCMCoordinates.inCM_n = 
	neuCMCoordinates.glob_n % p_lnp_par->neuronsPerCM;
      neuCMCoordinates.subPop = 
	getNeuralSubPop(neuCMCoordinates.inCM_n);
      neuCMCoordinates.neuralKind = 
        getNeuralKind(neuCMCoordinates.inCM_n);
      return(neuCMCoordinates);
};

simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::generateTargetNeu(
    simpleCM_neuCoordinatesStruct sourceNeuCMCoordinates, 
    uint32_t jSynIdInNeu)
{
  simpleCM_neuCoordinatesStruct targetNeu;
  synFiberEnum synFiberKind;
  if(jSynIdInNeu > (p_lnp_par->M - 1)) {
      printf("ERROR: invalid jSynIdInNeu=%d in generateTargetNeu\n",
	     jSynIdInNeu);
      fflush(stdout);
      exit(0);
  }
  
  switch(sourceNeuCMCoordinates.subPop) {
  case RS:
    synFiberKind = getSynFiberType(jSynIdInNeu,RS);
    if(isIntraCMSynFiber(synFiberKind))
      targetNeu = intraModule_targetNeu(
         jSynIdInNeu, p_lnp_par->synGen, 
	 sourceNeuCMCoordinates);
    else
    if(isInterCMSynFiber(synFiberKind))
      targetNeu = interModule_targetNeu(
	jSynIdInNeu, synFiberKind, p_lnp_par->synGen, 
	sourceNeuCMCoordinates);
    else
      {printf("ERROR: non implem. in generateTargetNeu-A\n");
	fflush(stdout);exit(0);}
    break;
  case FS: 
    synFiberKind = getSynFiberType(jSynIdInNeu,FS);
    switch(synFiberKind) {
    case intraCM_FS_RS:
      targetNeu = intraModule_targetNeu(
         jSynIdInNeu, p_lnp_par->synGen, 
	 sourceNeuCMCoordinates);
      break;
    default:
      printf("ERROR: non implem. in generateTargetNeu-B\n");
      fflush(stdout);exit(0);
      break;
    }
    break;
  default:
    printf("ERROR: non implem. in generateTargetNeu-C\n");
    fflush(stdout);exit(0);
    break;
  };
  return(targetNeu);
};

synFiberEnum simpleCM_connectomeClass::getSynFiberType(
  uint32_t synIdInSourceNeu, neuSubPopEnum sourceNeuSubPop)
{

  synFiberEnum synFiberKind;
  switch(sourceNeuSubPop) {
  case RS: 
    synFiberKind = intraCM_RS_RS;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = intraCM_RS_FS;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_xp;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_xm;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);  

    synFiberKind = interCM_yp;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_ym;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_xp_yp;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_xp_ym;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_xm_yp;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_xm_ym;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_2xp;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_2xm;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_2yp;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);

    synFiberKind = interCM_2ym;
    if(isThisSynFiberKind(synIdInSourceNeu,RS,synFiberKind))
      return(synFiberKind);
    { printf("ERROR should never get here in getSynFiber-A\n");
	printf("synIdInSourceNeu=%d sourceNeuSubPop=%d\n",
	       synIdInSourceNeu,sourceNeuSubPop); 
	fflush(stdout);exit(0);}
    break;
  case FS: 
    synFiberKind = intraCM_FS_RS;
    if(isThisSynFiberKind(synIdInSourceNeu,FS,intraCM_FS_RS))
      return(synFiberKind);
    else
      {printf("ERROR should never get here in getSynFiber-B\n");
       fflush(stdout);exit(0);};
    break;
  default:
    {printf("ERROR should never get here in getSynFiber-C\n");
     fflush(stdout);exit(0);};
    break;
  };
}; 

bool simpleCM_connectomeClass::isThisSynFiberKind(
  uint32_t jSynIdInNeu, neuSubPopEnum neuSubPopInCM, 
  synFiberEnum synFiberKind)
{  
  if((jSynIdInNeu >= synFiber_offset[neuSubPopInCM][synFiberKind]) && 
     (jSynIdInNeu < (synFiber_offset[neuSubPopInCM][synFiberKind] +
		     synFiber_count[neuSubPopInCM][synFiberKind])) &&
     (synFiber_count[neuSubPopInCM][synFiberKind] != 0)) 
    return(true);
  else
    return(false);
};

neuSubPopEnum simpleCM_connectomeClass::getNeuralSubPop(uint32_t iNeuInCM) 
{
  //  if(neuSubPopInCM_count[RS] != 0) // Check not needed now. To be reintroduced when more subPop available
    if((iNeuInCM >= neuSubPopInCM_offset[RS])  && 
       (iNeuInCM < (neuSubPopInCM_count[RS] + 
		    neuSubPopInCM_offset[RS]) ) )
       return(RS);
  else 
    return(FS);

#warning "getNeuralSubPop should be generalized"

};

neuralKindEnum simpleCM_connectomeClass::getNeuralKind(uint32_t iNeuInCM) 
{
  neuSubPopEnum subPopInCM;
  neuralKindEnum neuralKind;
  subPopInCM = getNeuralSubPop(iNeuInCM);
  switch(subPopInCM) {
  case RS:
    neuralKind = excitatoryRS;
    break;
  case FS:
    neuralKind = inhibitoryFS;
    break;
  default:
    printf("ERROR: unrecogn. subPop->neural kind\n");
    fflush(stdout);exit(0);
    break;
  };
return(neuralKind);
};

uint32_t simpleCM_connectomeClass::conv_neuIdInCM_to_glob_n(
  uint32_t iNeuIdInCM, uint32_t icfx, uint32_t icfy) {
  uint32_t neuGlobId;
  if((iNeuIdInCM >= p_lnp_par->neuronsPerCM) ||
     (icfx >= p_lnp_par->globCFX)            ||
     (icfy >= p_lnp_par->globCFY))
    {printf(
     "ERROR in conv_neuIdInCM_to_glob_n, iNeuIdInCM=%d,icfx=%d,icfy=%d\n",  
     iNeuIdInCM,icfx,icfy);fflush(stdout);exit(0);};

  neuGlobId = iNeuIdInCM + 
      icfx * p_lnp_par->neuronsPerCM +
      icfy * p_lnp_par->neuronsPerCM *p_lnp_par->globCFX;

  if(neuGlobId>=p_lnp_par->globN) {printf(
   "ERROR conv_neuIdInCM_to_glob_n, would ret. neuGlobId=%d out of range\n",  
     neuGlobId);fflush(stdout);exit(0);};

  return(neuGlobId);
};

simpleCM_neuCoordinatesStruct
  simpleCM_connectomeClass::intraModule_targetNeu(
  uint32_t mId, 
  synGenEnum synGenKind, 
  simpleCM_neuCoordinatesStruct sourceNeu) 
{
  simpleCM_neuCoordinatesStruct targetNeu;
  
  switch(synGenKind) {
  case default_random_synGen_1: 
      targetNeu = 
	intraModule_random_targetNeu(mId, sourceNeu);
    break;
  case simpleCorticalModule_synGen_2: 
      targetNeu = 
	intraModule_nonRand_targetNeu(mId, sourceNeu);
    break;
  case randTable_simpleCorticalModule_synGen_3: 
      targetNeu = 
	intraModule_randTable_targetNeu(mId, sourceNeu);
    break;
  default:
    printf("ERROR: non implemented synGen\n");fflush(stdout);exit(0);
    break;
  };
  return(targetNeu);
};

simpleCM_neuCoordinatesStruct
  simpleCM_connectomeClass::interModule_targetNeu(
    uint32_t mId,
    synFiberEnum synFiberKind,
    synGenEnum synGenKind,
    simpleCM_neuCoordinatesStruct sourceNeu)
{
  simpleCM_neuCoordinatesStruct targetNeu;
  switch(synGenKind) {
  case default_random_synGen_1: 
      targetNeu = 
	interModule_random_targetNeu(mId, synFiberKind, sourceNeu);
    break;
  case simpleCorticalModule_synGen_2: 
      targetNeu = 
	interModule_nonRand_targetNeu(mId, synFiberKind, sourceNeu);
    break;
  case randTable_simpleCorticalModule_synGen_3: 
      targetNeu = 
	interModule_randTable_targetNeu(mId, synFiberKind, sourceNeu);
    break;
  default:
    printf("ERROR: non implemented synGen\n");fflush(stdout);exit(0);
    break;
  };
  return(targetNeu);
};

simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::intraModule_random_targetNeu(
  uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu) 
{
  simpleCM_neuCoordinatesStruct targetNeu;
  uint32_t k, r, exists;

  targetNeu=sourceNeu;
  
  do { 
    //search for a valid post-synaptic neuron index r
    if(countRandSynGen++ > (p_lnp_par->M*8)) {
      printf("ERROR: in gen.For.Conn. more than M*8 random attempts on neu=%d syn=%d \n",sourceNeu.glob_n,mId); 
      fflush(stdout);exit(0);
    };
    exists = 0;// to avoid multiple synapses
    switch(sourceNeu.subPop) {
    case RS: //RS projects to both FS and RS
      //r = (sourceNeu.inCM_n + getRandom(p_lnp_par->neuronsPerCM))%p_lnp_par->neuronsPerCM;
      r = (sourceNeu.inCM_n + getRandom_r(&seedForSynapses,p_lnp_par->neuronsPerCM))%p_lnp_par->neuronsPerCM;
      break;
    case FS: 
#ifdef LIFCAneuron
      //r = (sourceNeu.inCM_n + getRandom(p_lnp_par->neuronsPerCM))%p_lnp_par->neuronsPerCM;
      r = (sourceNeu.inCM_n + getRandom_r(&seedForSynapses,p_lnp_par->neuronsPerCM))%p_lnp_par->neuronsPerCM;
#else //Izhikevich FS projects to RS!!
      //r = (getRandom(neuSubPopInCM_count[RS]));
      r = getRandom_r(&seedForSynapses,neuSubPopInCM_count[RS]);
#endif
      // inh -> exc only
      break;
    default:
      printf("ERROR: unknown Neural Kind in synGen\n");fflush(stdout); exit(0);
      break;
    };
    if (sourceNeu.glob_n==r) { // no self-synapses allowed
      exists=1;
      DPSNNverboseStart(false,0,1);
      printf("gen.Forw.Conn. R-h%d REJECTED source glob_n %d == target %d\n",
	     p_lnp_par->loc_h, sourceNeu.glob_n, r);
      DPSNNverboseEnd();
    }; 
    for (k=0;k<mId;k++) {
      if (synListOfThisNeu[k]==r) {
	//no duplicated synapses with same target allowed
	exists = 1; //already existing found
	DPSNNverboseStart(false,0,1);
	printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
	       p_lnp_par->loc_h, sourceNeu.glob_n, r);
	DPSNNverboseEnd();
      }; // synapse already exists
    };  
  } while (exists == 1); // if exist==1 try using another random r
  
  synListOfThisNeu[mId] = r;
  targetNeu.inCM_n = r;
  targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n,targetNeu.cfx_n,targetNeu.cfy_n);
  targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
  return(targetNeu);
};
 
simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::interModule_random_targetNeu(
    uint32_t mId,  
    synFiberEnum synFiberKind, 
    simpleCM_neuCoordinatesStruct sourceNeu) 
{
  int32_t naive_icfx_target, naive_icfy_target;
  simpleCM_neuCoordinatesStruct targetNeu;
  uint32_t k, r, exists;
 
  naive_icfx_target = deltaCM[synFiberKind].cfx + sourceNeu.cfx_n;
  naive_icfy_target = deltaCM[synFiberKind].cfy + sourceNeu.cfy_n;
    
  //naive, because there should be exceptions at the boundaries
  if((naive_icfx_target >= 0) &&
     ((uint32_t)naive_icfx_target < p_lnp_par->globCFX))
    targetNeu.cfx_n = (uint32_t)naive_icfx_target;
  else 
#warning " periodic boundary conditions on CFX"
    targetNeu.cfx_n = (uint32_t)(naive_icfx_target % p_lnp_par->globCFX);

  //naive, because there should be exceptions at the boundaries
  if((naive_icfy_target >= 0) &&
     ((uint32_t)naive_icfy_target < p_lnp_par->globCFY-1))
    targetNeu.cfy_n = (uint32_t) naive_icfy_target;
  else 
#warning " periodic boundary conditions on CFY"
    targetNeu.cfy_n = (uint32_t)(naive_icfy_target % p_lnp_par->globCFY);
     
  //generation of target neural id in target CM
  do { 
    //search for a valid post-synaptic neuron index r
    if(countRandSynGen++ > (p_lnp_par->M*8)) {
      printf("ERROR: in gen.For.Conn. more than M*8 random attempts on neu=%d syn=%d \n",sourceNeu.glob_n,mId); 
      fflush(stdout);exit(0);
    };
    exists = 0;// to avoid multiple synapses

    switch(sourceNeu.subPop) {
    case RS:  
      //r = (sourceNeu.inCM_n + getRandom(p_lnp_par->neuronsPerCM))%p_lnp_par->neuronsPerCM;
      r = (sourceNeu.inCM_n + getRandom_r(&seedForSynapses,p_lnp_par->neuronsPerCM))%p_lnp_par->neuronsPerCM;
      break;
    case FS: 
      printf("ERROR: wrong sourceNeuKind in interMod synGen\n");
      fflush(stdout); exit(0);
      break;
    default:
      printf("ERROR: unknown Neural Kind in synGen\n");fflush(stdout); exit(0);
      break;
    };
    if (sourceNeu.glob_n==r) { // no self-synapses allowed
      exists=1;
      DPSNNverboseStart(false,0,1);
      printf("gen.Forw.Conn. R-h%d REJECTED source glob_n %d == target %d\n",
	     p_lnp_par->loc_h, sourceNeu.glob_n, r);
      DPSNNverboseEnd();
    }; 
    for (k=0;k<mId;k++) {
      if (synListOfThisNeu[k]==r) {
	//no duplicated synapses with same target allowed
	exists = 1; //already existing found
	DPSNNverboseStart(false,0,1);
	printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
	       p_lnp_par->loc_h, sourceNeu.glob_n, r);
	DPSNNverboseEnd();
      }; // synapse already exists
    };  
  } while (exists == 1); // if exist==1 try using another random r

  synListOfThisNeu[mId] = r;
  targetNeu.inCM_n = r;
  targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n, targetNeu.cfx_n, targetNeu.cfy_n);
  targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
  return(targetNeu);
};

simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::intraModule_nonRand_targetNeu(
  uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu) 
{
  simpleCM_neuCoordinatesStruct targetNeu;
  unsigned int detOffset,detSpan,detSum;
  //generation of target id
  targetNeu=sourceNeu;
  switch(sourceNeu.subPop) {
  case RS: //RS projects to both FS and RS
    detSpan = (uint32_t) 
      (mId * (1 + (uint32_t) ((float) (p_lnp_par->neuronsPerCM/8) / 
		              (float) p_lnp_par->M) )); 
    targetNeu.inCM_n = 
    (sourceNeu.inCM_n + 1 + detSpan) % p_lnp_par->neuronsPerCM;
    /*
    if((mId%2)==0) {
      targetNeu.inCM_n = 
      (sourceNeu.inCM_n + 1 + detSpan) % p_lnp_par->neuronsPerCM;
    }else{
      targetNeu.inCM_n = 
      (sourceNeu.inCM_n - detSpan) % p_lnp_par->neuronsPerCM;
    }
    */
    break;
  case FS: //FS projects to RS!!
    
    detOffset = sourceNeu.inCM_n - neuSubPopInCM_offset[FS];

    detSpan = (unsigned int) 
	(mId * ((float) (neuSubPopInCM_count[RS]) / 
		(float) p_lnp_par->M));
    detSum = detOffset + detSpan;
    targetNeu.inCM_n = detSum % neuSubPopInCM_count[RS]; 
    // inh -> exc only
    break;
  default:
    printf("ERROR: unknown Neural Kind in synGen\n");fflush(stdout); exit(0);
    break;
  };
  targetNeu.glob_n = 
     conv_neuIdInCM_to_glob_n(
       targetNeu.inCM_n,targetNeu.cfx_n,targetNeu.cfy_n);
  targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
  return(targetNeu);
};
 
simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::interModule_nonRand_targetNeu(
    uint32_t mId,  
    synFiberEnum synFiberKind, 
    simpleCM_neuCoordinatesStruct sourceNeu) 
{


  int32_t naive_icfx_target, naive_icfy_target;
  uint32_t detSpan;
  simpleCM_neuCoordinatesStruct targetNeu;

  naive_icfx_target = deltaCM[synFiberKind].cfx + sourceNeu.cfx_n;
  naive_icfy_target = deltaCM[synFiberKind].cfy + sourceNeu.cfy_n;
    
  //naive, because there should be exceptions at the boundaries
  if((naive_icfx_target >= 0) &&
     ((uint32_t)naive_icfx_target < p_lnp_par->globCFX))
    targetNeu.cfx_n = (uint32_t)naive_icfx_target;
  else 
    #warning " periodic boundary conditions on CFX"
    targetNeu.cfx_n = (uint32_t)(naive_icfx_target % p_lnp_par->globCFX);

  //naive, because there should be exceptions at the boundaries
  if((naive_icfy_target >= 0) &&
     ((uint32_t)naive_icfy_target < p_lnp_par->globCFY-1))
    targetNeu.cfy_n = (uint32_t) naive_icfy_target;
  else 
    #warning " periodic boundary conditions on CFY"
    targetNeu.cfy_n = (uint32_t)(naive_icfy_target % p_lnp_par->globCFY);
     
  //generation of target neural id in target CM
  switch(sourceNeu.subPop) {
  case RS:  
    detSpan = (uint32_t) 
      (mId * (1 + (uint32_t) ((float) (p_lnp_par->neuronsPerCM/8) / 
		              (float) p_lnp_par->M) )); 
    targetNeu.inCM_n = 
      (sourceNeu.inCM_n + 1 + detSpan) % p_lnp_par->neuronsPerCM;
    break;
  case FS: 
    printf("ERROR: wrong sourceNeuKind in interMod synGen\n");
    fflush(stdout); exit(0);
    break;
  default:
    printf("ERROR: unknown Neural Kind in synGen\n");fflush(stdout); exit(0);
    break;
  };
  targetNeu.glob_n = 
     conv_neuIdInCM_to_glob_n(
       targetNeu.inCM_n, targetNeu.cfx_n, targetNeu.cfy_n);
  targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
  return(targetNeu);
};

simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::intraModule_randTable_targetNeu(
  uint32_t mId, simpleCM_neuCoordinatesStruct sourceNeu) 
{
  simpleCM_neuCoordinatesStruct targetNeu;
  uint32_t k, r, exists;

  targetNeu=sourceNeu;
  
  do { 
    //search for a valid post-synaptic neuron index r
    if(countRandSynGen++ > (p_lnp_par->M*4)) {
      printf("ERROR: in gen.For.Conn. more than M*4 random attempts on neu=%d syn=%d \n",sourceNeu.glob_n,mId); 
      fflush(stdout);exit(0);
    };
    exists = 0;// to avoid multiple synapses
    switch(sourceNeu.subPop) {
    case RS: //RS projects to both FS and RS
      r = (sourceNeu.inCM_n + pRandTable[(countRandSynGen+ sourceNeu.inCM_n + 7*mId)%DSD__randTable])%p_lnp_par->neuronsPerCM;
      break;
    case FS: //FS projects to RS!!
      r = (pRandTable[(countRandSynGen + sourceNeu.inCM_n + 7*mId)%DSD__randTable])%neuSubPopInCM_count[RS];
      // inh -> exc only
      break;
    default:
      printf("ERROR: unknown Neural Kind in synGen\n");fflush(stdout); exit(0);
      break;
    };
    if (sourceNeu.glob_n==r) { // no self-synapses allowed
      exists=1;
      DPSNNverboseStart(false,0,1);
      printf("gen.Forw.Conn. R-h%d REJECTED source glob_n %d == target %d\n",
	     p_lnp_par->loc_h, sourceNeu.glob_n, r);
      DPSNNverboseEnd();
    }; 
    for (k=0;k<mId;k++) {
      if (synListOfThisNeu[k]==r) {
	//no duplicated synapses with same target allowed
	exists = 1; //already existing found
	DPSNNverboseStart(false,0,1);
	printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
	       p_lnp_par->loc_h, sourceNeu.glob_n, r);
	DPSNNverboseEnd();
      }; // synapse already exists
    };  
  } while (exists == 1); // if exist==1 try using another random r
  
  synListOfThisNeu[mId] = r;
  targetNeu.inCM_n = r;
  targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n,targetNeu.cfx_n,targetNeu.cfy_n);
  targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
  return(targetNeu);
};
 
simpleCM_neuCoordinatesStruct 
  simpleCM_connectomeClass::interModule_randTable_targetNeu(
    uint32_t mId,  
    synFiberEnum synFiberKind, 
    simpleCM_neuCoordinatesStruct sourceNeu) 
{
  int32_t naive_icfx_target, naive_icfy_target;
  simpleCM_neuCoordinatesStruct targetNeu;
  uint32_t k, r, exists;
 
  naive_icfx_target = deltaCM[synFiberKind].cfx + sourceNeu.cfx_n;
  naive_icfy_target = deltaCM[synFiberKind].cfy + sourceNeu.cfy_n;
    
  //naive, because there should be exceptions at the boundaries
  if((naive_icfx_target >= 0) &&
     ((uint32_t)naive_icfx_target < p_lnp_par->globCFX))
    targetNeu.cfx_n = (uint32_t)naive_icfx_target;
  else 
#warning " periodic boundary conditions on CFX"
    targetNeu.cfx_n = (uint32_t)(naive_icfx_target % p_lnp_par->globCFX);

  //naive, because there should be exceptions at the boundaries
  if((naive_icfy_target >= 0) &&
     ((uint32_t)naive_icfy_target < p_lnp_par->globCFY-1))
    targetNeu.cfy_n = (uint32_t) naive_icfy_target;
  else 
#warning " periodic boundary conditions on CFY"
    targetNeu.cfy_n = (uint32_t)(naive_icfy_target % p_lnp_par->globCFY);
     
  //generation of target neural id in target CM
  do { 
    //search for a valid post-synaptic neuron index r
    if(countRandSynGen++ > (p_lnp_par->M*4)) {
      printf("ERROR: in gen.For.Conn. more than M*4 random attempts on neu=%d syn=%d \n",sourceNeu.glob_n,mId); 
      fflush(stdout);exit(0);
    };
    exists = 0;// to avoid multiple synapses

    switch(sourceNeu.subPop) {
    case RS:  
      r = (sourceNeu.inCM_n + pRandTable[(countRandSynGen + sourceNeu.inCM_n + 7*mId)%DSD__randTable])%p_lnp_par->neuronsPerCM;
      break;
    case FS: 
      printf("ERROR: wrong sourceNeuKind in interMod synGen\n");
      fflush(stdout); exit(0);
      break;
    default:
      printf("ERROR: unknown Neural Kind in synGen\n");fflush(stdout); exit(0);
      break;
    };
    if (sourceNeu.glob_n==r) { // no self-synapses allowed
      exists=1;
      DPSNNverboseStart(false,0,1);
      printf("gen.Forw.Conn. R-h%d REJECTED source glob_n %d == target %d\n",
	     p_lnp_par->loc_h, sourceNeu.glob_n, r);
      DPSNNverboseEnd();
    }; 
    for (k=0;k<mId;k++) {
      if (synListOfThisNeu[k]==r) {
	//no duplicated synapses with same target allowed
	exists = 1; //already existing found
	DPSNNverboseStart(false,0,1);
	printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
	       p_lnp_par->loc_h, sourceNeu.glob_n, r);
	DPSNNverboseEnd();
      }; // synapse already exists
    };  
  } while (exists == 1); // if exist==1 try using another random r

  synListOfThisNeu[mId] = r;
  targetNeu.inCM_n = r;
  targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n, targetNeu.cfx_n, targetNeu.cfy_n);
  targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
  return(targetNeu);
};

void simpleCM_connectomeClass::report() {
};

