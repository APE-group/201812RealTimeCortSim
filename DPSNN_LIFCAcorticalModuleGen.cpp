// DPSNN_LIFCAcorticalModuleGen.cpp
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
#include "DPSNN_LIFCAconnectome.h"
#ifdef LIFCAneuron
#include "randdev.h"
#endif

void localNetClass::simpleCM_prepareForwardSynapses() {
  chronoClass chronoDescribeConnectome;  
  chronoClass chronoInvokeIndivSynGen;
  DPSNNverboseStart(true,1,0);
    chronoDescribeConnectome.clearAndStartChrono();
  DPSNNverboseEnd();

  simpleCM_connectome.describeConnectome(&lnp_par,neuSubPopParam);

  DPSNNverboseStart(true,1,0);
    chronoDescribeConnectome.stopChrono(); 
    if(lnp_par.loc_h <= 1 ||  lnp_par.loc_h>=(lnp_par.globH-2)) {
      printf("CHRONO:describeConnectome h=%d = %f sec \n",
        lnp_par.loc_h, chronoDescribeConnectome.getAccumulatedChrono());
      fflush(stdout);};
  DPSNNverboseEnd();

  initLUT();

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
  uint32_t totSynNumPerThisNeu;
  uint64_t totSynNum;
  simpleCM_neuCoordinatesStruct sourceNeuIds;
  float thisNeuNuExt;
  float thisNeuCExt;
  uint32_t thisSubPopTotN;
  neuSubPopEnum sourceSubPop;
  neuSubPopEnum targetSubPop;
  uint32_t thisTab;

  totSynNumPerThisNeu = 0;
  totSynNum = 0;

  for(i=0;i < lnp_par.locN; i++) {
    uint32_t localSeed;
    sourceNeuIds = simpleCM_connectome.convert_loc_n_h_to_neuCMCoordinates(i,lnp_par.loc_h, &lnp_par);
    sourceSubPop = sourceNeuIds.subPop;
    n[i].initGlobalAwareness(lnp_par);
    n[i].initStat(pStat);
    n[i].set_glob_n(sourceNeuIds.glob_n);
    n[i].set_loc_n(i);
    n[i].set_loc_h(lnp_par.loc_h);
    //n[i].initM(lnp_par.M);
    //n[i].initD(lnp_par.D);
    n[i].clearForwardConnections();
    n[i].initNeuralKind(sourceNeuIds.neuralKind,sourceNeuIds.subPop,neuSubPopParam); 
    //n[i].setTableLUT_synBath(tableLUT_synBath_Exc,tableLUT_synBath_Inh);
    n[i].setTableLUT_synBathMatrix((uint64_t)&(tableLUT_synBathMatrix[sourceSubPop][0]));
    n[i].didItFire = false;
      
    if ((sourceNeuIds.glob_n % lnp_par.tabNeuron) == 0){
    thisTab = (uint32_t)(sourceNeuIds.glob_n / lnp_par.tabNeuron);
    //localSeed = sourceNeuIds.glob_n + lnp_par.globN * lnp_par.globalSeed;
    //localSeed = ((0xAFDE75 + thisTab * DSD__maxTabNumber) % DSD__INT_MAX +
    //		 lnp_par.globN * lnp_par.globalSeed) % DSD__INT_MAX;
    localSeed = thisTab;//getRandom(DSD__INT_MAX);
    if(localSeed >= DSD__INT_MAX) {
      printf("ERROR in CorticalModGen: localSeed out of range, > %d\n",DSD__INT_MAX); fflush(stdout); exit(0);
    }
    localNetRandDev.SetRandomSeed(localSeed);
    //localNetRandDev.SetRandomSeed(uint32_t(time(NULL)));
    }
      
    simpleCM_connectome.countRandSynGen = 0;
    
    totSynNumPerThisNeu=simpleCM_connectome.generateTargetNeuList(sourceNeuIds,simpleCM_connectome.targetNeuList, &lnp_par);

    if(totSynNumPerThisNeu >= DSD__maxM){
      printf(
	     "Error: attempting to generate %d synapses larger then DSD__maxM=%d\n",totSynNumPerThisNeu,DSD__maxM);
      fflush(stdout);exit(0);
    }
    
    for(jSynIdInNeu=0;jSynIdInNeu < totSynNumPerThisNeu;jSynIdInNeu++) 
      {
	synapseClass generatedSyn;
	uint32_t target_h;
	target_h =  simpleCM_connectome.targetNeuList[jSynIdInNeu].loc_h;

	generatedSyn.pre_glob_n = sourceNeuIds.glob_n;
	generatedSyn.post_glob_n = 
	  simpleCM_connectome.targetNeuList[jSynIdInNeu].glob_n;

	n[i].forwardNeuralTargetDistrInHost[target_h] ++;

	targetSubPop = simpleCM_connectome.targetNeuList[jSynIdInNeu].subPop;

        #ifdef DelayUniform
	   generatedSyn.delay = getRandomDelay_UNI(neuSubPopParam,sourceSubPop,targetSubPop);
        #else
	   generatedSyn.delay = getRandomDelay_EXP(neuSubPopParam,sourceSubPop,targetSubPop);
        #endif

	generatedSyn.preSynNeuralKind = (uint8_t) sourceNeuIds.subPop;
	//EPA this could be separated from the previous one because it seems to go faster
	//still to be verified for large configurations

	//generatedSyn.weight=
	//  tableLUT_synMatrix[sourceSubPop][targetSubPop]
	//                    [(uint32_t)(localNetRandDev.Random()*ANALOG_DEPTH)];
	
	// The following define allow the usage of the Perseo-like synaptic weights generation
	// based on discrete values stored in a LUT, or, if not defined, in "continuous" mode
	// with a gaussian distribution as in NEST
	//#define synWeightWithTableLUT
	{
	  float weightF;
	  #ifdef synWeightWithTableLUT
	  weightF =
	    tableLUT_synMatrix[sourceSubPop][targetSubPop]
	                      [(uint32_t)(localNetRandDev.Random()*ANALOG_DEPTH)];
	  #else
	  float mean;
	  float sigma;
	  mean = neuSubPopParam[sourceSubPop].J[targetSubPop];
	  sigma = neuSubPopParam[sourceSubPop].DJ[targetSubPop]*mean;

	  if(mean > 0.0)
	    do{
	      weightF = mean + localNetRandDev.NormDev() * sigma;
	    } while ( ( weightF < 0.0 ) or ( weightF >= 2.0*mean ) );
	  else
	    do{
	      weightF = mean + localNetRandDev.NormDev() * sigma;
	    } while ( ( weightF <= 2.0*mean ) or ( weightF > 0 ) );
	  #endif

	  DPSNNverboseStart(true,1,0);
	  if((weightF<lnp_par.minSynWeight_f)||(weightF>(-lnp_par.minSynWeight_f)))
	    {
	    printf(
	    "Error: requested syn weight %f out of minSynWeight_f range %f\n",
	    weightF,lnp_par.minSynWeight_f);
	    fflush(stdout);exit(0);
	    }
	  DPSNNverboseEnd();
	  //conversion factor from float to int16
	  generatedSyn.weight = (int16_t)(weightF * lnp_par.factorFloat_2_weightType);
	}
	

	if(target_h == lnp_par.loc_h) {
	  localSynList[localSynCount] = generatedSyn;
	  localSynCount++;
	  
	} else{
	  forwardSynList[forwardSynCount] = generatedSyn;
	  forwardSynCount++;
	} 
      };


    n[i].initM(totSynNumPerThisNeu);
    totSynNum += totSynNumPerThisNeu;
  
    thisNeuNuExt = neuSubPopParam[sourceSubPop].NuExt;
    thisNeuCExt = simpleCM_connectome.generateSourceNeuCExt(sourceNeuIds.cfx_n,sourceNeuIds.cfy_n, &lnp_par);
    thisSubPopTotN = neuSubPopParam[sourceSubPop].count;
    //if(thisNeuCExt==0) {printf("ERROR in NeuCExt calculation for neu=%d \n",i);fflush(stdout);exit(0);}
    if(thisNeuCExt!=0){
      n[i].invNuExt = 1000.0 / (thisNeuNuExt * thisNeuCExt /* * thisSubPopTotN */);
      n[i].maxExternCurrentCount = (uint32_t)ceil(20 / n[i].invNuExt);
    } else {
      n[i].invNuExt = 0;
      n[i].maxExternCurrentCount = 0;
    }
  };//END loop on all source neurons 

  projectedSynCount = localSynCount + forwardSynCount;
  DPSNNverboseStart(false,1,0);
  printf("Count of generated synapses on h=%d: projectedSynCount=%d of which localSynCount=%d and forwardSynCount=%d \n",lnp_par.loc_h,projectedSynCount,localSynCount,forwardSynCount);
      fflush(stdout);
  DPSNNverboseEnd();

};

/* The following code has been imported from Perseo */

/*--------------------------------------------*
 *                                            *
 *   EXP (Exponential) delays distribution.   *
 *                                            *
 *--------------------------------------------*/


/** 
 *  Returns the number of layers from an 
 *  EXPONENTIAL delay distribution starting from
 *  DMin and cut at DMax in order to neglect a 
 *  portion TailNeglect to the exponential p.d.f.
 *  The mean transmission delay is given by
 *
 *   DMax TD - DMin     DMin - DMax
 *  ---------------- + -------------
 *       TD - 1           Log[TD]
 *
 *  where TD = TailNeglected.
 */

#define TailNeglected 0.05 /* p.d.f. tail of the exponential distribution to neglect. */

uint32_t localNetClass::getRandomDelay_EXP(neuSubPopParamStruct *neuSubPopParam,neuSubPopEnum sourceSubPop,neuSubPopEnum targetSubPop)
{
   float InvLogTN;
   float outfuncF;
   float DMin,DMax;
   float DelayMin,DelayMax;
   float DelayStep,DelayNum;
   uint32_t outfunc;

   DMin = neuSubPopParam[sourceSubPop].DMin[targetSubPop];
   DMax = neuSubPopParam[sourceSubPop].DMax[targetSubPop];
   DelayMin = lnp_par.delayMin;
   DelayMax = lnp_par.delayMax;
   DelayNum = DSD__delayNum;
   InvLogTN = 1.0 / log(TailNeglected);

   DelayStep = (DelayMax - DelayMin)/DelayNum;
   DelayMax -= DelayStep/2;
   DelayMin += DelayStep/2;
   DMax = pErfLib->roundr2i((DMax - DelayMin)/DelayStep)*DelayStep+DelayMin;
   DMin = pErfLib->roundr2i((DMin - DelayMin)/DelayStep)*DelayStep+DelayMin;
   // The following code used in Perseo can't be used
   // otherwise the higher delay will be DMax-2 !!
   //if (DMax > DelayMax) DMax = DelayMax;
   //else if (DMax < DelayMin) DMax = DelayMin;
   //if (DMin > DelayMax) DMin = DelayMax;
   //else if (DMin < DelayMin) DMin = DelayMin;

   outfuncF = (((DMin + (DMax - DMin - DelayStep) * log(1.0 - localNetRandDev.Random() * (1.0 - TailNeglected)) 
   		 * InvLogTN) - DelayMin) / DelayStep);

   DPSNNverboseStart(true,1,0);
   if(outfuncF >= DelayNum){
     printf("ERROR: delay level generated %f > %f max allowed delay levels\n",outfuncF,DelayNum);
     fflush(stdout); exit(0);
   }
   DPSNNverboseEnd();

   outfunc = pErfLib->roundr2i(outfuncF*DelayStep);

   DPSNNverboseStart(true,1,0);
   if(outfunc >= neuSubPopParam[sourceSubPop].DMax[targetSubPop]){
     printf("ERROR: delay generated %d >= %f \n",outfunc,neuSubPopParam[sourceSubPop].DMax[targetSubPop]);
     fflush(stdout); exit(0);
   }
   DPSNNverboseEnd();

   return (outfunc);
}

#undef TailNeglected

uint32_t localNetClass::getRandomDelay_UNI(neuSubPopParamStruct *neuSubPopParam,neuSubPopEnum sourceSubPop,neuSubPopEnum targetSubPop)
{
   float outfuncF;
   float DMin,DMax;
   float DelayMin,DelayMax;
   float DelayStep,DelayNum;
   uint32_t outfunc;

   DMin = neuSubPopParam[sourceSubPop].DMin[targetSubPop];
   DMax = neuSubPopParam[sourceSubPop].DMax[targetSubPop];
   DelayMin = lnp_par.delayMin;
   DelayMax = lnp_par.delayMax;
   DelayNum = DSD__delayNum;

   DelayStep = (DelayMax - DelayMin)/DelayNum;
   DelayMax -= DelayStep/2;
   DelayMin += DelayStep/2;
   DMax = pErfLib->roundr2i((DMax - DelayMin)/DelayStep)*DelayStep+DelayMin;
   DMin = pErfLib->roundr2i((DMin - DelayMin)/DelayStep)*DelayStep+DelayMin;

   outfuncF = ((localNetRandDev.Random()*(DMax - DMin) + DMin - DelayMin) / DelayStep);

   DPSNNverboseStart(true,1,0);
   if(outfuncF >= DelayNum){
     printf("ERROR: delay level generated %f > %f max allowed delay levels\n",outfuncF,DelayNum);
     fflush(stdout); exit(0);
   }
   DPSNNverboseEnd();

   outfunc = pErfLib->roundr2i(outfuncF*DelayStep);

   DPSNNverboseStart(true,1,0);
   if(outfunc > neuSubPopParam[sourceSubPop].DMax[targetSubPop]){
     printf("ERROR: delay generated %d >= %f \n",outfunc,neuSubPopParam[sourceSubPop].DMax[targetSubPop]);
     fflush(stdout); exit(0);
   }
   DPSNNverboseEnd();

   return (outfunc);
}








void localNetClass::initLUT(){

  uint32_t sourceSubPop,targetSubPop,subPop;
  float J,DJ;

  for(subPop = 0; subPop < lnp_par.subPopNumber; subPop++){
    J=neuSubPopParam[subPop].JExt;
    DJ=neuSubPopParam[subPop].DJExt;
    pErfLib->makeGaussianLUT(tableLUT_synBathMatrix[subPop], ANALOG_DEPTH, J, J*DJ, 0.0, 2.0*J);
  }
  DPSNNverboseStart(false,1,0);
  {
    char *subPopName;
    char fileName[10];
    FILE *fp;
    uint32_t i;

    for(subPop = 0; subPop < lnp_par.subPopNumber; subPop++){
      subPopName = simpleCM_connectome.getSubPopName(simpleCM_connectome.getSubPopEnum(subPop));
      sprintf(fileName,"LUT_synBath_%s.txt",subPopName);
      fp = fopen(fileName,"w");
      for(i = 0; i < ANALOG_DEPTH; i++)
	fprintf(fp,"Table[%d] = %f \n",i,tableLUT_synBathMatrix[subPop][i]);
      fclose(fp);
    }
  }
  DPSNNverboseEnd();
  
  // Init tableLUT for the synaptic matrix
  for(sourceSubPop = 0; sourceSubPop < lnp_par.subPopNumber; sourceSubPop++)
    for(targetSubPop = 0; targetSubPop < lnp_par.subPopNumber; targetSubPop++){
      J=neuSubPopParam[sourceSubPop].J[targetSubPop];
      DJ=neuSubPopParam[sourceSubPop].DJ[targetSubPop];
      if(J>0)
	pErfLib->makeGaussianLUT(tableLUT_synMatrix[sourceSubPop][targetSubPop],ANALOG_DEPTH,J,J*DJ,0.0,2.0*J);
      else
	pErfLib->makeGaussianLUT(tableLUT_synMatrix[sourceSubPop][targetSubPop],ANALOG_DEPTH,J,-J*DJ,2.0*J,0.0);
    }

   DPSNNverboseStart(false,1,0);
   {
     char *sourceSubPopName;
     char *targetSubPopName;
     char fileName[10];
     FILE *fp;
     uint32_t i;

     for(sourceSubPop = 0; sourceSubPop < lnp_par.subPopNumber; sourceSubPop++)
       for(targetSubPop = 0; targetSubPop < lnp_par.subPopNumber; targetSubPop++){
	 sourceSubPopName = simpleCM_connectome.getSubPopName(simpleCM_connectome.getSubPopEnum(sourceSubPop));
	 targetSubPopName = simpleCM_connectome.getSubPopName(simpleCM_connectome.getSubPopEnum(targetSubPop));
	 sprintf(fileName,"LUT_%s%s.txt",sourceSubPopName,targetSubPopName);
	 fp = fopen(fileName,"w");
	 for(i = 0; i < ANALOG_DEPTH; i++)
	   fprintf(fp,"Table[%d] = %f \n",i,tableLUT_synMatrix[sourceSubPop][targetSubPop][i]);
	 fclose(fp);
       }
   }
   DPSNNverboseEnd();
  
}
