// DPSNN_writeIniFiles.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/sysinfo.h>

#include "DPSNN_debug.h"
#include "DPSNN_random.h"
#include "DPSNN_neuron.h"
#include "DPSNN_messagePassing.h"
#include "DPSNN_spike.h"
#include "DPSNN_chrono.h"
#include "DPSNN_memMeasure.h"
#include "DPSNN_localNet.h"
#include "DPSNN_LIFCAconnectome.h"
#ifdef LIFCAneuron
#include "erflib.h"
#include "randdev.h"
#endif


void localNetClass::writeIniFiles() {
    // This method generates the .ini files 
    // for Perseo simulator, both connections and modules files

  char fileNameC[30];
  char fileNameM[30];
  FILE *fpc,*fpm;
  uint32_t totalLocalCFT;
  uint32_t i_loc_cfx,i_loc_cfy;
  uint32_t this_cf,loc_cf,glob_source_cf,glob_target_cf;
  uint32_t stx,sty;
  uint32_t sourceGlobPop,targetGlobPop;
  float connectProb;
  double dmin,dmax,synJ,synDJ;
     
  int32_t target_cfx,target_cfy;
  uint32_t sourceSubPop,targetSubPop;
  uint32_t startSourcePop,startTargetPop;
  uint32_t count;
  double JExt,DJExt,CExt,NuExt,Tau,Theta,H,Tarp,AlphaC,TauC,gC;
  uint16_t NeuronInitType;
  char synType [6];

  if(floor(lnp_par.locCFT) == lnp_par.locCFT){
    printf("-650- init()  h=%03d globH=%d  locN=%d globCFT=%d locCFT=%f first_locCFX=%d last_locCFX=%d first_locCFY=%d last_locCFY=%d \n", 
	   lnp_par.loc_h, lnp_par.globH, lnp_par.locN,
	   lnp_par.globCFT, lnp_par.locCFT, 
	   lnp_par.first_locCFX, lnp_par.last_locCFX, 
	   lnp_par.first_locCFY, lnp_par.last_locCFY);
      
    totalLocalCFT = (uint32_t)floor(lnp_par.locCFT);
    sprintf(fileNameC,"PerseoConFile_h%d.txt",lnp_par.loc_h);
    sprintf(fileNameM,"PerseoModFile_h%d.txt",lnp_par.loc_h);
    fpc = fopen(fileNameC,"w");
    fpm = fopen(fileNameM,"w");

    for(i_loc_cfx = lnp_par.first_locCFX; 
	i_loc_cfx <= lnp_par.last_locCFX; i_loc_cfx++) {
      for(i_loc_cfy  = lnp_par.first_locCFY; 
	  i_loc_cfy <= lnp_par.last_locCFY; i_loc_cfy++) {

	this_cf = i_loc_cfx + i_loc_cfy * lnp_par.globCFX;
	loc_cf = this_cf % (uint32_t) lnp_par.locCFT;
	glob_source_cf = loc_cf + (uint32_t)lnp_par.locCFT * lnp_par.loc_h;
	  
	for(stx=0;stx<DSD__stencilX_Max;stx++) {	      
	  for(sty=0;sty<DSD__stencilY_Max;sty++){	      

	    switch(lnp_par.overallConnectivity) {
	      case explicitStencil: {
	        target_cfx = i_loc_cfx + stx - lnp_par.maxModDeltaX;
	        target_cfy = i_loc_cfy + sty - lnp_par.maxModDeltaY; };
	      break;
	      case homogeneous: {
	        target_cfx = stx;
	        target_cfy = sty; };
	      break;
	      default: {
	        printf("ERROR in writeIniFiles, unknown overallTopology\n");
	        fflush(stdout);	exit(0);
	      }
	    };


	      if((target_cfx>=0) && (target_cfx<(int32_t)lnp_par.globCFX) &&
	       (target_cfy>=0) && (target_cfy<(int32_t)lnp_par.globCFY)){

	      glob_target_cf = target_cfx + lnp_par.globCFX * target_cfy;

	      if((glob_target_cf<0) || (glob_target_cf>=lnp_par.globCFT) ||
		 (glob_source_cf<0) || (glob_source_cf>=lnp_par.globCFT)){
		printf("ERROR in Perseo .ini file generation: glob_target_cf not valid");
		fflush(stdout);exit(1);
	      }
		      
	      startSourcePop = lnp_par.subPopNumber*glob_source_cf;
	      startTargetPop = lnp_par.subPopNumber*glob_target_cf;

	      //for(sourceSubPop=0;sourceSubPop<lnp_par.subPopNumber;sourceSubPop++){
	        //for(targetSubPop=0;targetSubPop<lnp_par.subPopNumber;targetSubPop++){
	      for(targetSubPop=0;targetSubPop<lnp_par.subPopNumber;targetSubPop++){
		for(sourceSubPop=0;sourceSubPop<lnp_par.subPopNumber;sourceSubPop++){
		  sourceGlobPop = startSourcePop+sourceSubPop;
		  targetGlobPop = startTargetPop+targetSubPop;
		  connectProb = simpleCM_connectome.getConnectProb(sourceSubPop,targetSubPop,stx,sty);
		  dmin = neuSubPopParam[sourceSubPop].DMin[targetSubPop]+0.5;
		  dmax = neuSubPopParam[sourceSubPop].DMax[targetSubPop]+0.5;
		  synJ = neuSubPopParam[sourceSubPop].J[targetSubPop];
		  synDJ = neuSubPopParam[sourceSubPop].DJ[targetSubPop];
		  strcpy(synType,"Fixed");
		    
		  if(connectProb > 0.0){
		    fprintf(fpc," %5d %5d     %g   %g   %g   '%s'   %8g   %g\n",
			    targetGlobPop,sourceGlobPop,connectProb,dmin,dmax,synType,synJ,synDJ);
		  }
		}
	      }
	    }
	  }
	}

	//fprintf(fpm,"\n\n\nCFT: %d %d %d\n",glob_source_cf,i_loc_cfx,i_loc_cfy);
	for(sourceSubPop=0;sourceSubPop<lnp_par.subPopNumber;sourceSubPop++){
	  count = neuSubPopParam[sourceSubPop].count;
	  JExt = neuSubPopParam[sourceSubPop].JExt;
	  DJExt = neuSubPopParam[sourceSubPop].DJExt;
	  CExt = simpleCM_connectome.generateSourceNeuCExt(i_loc_cfx,i_loc_cfy, &lnp_par);	   
	  NuExt = neuSubPopParam[sourceSubPop].NuExt;
	  Tau = neuSubPopParam[sourceSubPop].Tau;
	  Theta = neuSubPopParam[sourceSubPop].Theta;
	  H = neuSubPopParam[sourceSubPop].H;
	  Tarp = neuSubPopParam[sourceSubPop].Tarp;
	  AlphaC = neuSubPopParam[sourceSubPop].AlphaC;
	  TauC = neuSubPopParam[sourceSubPop].TauC;
	  gC = neuSubPopParam[sourceSubPop].gC;
	  NeuronInitType = 0;
		    
	  fprintf(fpm," %6d   %g   %g   %g   %7g   %g   %g   %g   %g   %g   %4g   %4g   %d \n",
		  count,JExt,DJExt,CExt,NuExt,Tau,Theta,H,Tarp,AlphaC,TauC,gC,NeuronInitType);
	}
      }
    }


    fclose(fpc);
    fclose(fpm);
  }
};
