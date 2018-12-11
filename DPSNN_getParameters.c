// DPSNN_getParameters.c
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

//#include <string.h>
#include <time.h>
#include <unistd.h>

#include "DPSNN_getParameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_debug.h"
#include "DPSNN_LIFCAconnectome.h"

struct DPSNN_parameters getParameters(const uint32_t loc_h_initValue) {
  char *ps_env_globH; 
  char *ps_env_synGen;
  char *ps_env_CFX;
  char *ps_env_CFY;
  char *ps_env_D;
  char *ps_env_chrono;
  char *ps_env_totalSimTime_ms;
  char *ps_env_startPlasticity_ms;
  char *ps_env_stopPlasticity_ms;

  char *ps_env_moduloSec;
  char *ps_env_ratesSampling;
  char *ps_env_fastDebugDyn;
  char *ps_env_debugPrintEnable_ms;

  char *ps_env_globalSeed;
  char *ps_env_howManyOutputFiles;

  char *ps_env_startStatFiles_s;
  char *ps_env_stopStatFiles_s;

  char *ps_env_spikes_writeDisable;

  char *ps_env_statRatesPerPop_writeDisable;
  char *ps_env_startRatesPerPop_s;
  char *ps_env_stopRatesPerPop_s;
  
  char *ps_env_statSynPeriodicProbe_writeDisable;
  char *ps_env_startSynPeriodicProbe_s;
  char *ps_env_stopSynPeriodicProbe_s;

  char *ps_env_statSTDPevent_writeDisable;
  char *ps_env_startSTDPevent_s;
  char *ps_env_stopSTDPevent_s;

  char *ps_env_minSynWeight_f;

  char *ps_env_overallConnectivity;
  char *ps_env_stencilDim;

  char *ps_env_collectivesAlgo;
  char *ps_env_pktLength;

  //the following thalamic... is obsolete, should be removed
  /* char *ps_env_thalamicInputFreq; */
  
  struct DPSNN_parameters lnp_par;

   lnp_par.sequentialCode=false;

   lnp_par.loc_h = loc_h_initValue;

   //strcpy(lnp_par.hostName,getenv("HOSTNAME"));//name of computer
   gethostname(lnp_par.hostName,sizeof(lnp_par.hostName));
   DPSNNverboseStart(false,1,0);
   printf("Process H=%d runs on node %s \n",lnp_par.loc_h,lnp_par.hostName);
   DPSNNverboseEnd();

   //kind of overallConnectivity
   lnp_par.overallConnectivity = undefinedConnectivity;
   ps_env_overallConnectivity = getenv("env_overallConnectivity");
   if (ps_env_overallConnectivity!=NULL)
   {
     overallConnectivityEnum overallConnectivity;
     overallConnectivity=(overallConnectivityEnum)atoi(ps_env_overallConnectivity);
     switch (overallConnectivity) {
     case explicitStencil: {
       lnp_par.overallConnectivity = explicitStencil; }
       break;
     case homogeneous: {
       lnp_par.overallConnectivity = homogeneous;
       
       }
       break;
     case undefinedConnectivity:
     default: {
         if(lnp_par.loc_h==0) {
           printf(
	     "ERROR: invalid env_overallConnectivity param in DPSNN_script\n");
           fflush(stdout);exit(0);
	 };
       }
     }
   } else {
      if(lnp_par.loc_h==0) {
         printf("ERROR: missing env_overallConnectivity param in DPSNN_script\n");
         fflush(stdout);exit(0);
      };
   };

   //stencil dim for overallConnectivity
   ps_env_stencilDim = getenv("env_stencilDim");
   if (ps_env_stencilDim!=NULL) { 
     lnp_par.stencilDim = atoi(ps_env_stencilDim);
     if(lnp_par.stencilDim<=0 || lnp_par.stencilDim>DSD__stencilX_Max) {
        if(lnp_par.loc_h==0) {
          printf("ERROR: env_stencilDim = %d out of %d range in DPSNN_script \n",
		lnp_par.stencilDim, DSD__stencilX_Max);
           fflush(stdout);exit(0);
	};
     }
   }else{
      if(lnp_par.loc_h==0) {
        printf("ERROR: missing env_stencilDim param in DPSNN_script \n");
        fflush(stdout);exit(0);
      };
   };
  
   switch(lnp_par.overallConnectivity) {
      case explicitStencil: {
        lnp_par.maxModDeltaX = lnp_par.stencilDim-1; 
        lnp_par.maxModDeltaY = lnp_par.stencilDim-1;
        lnp_par.stencilX_Max = 2 * lnp_par.maxModDeltaX + 1;
        lnp_par.stencilY_Max = 2 * lnp_par.maxModDeltaY + 1;
        lnp_par.bathEfficacyTemplateX_Max=lnp_par.stencilDim;
        lnp_par.bathEfficacyTemplateY_Max=lnp_par.stencilDim;
        if(lnp_par.stencilX_Max>DSD__stencilX_Max) {
	  if(lnp_par.loc_h==0) {
	     printf("ERROR: stencilX_Max = %d out of %d range, due to env_stencilDim set in DPSNN_script \n",
		lnp_par.stencilX_Max, DSD__stencilX_Max);
	     fflush(stdout);exit(0);
	  };
        };
        if(lnp_par.stencilY_Max>DSD__stencilY_Max) {
	  if(lnp_par.loc_h==0) {
	     printf("ERROR: env_stencilY_Max = %d out of %d range, due to env_stencilDim set in DPSNN_script \n",
		lnp_par.stencilY_Max, DSD__stencilY_Max);
	     fflush(stdout);exit(0);
	  };
        };	  
      };
      break;
      case homogeneous: {
        lnp_par.maxModDeltaX = 0; 
        lnp_par.maxModDeltaY = 0;
        lnp_par.stencilX_Max = lnp_par.stencilDim;
        lnp_par.stencilY_Max = lnp_par.stencilDim;
        lnp_par.bathEfficacyTemplateX_Max=lnp_par.stencilDim;
        lnp_par.bathEfficacyTemplateY_Max=lnp_par.stencilDim;
       
        if(lnp_par.stencilX_Max>DSD__stencilX_Max) {
	   if(lnp_par.loc_h==0) {
	      printf("ERROR: stencilX_Max = %d out of %d range, due to env_stencilDim in DPSNN_script \n",
		lnp_par.stencilX_Max, DSD__stencilX_Max);
	      fflush(stdout);exit(0);
	   };
        };
        if(lnp_par.stencilY_Max>DSD__stencilY_Max) {
	  if(lnp_par.loc_h==0) {
	    printf("ERROR: env_stencilY_Max = %d out of %d range, due to env_stencilDim in DPSNN_script \n",
		lnp_par.stencilY_Max, DSD__stencilY_Max);
	    fflush(stdout);exit(0);
	  };
        };
      };
      break;
      case undefinedConnectivity:
      default: {
	  if(lnp_par.loc_h==0) {      
            printf(
	    "ERROR: invalid env_overallConnectivity param in DPSNN_script\n");
            fflush(stdout);exit(0);
	  };
      };
    };
   
   ps_env_globH = getenv("env_globH");//number of processes
   if (ps_env_globH != NULL) { 
        lnp_par.globH = atoi(ps_env_globH);
	if(lnp_par.globH<=0 || lnp_par.globH> DSD__maxGlobH ) {
	  if(lnp_par.loc_h==0) {
	     printf("ERROR: env_globH = %d out of range: [from 1 to %d] in DPSNN_script \n",
		  lnp_par.globH, DSD__maxGlobH);
              fflush(stdout);exit(0);
	  };
        }
   }else{
      if(lnp_par.loc_h==0) {
        printf("ERROR: undefined env_globH param in DPSNN_script\n");
        fflush(stdout);exit(0);
      };
   };

   //what kind of generation of synapses?
   lnp_par.synGen = undefinedSynGen;
   ps_env_synGen = getenv("env_synGen");
   if (ps_env_synGen!=NULL) {
     synGenEnum synGen;
     synGen = (synGenEnum) atoi(ps_env_synGen);
     switch (synGen) {
        case simpleCorticalModule_synGen_2: {
           lnp_par.synGen = simpleCorticalModule_synGen_2;
        }
        break;
        default: {
	  if(lnp_par.loc_h==0) {
	     printf("ERROR: unsupported synGen %d requested by DPSNN_script\n",
		 (int)synGen);
             fflush(stdout);exit(0);
	  };
        };
      };
   } else {
      if(lnp_par.loc_h==0) {
         printf("ERROR: undefined env_synGen param in DPSNN_script\n");
         fflush(stdout);exit(0);
      };
   };
     
   readColumnFile(&lnp_par.subPopNumber,&lnp_par.neuronsPerCM);

   if (lnp_par.neuronsPerCM<=0) {
      if(lnp_par.loc_h==0) {
        printf("ERROR: lnp_par.neuronsPerCM <=0 in DPSNN_getParameters\n");
         fflush(stdout);exit(0);
      };
   };
   
   ps_env_CFX = getenv("env_CFX");
   ps_env_CFY = getenv("env_CFY");
   
   if ((ps_env_CFX!=NULL) && (ps_env_CFY!=NULL)) { 
      lnp_par.globCFX = atoi(ps_env_CFX);
      lnp_par.globCFY = atoi(ps_env_CFY);
      lnp_par.globCFT = lnp_par.globCFX * lnp_par.globCFY;
      if(lnp_par.globCFT > DSD__maxGlobCFT){
	if(lnp_par.loc_h==0) {
	  printf(
		  "ERROR: globCFT more then %d\n",DSD__maxGlobCFT);
	   fflush(stdout);
	   exit(0);
	};
       }
    }else{
       if(lnp_par.loc_h==0) {
         printf(
	 "ERROR: please set the env_CFX and env_CFY values in DPSNN_script (grid of Cortical Modules)\n");
         fflush(stdout);
         exit(0);
       };
    };
   
    if (((lnp_par.globH > lnp_par.globCFT)         &&
           ((lnp_par.globH % lnp_par.globCFT) != 0 )
	       ) || (
	    (lnp_par.globCFT > lnp_par.globH)         &&
	   ((lnp_par.globCFT % lnp_par.globH) != 0 )
	  ))
    {
       if(lnp_par.loc_h==0) {
         printf(
         "ERROR:#processes and #cort. columns must be  multiple OR equal\n");
         fflush(stdout);exit(0);
       };
    };

    lnp_par.locCFT = ((float)lnp_par.globCFT/(float)lnp_par.globH);

    if(lnp_par.locCFT < 1.0) {
      if(lnp_par.globCFT >= 64) {
        if(lnp_par.loc_h==0) {
	    printf(
	    "ERROR:not foreseen (?why?) fractioning of system with more than  64 CFT\n");
	       fflush(stdout);exit(0);
	};
      };
      
      if(lnp_par.globCFT < 64) {
	if (lnp_par.locCFT < 1.0/DSD__maxCMFract) {
	   if(lnp_par.loc_h==0) {
	       printf(
		    "ERROR:a Cortical Module can be divided by max %d processes\n",DSD__maxCMFract);
	       fflush(stdout);exit(0);
	   };
	};
      };
	   
      if(!(((lnp_par.locCFT==(1.0/2.0))&&((lnp_par.neuronsPerCM%2)==0)) ||
		((lnp_par.locCFT==(1.0/4.0))&&((lnp_par.neuronsPerCM%4)==0))||
		((lnp_par.locCFT==(1.0/8.0))&&((lnp_par.neuronsPerCM%8)==0))||
		((lnp_par.locCFT==(1.0/16.0))&&((lnp_par.neuronsPerCM%16)==0))||
		((lnp_par.locCFT==(1.0/32.0))&&((lnp_par.neuronsPerCM%32)==0))||
		((lnp_par.locCFT==(1.0/64.0))&&((lnp_par.neuronsPerCM%64)==0))||
		((lnp_par.locCFT==(1.0/128.0))&&((lnp_par.neuronsPerCM%128)==0)))) 
      {
	 if(lnp_par.loc_h==0) {
	      printf(
		 "ERROR if globH>globCFT,neuronsPerCM and locCFT must be (sub)multiple of 2^n \n");
	          fflush(stdout);exit(0);
	 };
       };
    };// end of code managing fractioning  columns among processes 

    lnp_par.locN = 
	 (uint32_t)(lnp_par.neuronsPerCM * lnp_par.locCFT);
    lnp_par.globN = lnp_par.locN * lnp_par.globH;

    DPSNNverboseStart(true,1,0);
       if(lnp_par.loc_h==0) {
         printf("getParameters(): total neurons (globN)       =%d \n",
		lnp_par.globN);
	 printf("getParameters(): local neurons (locN)        =%d \n",
		lnp_par.locN);
	 printf("getParameters(): neurons per cortical module =%d \n",
		lnp_par.neuronsPerCM);
	 printf("getParameter():  overallConnectivity         =%d \n",
		lnp_par.overallConnectivity);
	 printf("getParameter():  stencilDim                  =%d \n",
		lnp_par.stencilDim);
       };
    DPSNNverboseEnd();
       
    {//first and last neuron on this process
	 lnp_par.first_glob_n = lnp_par.locN * lnp_par.loc_h;
	 lnp_par.last_glob_n = lnp_par.first_glob_n + lnp_par.locN - 1;
    };
       
    { //computes the indexes of first and last CM on this process
         lnp_par.first_locCFX = 
	   (lnp_par.first_glob_n / lnp_par.neuronsPerCM) % 
	    lnp_par.globCFX;  
         lnp_par.first_locCFY = 
           (lnp_par.first_glob_n / lnp_par.neuronsPerCM) / 
	   lnp_par.globCFX;
         lnp_par.last_locCFX = 
	   (lnp_par.last_glob_n / lnp_par.neuronsPerCM) % 
	    lnp_par.globCFX;  
         lnp_par.last_locCFY = 
           (lnp_par.last_glob_n / lnp_par.neuronsPerCM) / 
	   lnp_par.globCFX;
   };

   if(lnp_par.locN > DSD__maxLocN) {
     if(lnp_par.loc_h==0) {
        printf("ERROR too many neurons requested \n");
	  fflush(stdout);
	  exit(0);
     };
   };

   //minimum value of synapses (for conversion/storage in int16)
   ps_env_minSynWeight_f = getenv("env_minSynWeight_f");
   if (ps_env_minSynWeight_f!=NULL)
   {
     lnp_par.minSynWeight_f = atof(ps_env_minSynWeight_f);
     if((lnp_par.minSynWeight_f >=0) ||
	(lnp_par.minSynWeight_f < DSD__acceptable_minSynWeigth_f))
     {
        if(lnp_par.loc_h==0) {
          printf("env_minSynWeigth out of range %f - 0.0\n",
	      DSD__acceptable_minSynWeigth_f);
          fflush(stdout); exit(0);
	};
     }else{
       lnp_par.factorFloat_2_weightType = DSD__INT16_MIN / lnp_par.minSynWeight_f;
       lnp_par.factorWeightType_2_Float = lnp_par.minSynWeight_f / DSD__INT16_MIN;
     };
   }else{
       if(lnp_par.loc_h==0) {
          printf("MISSING env_minSynWeigth in DPSNN launch script\n");
          fflush(stdout);
          exit(0);
       };
   };

   DPSNNverboseStart(true,1,0);
       if(lnp_par.loc_h==0) {
         printf("getParameters(): accepted synrange (minSynWeight_f)=%f \n",
		lnp_par.minSynWeight_f);
       };
   DPSNNverboseEnd();

   if(lnp_par.globCFT * lnp_par.neuronsPerCM !=
      lnp_par.globN) {
         if(lnp_par.loc_h==0) {
            printf("ERROR globN=%d not equal globCFT=%d * neuronsPerCM=%d\n",
	    lnp_par.globN, lnp_par.globCFT, lnp_par.neuronsPerCM);
	 };
   };

   DPSNNverboseStart(false,1,0);
     if(lnp_par.loc_h == 0) {
       printf(
       "getenv:lnp_par.globCFT=%d,globCFX=%d,globCFY=%d\n", 
	  lnp_par.globCFT, lnp_par.globCFX, lnp_par.globCFY);
       printf(
       "getenv:lnp_par.neuronsPerCM=%d,locN=%d,locCFT=%f\n", 
	  lnp_par.neuronsPerCM,lnp_par.locN,lnp_par.locCFT);
       printf(
       "getenv:lnp_par.globH=%d\n", 
	  lnp_par.globH);
     };
   DPSNNverboseEnd();
   
   if((lnp_par.neuronsPerCM * 
      lnp_par.globCFX * lnp_par.globCFY) !=
      lnp_par.globN) {
     	  if(lnp_par.loc_h==0) {
             printf("ERROR neuronsPerCM * globCFX * globCFY != globN)\n");
             fflush(stdout);
             exit(0);
	  };
   };

   //set lnp_par.tabNeuron
   //used in random generation reset 
   //during neuronal dynamic
   {
     uint32_t maxNumOfProc;
     uint32_t CMDivisor;
     
     if(lnp_par.globCFT<64){
       //Finds the max 2^N-->CMDivisor that is a divider of neuronsPerCM
       CMDivisor = lnp_par.neuronsPerCM & (-lnp_par.neuronsPerCM);
       if(CMDivisor > DSD__maxCMFract)
	 CMDivisor = DSD__maxCMFract;
     } else {
       //if there are more than 64 colums, we do not divide colums
    	 CMDivisor = 1;
     };

     maxNumOfProc = lnp_par.globCFT * CMDivisor;
     lnp_par.tabNeuron = lnp_par.globN/maxNumOfProc;
     
     
     //lnp_par.tabNeuron = lnp_par.globN/DSD__maxGlobH;

     DPSNNverboseStart(true,1,0);
     if(lnp_par.loc_h==0) {
       printf("getParameters(): on H=%d tabNeuron=%d \n",
	      lnp_par.loc_h,lnp_par.tabNeuron);
     };
     DPSNNverboseEnd();
   };

   //what is range of possible synaptic delays 
   //(default 1-20 ms) 
   ps_env_D = getenv("env_D");
   if (ps_env_D!=NULL) { 
     lnp_par.D = atoi(ps_env_D);
   } else {
          if(lnp_par.loc_h==0) {
            printf("ERROR: missing env_D param in DPSNN_script mpirun invocation\n");
            fflush(stdout);exit(0);
	  };
   };

   //performance chrono activated?
   ps_env_chrono = getenv("env_chrono");
   if (ps_env_chrono!=NULL) { 
     lnp_par.chrono = atoi(ps_env_chrono);
   } else {  
     lnp_par.chrono = 1;
   }
  
   //how many ms of simulated time? (0->synapse generation only)
   ps_env_totalSimTime_ms = getenv("env_totalSimTime_ms");
   if (ps_env_totalSimTime_ms!=NULL) { 
     lnp_par.totalSimTime_ms = atoi(ps_env_totalSimTime_ms);
   }else{
       if(lnp_par.loc_h==0) {
         printf("ERROR: missing env_totalSimTime_ms param in DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };

   //start time for plasticity
   ps_env_startPlasticity_ms = getenv("env_startPlasticity_ms");
   if(ps_env_startPlasticity_ms!=NULL) {
      lnp_par.startPlasticity_ms = atoi(ps_env_startPlasticity_ms);
   }else{
     printf("ERROR: missing env_startPlasticity_ms param in DPSNN_script mpirun invocation\n");
     fflush(stdout);exit(0);
   }

   //stop time for plasticity
   ps_env_stopPlasticity_ms = getenv("env_stopPlasticity_ms");
   if(ps_env_stopPlasticity_ms!=NULL) {
      lnp_par.stopPlasticity_ms = atoi(ps_env_stopPlasticity_ms);
   }else{
       if(lnp_par.loc_h==0) {
         printf("ERROR: missing env_stopPlasticity_ms param in DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   }
   
   DPSNNverboseStart(true,1,0);
   if(lnp_par.loc_h==0) {
       printf("getParameters():         startPlasticity_ms=%d\n",
	      lnp_par.startPlasticity_ms);
       printf("getParameters():         stopPlasticity_ms=%d\n",
	      lnp_par.stopPlasticity_ms);
   };
   DPSNNverboseEnd();

   //what is the interval in seconds between statistic printing 
   //(default could be 60 s, now is set to 1) 
   ps_env_moduloSec = getenv("env_moduloSec");
   if (ps_env_moduloSec!=NULL) { 
     lnp_par.moduloSec = atoi(ps_env_moduloSec);
   } else {
       if(lnp_par.loc_h==0) {
         printf("ERROR: missing env_moduloSec param in DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };

   //when, in general, the print of statistic files should start?
   //this can be further restricted by other options
   ps_env_startStatFiles_s = getenv("env_startStatFiles_s");
   if (ps_env_startStatFiles_s!=NULL) { 
     lnp_par.startStatFiles_s = atoi(ps_env_startStatFiles_s);
   } else {
        if(lnp_par.loc_h==0) {
          printf("ERROR: missing env_startStatFiles_s param in DPSNN_script mpirun invocation\n");
          fflush(stdout);exit(0);
	};
   };

   //when, in general, the print of statistic files should end?
   //this can be further restricted by other options
   ps_env_stopStatFiles_s = getenv("env_stopStatFiles_s");
   if (ps_env_stopStatFiles_s!=NULL) { 
     lnp_par.stopStatFiles_s = atoi(ps_env_stopStatFiles_s);
   } else {
        if(lnp_par.loc_h==0) {

           printf("ERROR: missing env_stopStatFiles_s param in DPSNN_script mpirun invocation\n");
           fflush(stdout);exit(0);
	};
   };

   //should we print any file of spikes
   lnp_par.spikes_writeDisable = true;
   ps_env_spikes_writeDisable = getenv("env_spikes_writeDisable");
   if (ps_env_spikes_writeDisable!=NULL)
   {
     uint32_t truth;
     truth=atoi(ps_env_spikes_writeDisable);
     switch (truth) {
       case 0: lnp_par.spikes_writeDisable = false; break;
       case 1: lnp_par.spikes_writeDisable = true; break;
       default: {
          if(lnp_par.loc_h==0) {
            printf(
	       "ERROR: invalid env_spikes_writeDisable param in DPSNN_script mpirun invocation\n");
            fflush(stdout);exit(0);
          };
       };
     };
   } else {
      if(lnp_par.loc_h==0) {
        printf("ERROR: missing env_spikes_writeDisable param in DPSNN_script mpirun invocation\n");
        fflush(stdout);exit(0);
      };
   };

   //should we print any file of STDP event lists
   lnp_par.statSTDPevent_writeDisable = true;
   ps_env_statSTDPevent_writeDisable = getenv("env_statSTDPevent_writeDisable");
   if (ps_env_statSTDPevent_writeDisable!=NULL)
   {
     uint32_t truth;
     truth=atoi(ps_env_statSTDPevent_writeDisable);
     switch (truth) {
       case 0: lnp_par.statSTDPevent_writeDisable = false; break;
       case 1: lnp_par.statSTDPevent_writeDisable = true; break;
       default: {
           if(lnp_par.loc_h==0) {
              printf(
	         "ERROR: invalid env_statSTDPevent_writeDisable param in DPSNN_script mpirun invocation\n");
             fflush(stdout);exit(0);
	   };
       };
     };
   } else {
      if(lnp_par.loc_h==0) {
        printf("ERROR: missing env_statSTDPevent_writeDisable param DPSNN_script mpirun invocation\n");
        fflush(stdout);exit(0);
      };
   };

   //should we print any file of syn periodic probe
   lnp_par.statSynPeriodicProbe_writeDisable = true;
   ps_env_statSynPeriodicProbe_writeDisable = getenv("env_statSynPeriodicProbe_writeDisable");
   
   if (ps_env_statSynPeriodicProbe_writeDisable!=NULL)
   {
     uint32_t truth;
     truth=atoi(ps_env_statSynPeriodicProbe_writeDisable);
     switch (truth) {
       case 0: lnp_par.statSynPeriodicProbe_writeDisable = false; break;
       case 1: lnp_par.statSynPeriodicProbe_writeDisable = true; break;
       default: {
	   if(lnp_par.loc_h==0) {
             printf(
	        "ERROR: invalid env_statSynPeriodicProbe_writeDisable param in DPSNN_script mpirun invoc\n");
             fflush(stdout);exit(0);
	   };
       };
     };
   } else {
     if(lnp_par.loc_h==0) {
       printf(
         "ERROR: missing env_statSynPeriodicProbe_writeDisable param DPSNN_script mpirun invocation\n");
       fflush(stdout);exit(0);
     };
   };

  //enabling print of ratesPerPop files
   ps_env_statRatesPerPop_writeDisable = getenv("env_statRatesPerPop_writeDisable");
   if (ps_env_statRatesPerPop_writeDisable!=NULL) { 
     uint32_t truth;
     truth=atoi(ps_env_statRatesPerPop_writeDisable);
     switch (truth) {
       case 0: lnp_par.statRatesPerPop_writeDisable = false; break;
       case 1: lnp_par.statRatesPerPop_writeDisable = true; break;
       default: {
	   if(lnp_par.loc_h==0) {
             printf(
	      "ERROR: invalid env_statRatesPerPop_writeDisable param in DPSNN_script mpirun invoc\n");
              fflush(stdout);exit(0);
	   };
       };
     }
   } else {
       if(lnp_par.loc_h==0) {
         printf(
         "ERROR: missing env_statRatesPerPop_writeDisable param DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };

  // env_startRatesPerPop_s 
   ps_env_startRatesPerPop_s = getenv("env_startRatesPerPop_s");
   if (ps_env_startRatesPerPop_s!=NULL) { 
     lnp_par.startRatesPerPop_s = atoi(ps_env_startRatesPerPop_s);
   } else {
       if(lnp_par.loc_h==0) {
         printf(
         "ERROR: missing env_startRatesPerPop_s param DPSNN_script mpirun invocation\n");
          fflush(stdout);exit(0);
       };
   };
  
   // env_stopRatesPerPop_s 
   ps_env_stopRatesPerPop_s = getenv("env_stopRatesPerPop_s");
   if (ps_env_stopRatesPerPop_s!=NULL) { 
     lnp_par.stopRatesPerPop_s = atoi(ps_env_stopRatesPerPop_s);
   } else {
       if(lnp_par.loc_h==0) {
          printf(
             "ERROR: missing env_stopRatesPerPop_s param DPSNN_script mpirun invocation\n");
          fflush(stdout);exit(0);
       };
   };

   
   // env_startSynPeriodicProbe_s 
   ps_env_startSynPeriodicProbe_s = getenv("env_startSynPeriodicProbe_s");
   if (ps_env_startSynPeriodicProbe_s!=NULL) { 
     lnp_par.startSynPeriodicProbe_s = atoi(ps_env_startSynPeriodicProbe_s);
   } else {
      if(lnp_par.loc_h==0) {
        printf(
        "ERROR: missing env_startSynPeriodicProbe_s param DPSNN_script mpirun invocation\n");
        fflush(stdout);exit(0);
      };
   };
  
   // env_stopSynPeriodicProbe_s 
   ps_env_stopSynPeriodicProbe_s = getenv("env_stopSynPeriodicProbe_s");
   if (ps_env_stopSynPeriodicProbe_s!=NULL) { 
     lnp_par.stopSynPeriodicProbe_s = atoi(ps_env_stopSynPeriodicProbe_s);
   } else {
       if(lnp_par.loc_h==0) {
         printf(
         "ERROR: missing env_stopSynPeriodicProbe_s param DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };

   // env_startSTDPevent_s 
   ps_env_startSTDPevent_s = getenv("env_startSTDPevent_s");
   if (ps_env_startSTDPevent_s!=NULL) { 
     lnp_par.startSTDPevent_s = atoi(ps_env_startSTDPevent_s);
   } else {
       if(lnp_par.loc_h==0) {
         printf(
         "ERROR: missing env_startSTDPevent_s param DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };
  
   // env_stopSTDPevent_s 
   ps_env_stopSTDPevent_s = getenv("env_stopSTDPevent_s");
   if (ps_env_stopSTDPevent_s!=NULL) { 
     lnp_par.stopSTDPevent_s = atoi(ps_env_stopSTDPevent_s);
   } else {
       if(lnp_par.loc_h==0) {
         printf(
         "ERROR: missing env_stopSTDPevent_s param DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };

   //what is the interval in milliseconds between ratesPerPopulations printing 
   //(default 5 ms) 
   ps_env_ratesSampling = getenv("env_ratesSampling");
   if (ps_env_ratesSampling!=NULL) { 
     lnp_par.ratesSampling = atoi(ps_env_ratesSampling);
   } else {
      if(lnp_par.loc_h==0) {
         printf(
	 "ERROR: missing env_ratesSampling param DPSNN_script mpirun invocation\n");
         fflush(stdout);exit(0);
       };
   };

   //triggers a fast dynamic behaviour for easier debugging
   //(default 0) 
   ps_env_fastDebugDyn = getenv("env_fastDebugDyn");
   if (ps_env_fastDebugDyn!=NULL) { 
     lnp_par.fastDebugDyn = (fastDebugDynEnum)atoi(ps_env_fastDebugDyn);
   } else {
     lnp_par.fastDebugDyn = default_nonActiveFastDebugDyn_0;
   };

   //enabling debug prints starting from a given ms
   //(default no 0) 
   ps_env_debugPrintEnable_ms = getenv("env_debugPrintEnable_ms");
   if (ps_env_debugPrintEnable_ms!=NULL) { 
     lnp_par.debugPrintEnable_ms = 
      atoi(ps_env_debugPrintEnable_ms);
   } else {
     lnp_par.debugPrintEnable_ms = 0;
   };

   
   //set the globalSeed for random number generator initalization
   //lnp_par.globalSeed = lnp_par.loc_h+1;
   {
     int32_t tempSeed;
     ps_env_globalSeed = getenv("env_globalSeed");
     tempSeed = atoi(ps_env_globalSeed);
     if (tempSeed < 0) { 
       lnp_par.globalSeed = uint32_t(time(NULL));
     } else {
       lnp_par.globalSeed = tempSeed;
     }     
   }
   
   //set output mode 
   //0-->short output, only 3 files, start, lasteven, lastodd sec
   //1-->long output, 1 file/sec + 1 global file
   ps_env_howManyOutputFiles = getenv("env_howManyOutputFiles");
   lnp_par.howManyOutputFiles =
     (howManyOutputFilesEnum)atoi(ps_env_howManyOutputFiles);
   if ((ps_env_howManyOutputFiles==NULL) ||
       (lnp_par.howManyOutputFiles > longOutput_1)) { 
     lnp_par.howManyOutputFiles = default_shortOutput_0;;
   };

   readPlasticityParametersFile(&lnp_par);

   //the following three lines are obsolete, should be removed
   /* lnp_par.thalamicInput = linear_neuronsPerCM_vs_time_thalamicInput_3; */
   /* lnp_par.thalamicInputFreq = 8; */
   /* lnp_par.M = 100; */

   // env_collectivesAlgo 
   ps_env_collectivesAlgo = getenv("env_collectivesAlgo");
   if (ps_env_collectivesAlgo!=NULL) { 
     lnp_par.collectivesAlgo = (collectivesAlgoEnum)atoi(ps_env_collectivesAlgo);
   } else {
      if(lnp_par.loc_h==0) {
        printf(
        "ERROR: missing env_collectivesAlgo param DPSNN_script mpirun invocation\n");
        fflush(stdout);exit(0);
      };
   };

   DPSNNverboseStart(true,1,0);
      if(lnp_par.loc_h==0) {	    
        printf("collectivesAlgo = %d\n", lnp_par.collectivesAlgo);
      };
   DPSNNverboseEnd();

   // env_pktLength
   ps_env_pktLength = getenv("env_pktLength");
   if (ps_env_pktLength!=NULL) {
     lnp_par.pktLength = (uint32_t)atoi(ps_env_pktLength);
   } else {
      if(lnp_par.loc_h==0) {
        printf(
        "ERROR: missing env_pktLength param DPSNN_script mpirun invocation\n");
        fflush(stdout);exit(0);
      };
   };

   DPSNNverboseStart(true,1,0);
      if(lnp_par.loc_h==0) {
        printf("pktLength = %d\n", lnp_par.pktLength);
      };
   DPSNNverboseEnd();

   return(lnp_par);
};


void readColumnFile(uint32_t *pSubPopNumber,uint32_t *pNeuronsPerCM) {
#define maxNumOfParam 2
#define bufferSize 200

  char fileName[30];
  FILE *fp;
  uint32_t x,y;
  char buffer[bufferSize];
  char *token;
  uint32_t numOfTokenInLine;
  float initValue[maxNumOfParam];
  uint32_t subPopNumber;
  uint32_t neuronsPerCM;

  //Do we need some initialization to default values??
  
  //sprintf(fileName,"column_h%d.txt",lnp_par.loc_h);
  sprintf(fileName,"column.txt");
  fp = fopen(fileName,"r");
  if (fp == NULL) {
    printf("ERROR opening file column.txt \n");
    fflush(stdout);exit(0);
  }

  x=0;
  subPopNumber = 0;
  neuronsPerCM = 0;
  while(fgets(buffer,bufferSize,fp)!=NULL) {
    initValue[0] = 0.0;
    initValue[1] = 0.0;
    numOfTokenInLine = 0;
    token = strtok(buffer," \t\n");
    y=0;
    while (token != NULL){
      if (token[0] == '#') break;
      if (numOfTokenInLine >= maxNumOfParam) break;
      numOfTokenInLine++;
      initValue[y] = atof(token);
      y++;
      token = strtok (NULL," \t\n");
    }
    if(numOfTokenInLine>0){
      x++;
      neuronsPerCM += (uint32_t) initValue[1];
    }
    if ((token[0] != '#')&&(numOfTokenInLine>0))
      subPopNumber++;
    if(subPopNumber > neuSubPopTotal){
      printf("ERROR in subPop initialization: wrong subPopulation total number in column.txt file \n");
      fflush(stdout);exit(0);
    } 
  }
  DPSNNverboseStart(false,1,0);
  if(x!=neuSubPopTotal){
    printf("WARNING in subPop initialization: number of subPops in column.txt file doesn't match the predefined number %d \n",neuSubPopTotal);
  }
  DPSNNverboseEnd();
  fclose(fp);
  *pSubPopNumber = subPopNumber;
  *pNeuronsPerCM = neuronsPerCM;
  DPSNNverboseStart(false,1,0);
  printf("From initialization file: %d subPop/column - %d neu/column \n",subPopNumber,neuronsPerCM);
  DPSNNverboseEnd();
#undef maxNumOfParam
}

void readPlasticityParametersFile(struct DPSNN_parameters *pLnp_par) {
  #define maxNumOfPlasticityParam 20
  #define bufferSize 200

  char fileName[30];
  FILE *fp;
  char buffer[bufferSize];
  char *token;
  uint32_t numOfTokenInLine;
  uint32_t numOfNonRemarkLines;
  float initValue[maxNumOfPlasticityParam];
  uint32_t y,j;
   
  sprintf(fileName,"plasticity.txt");
  fp = fopen(fileName,"r");
  if (fp == NULL) {
    if(pLnp_par->loc_h==0) {
      printf("WARNING: Plasticity.txt miss -set default Song STDP causal params \n");
      fflush(stdout);
    };
    pLnp_par->plasticityAlgo = SongSTDPenum;
    pLnp_par->SongSTDP_param.antiCausalCoeff=-1.0;
    pLnp_par->SongSTDP_param.STDPmultiplier=0.001;
    pLnp_par->SongSTDP_param.tauSTDP_ms=20.0;
    pLnp_par->SongSTDP_param.causalSTDPmaxAbsValue=0.1;
    pLnp_par->SongSTDP_param.antiCausalSTDPmaxAbsValue=0.12;
    pLnp_par->SongSTDP_param.tauDerivativeDecay_s=1.0;
    pLnp_par->SongSTDP_param.derivativeDrift=0.0;
    return;
  }
  else
  {
    for(j=1;j<maxNumOfPlasticityParam;j++)
      initValue[j] = 0.0;
    
    numOfNonRemarkLines=0;

    while(fgets(buffer,bufferSize,fp)!=NULL) {
      //initValue[0] = p_lnp_par->subPopNumber;
      numOfTokenInLine = 0;
      token = strtok(buffer," \t\n");
      y=0;
      while (token != NULL){
        if (token[0] == '#') break;
        numOfTokenInLine++;
        initValue[y] = atof(token);
        y++;
        token = strtok (NULL," \t\n");
      }
      if (numOfNonRemarkLines>0) break;
    }
    
    fclose(fp);

    if((numOfTokenInLine==0)||(numOfTokenInLine >= maxNumOfPlasticityParam)) {
      if(pLnp_par->loc_h==0) {
        printf("ERROR: Plasticity.txt %d columns (empty or too much params)\n",
	       numOfTokenInLine);
        fflush(stdout); exit(0);
      }
    }
      
    if(initValue[0] != SongSTDPenum) {
      if(pLnp_par->loc_h==0) {
        printf("ERROR: uknown algo in Plasticity.txt (first col)\n");
        fflush(stdout); exit(0);
      }
    }
      
    if(initValue[0]==SongSTDPenum) {
      if(numOfTokenInLine!=(numOfSTDPparamsEnum)) {
        if(pLnp_par->loc_h==0) {	  
          printf("ERROR: wrong number %d columns Plasticity.txt (expected %d)\n",
	     numOfTokenInLine, numOfSTDPparamsEnum+1);
             fflush(stdout); exit(0);
	}
      }
    }
	
    pLnp_par->plasticityAlgo=SongSTDPenum;
    {
      uint32_t counter;
      counter=0;
	 
      pLnp_par->SongSTDP_param.STDPmultiplier =
	               initValue[STDPmultiplierEnum];
      counter++;
	  
      pLnp_par->SongSTDP_param.antiCausalCoeff =
	               initValue[antiCausalCoeffEnum];
      counter++;
	  
      pLnp_par->SongSTDP_param.tauSTDP_ms =
	               initValue[tauSTDP_msEnum];
      counter++;
	  
      pLnp_par->SongSTDP_param.causalSTDPmaxAbsValue =
		       initValue[causalSTDPmaxAbsValueEnum];
      counter++;
	  
      pLnp_par->SongSTDP_param.antiCausalSTDPmaxAbsValue =
		       initValue[antiCausalSTDPmaxAbsValueEnum];
      counter++;
	  
      pLnp_par->SongSTDP_param.tauDerivativeDecay_s =
		       initValue[tauDerivativeDecay_sEnum];
      counter++;
	  
      pLnp_par->SongSTDP_param.derivativeDrift =
		       initValue[derivativeDriftEnum];
      counter++;
      //-1 because the first column is for algo	  
      if(counter != (numOfSTDPparamsEnum-1)) {
	if(pLnp_par->loc_h==0) {
	  printf("ERROR missing STDP params: initialized %d, expected %d\n",
		 counter,numOfSTDPparamsEnum);
	  fflush(stdout); exit(0);
	}
      }
    }
	
    DPSNNverboseStart(true,1,0);
      if(pLnp_par->loc_h==0) {
        printf("PLASTICITY PARAMS from file Plasticity.txt:\n");
	    
        printf("Plasticiy Algo = SongSTDPenum\n");
	
        printf("STDPmultiplier = %f\n",
            pLnp_par->SongSTDP_param.STDPmultiplier);

        printf("antiCausalCoeff = %f\n",
            pLnp_par->SongSTDP_param.antiCausalCoeff);

        printf("tauSTDP_ms = %f\n",
	    pLnp_par->SongSTDP_param.tauSTDP_ms);

        printf("causalSTDPmaxAbsValue = %f\n",
            pLnp_par->SongSTDP_param.causalSTDPmaxAbsValue);

        printf("antiCausalSTDPmaxAbsValue = %f\n",
            pLnp_par->SongSTDP_param.antiCausalSTDPmaxAbsValue);

        printf("tauDerivativeDecay_s = %f\n",
            pLnp_par->SongSTDP_param.tauDerivativeDecay_s);

        printf("derivativeDrift = %f\n",
            pLnp_par->SongSTDP_param.derivativeDrift);
      }
    DPSNNverboseEnd();
    return;
  }
  #undef maxNumOfPlasticityParam
  #undef bufferSize
}
