// DPSNN_localNet_init.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy Computation with Spikes" 

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
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


void localNetClass::prepareForwardSynapses()
{

  if((lnp_par.synGen == default_random_synGen_1)||
     (lnp_par.synGen == simpleCorticalModule_synGen_2)||
     (lnp_par.synGen == randTable_simpleCorticalModule_synGen_3)||
     (lnp_par.synGen == LIFCACorticalModule_synGen_4)) {
    simpleCM_prepareForwardSynapses();
  }else { 
    printf("ERROR unrecognized syn Gen option\n");fflush(stdout);exit(0);
  };
};

void localNetClass::init(
  const struct DPSNN_parameters lnp_par_initValue, 
  messagePassingClass * pMessagePassing_initValue,
  stopWatchClass * pStopWatch_initValue,
  statClass *pStat_initValue) 
{
    // initializes a network of locN neurons stored and managed by the host named loc_h
    uint32_t i;//,j,k,r,exists;
    lnp_par = lnp_par_initValue;
 
    struct sysinfo si;

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
      printf(
      "-100- init() loc_h=%d start,locN=%d,M=%d,H=%d,CFT=%d,CFX=%d,CFY=%d\n", 
      lnp_par.loc_h, lnp_par.locN, lnp_par.M, 
      lnp_par.globH, lnp_par.globCFT, lnp_par.globCFX, lnp_par.globCFY);
    DPSNNverboseEnd();

    DPSNNverboseStart(true,1,0);    

      //freeMemory=system ("free");
      sysinfo (&si);

      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();
      lnp_par.freeMemAtStart = si.freeram/1024;

      if(lnp_par.loc_h==0){
        //freeMemory=system ("free");
        printf ("MEMORY: At init start the free memory on %s is %d kB \n",
	      lnp_par.hostName, lnp_par.freeMemAtStart);
      }
    DPSNNverboseEnd();


    if(DSD__maxBackwardLocSyn != DSD__maxForwardLocSyn) {
      printf(
     "ERROR: after mem optim DSD__maxBackwardLocSyn should be equal to forward\n");
      fflush(stdout);exit(0);
    }

    memPoolA = (uint8_t *)new synapseClass [DSD__maxBackwardLocSyn];
    memPoolB = (uint8_t *)new synapseClass [DSD__maxBackwardLocSyn];
 
    DPSNNverboseStart(false,1,0);    
      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();

      if(lnp_par.loc_h==0){
        //freeMemory=system ("free");
        sysinfo (&si);
        printf ("MEMORY: AFTER memPoolA and B new the used memory on %s is %lu kB \n",
	      lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
    DPSNNverboseEnd();

    pMessagePassing = pMessagePassing_initValue;
    pStopWatch = pStopWatch_initValue;
    pStat = pStat_initValue;
    
    localSynCount=0;
    forwardSynCount=0;
    backwardSynCount=0;
    reportCount=0;

    N_firingsTotInFrame = 0;
    clearFiringsInChronoWindow();
    for(i=0;i<lnp_par.globCFT*lnp_par.subPopNumber;i++)
      N_firingsPerPop[i] = 0;
  
    DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h==0)
      printf ("Random seed for simulation: %d \n",lnp_par.globalSeed);
    DPSNNverboseEnd();

    localNetRandDev.SetRandomSeed(lnp_par.globalSeed);
    simpleCM_connectome.initLocalNetRandDevPointer(&localNetRandDev);
 
    DPSNNverboseStart(false,0,0);
    if(lnp_par.loc_h==0) {
      memMeasure.clear();
      memMeasure.measure(lnp_par);
      memMeasure.report(lnp_par);
    };
    DPSNNverboseEnd();

    clearAllChronometers();
    initChronoWindow();

    check_lnp_par_initValues();

    localSynList = (synapseClass*)&memPoolA[0];
    forwardSynList = (synapseClass*)&memPoolB[0];

    DPSNNverboseStart(false,1,0);    
      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();

      if(lnp_par.loc_h == 0){
        //freeMemory=system ("free");
        sysinfo (&si);
        printf ("MEMORY: BEFORE prepareForward Synapses the used memory on %s is %lu kB \n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
    DPSNNverboseEnd();

    prepareForwardSynapses();

    DPSNNverboseStart(false,1,0);    
      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();

      if(lnp_par.loc_h == 0){
        //freeMemory=system ("free");
        sysinfo (&si);
        printf ("MEMORY: AFTER prepareForward Synapses the used memory on %s is %lu kB \n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
    DPSNNverboseEnd();

    DPSNNverboseStart(false,1,0);    
      printf("-400- init()  h=%03d globH=%d  locN=%d\n", 
	   lnp_par.loc_h, lnp_par.globH, lnp_par.locN);
    DPSNNverboseEnd();

    for(i=0;i<lnp_par.subPopNumber;i++)
      neuSubPopCount[i] = neuSubPopParam[i].count;

    DPSNNverboseStart(false,1,0);
	reportLocalNetAfterInit();
    DPSNNverboseEnd();

    DPSNNverboseStart(true,1,0); 
      chronoSendReceiveSynList.clearAndStartChrono();
    DPSNNverboseEnd();

    sendReceiveSynListWithOtherLocalNets();
    //NOTE: backwardSynList is in memPoolA
    //memPoolB can be deallocated
    delete [] memPoolB;

    DPSNNverboseStart(false,1,0);    
      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();

      if(lnp_par.loc_h == 0){
        //freeMemory=system ("free");
        sysinfo (&si);
        printf ("MEMORY: AFTER delete memPoolB free on %s is %lu kB \n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
    DPSNNverboseEnd();

    DPSNNverboseStart(true,1,0); 
      chronoSendReceiveSynList.stopChrono();
    DPSNNverboseEnd();
    DPSNNverboseStart(false,1,0);    
      printf("-600- init()  h=%03d globH=%d  locN=%d\n", 
	   lnp_par.loc_h, lnp_par.globH, lnp_par.locN);
    DPSNNverboseEnd();


#if defined(makeActiveLTD) || defined (makeActiveLTP)
    DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h==0 && i==0) {
      printf("localNet_init: initializing Plasticity on h=0, local neuron 0\n)");fflush(stdout);
    }
    DPSNNverboseEnd();
    for (i=0;i<lnp_par.locN;i++) {
      //clear the counter of incoming syapses before adding the list of requested incoming syapses
      n[i].clearBackwardSynList();   
     //we have to put the last spike eons ago
      //otherwise at the sim start (t=0) LTP LTD and activity would be wrong
      n[i].setLastEmittedSpikeTime_ms(-100000);
      n[i].initLongTermPlasticity();
    };

    for (i=0;i<lnp_par.locN;i++) {
     n[i].initBackwardSynListPointer(backwardSynList);
    };

    {/* adding the incoming synapses */
      uint32_t synOffset, i;
      synapseClass backwardSyn;
      for(synOffset=0;synOffset<backwardSynCount;synOffset++) {
	backwardSyn=backwardSynList[synOffset];
	i = backwardSyn.post_glob_n % lnp_par.locN;
	if((backwardSyn.post_glob_n / lnp_par.locN) != lnp_par.loc_h) {
	  printf("ERROR wrong post_glob_n when adding incoming synapses\n");
	  fflush(stdout);exit(0);
	}
        n[i].addBackwardSyn(synOffset);
      };
    };

#endif


#ifndef LIFCAneuron
    for (i=0;i<lnp_par.locN;i++)
      n[i].set_hashId();
#endif

    DPSNNverboseStart(false,1,0);    
      writeIniFiles();
    DPSNNverboseEnd();

    DPSNNverboseStart(false,1,0);    
      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();

      if(lnp_par.loc_h == 0){
        //freeMemory=system ("free");
        sysinfo (&si);
        printf ("MEMORY: BEFORE axonalSpikeScheduler.reset() free on %s is %lu kB \n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
    DPSNNverboseEnd();

    prepareAxonalSpikeBuffers();

    //EPA - 2015-03-10
    //Allocate memory for inputCurrents buffer of each neuron
    for(i=0;i<lnp_par.locN;i++) {
      n[i].inputCurrents = (inputCurrentClass *)new inputCurrentClass [DSD__maxSimultaneousSpikesOnSameTarget];
    }

    for (i=0;i<lnp_par.locN;i++) {
      n[i].initRandDevPointer(&localNetRandDev);
    };

    DPSNNverboseStart(false,1,0);    
    {
      uint64_t counterRandGenerated;
      counterRandGenerated = localNetRandDev.Statistics();
      printf("On h=%d after init %lu random numbers have been generated \n",lnp_par.loc_h,counterRandGenerated);
    }
    DPSNNverboseEnd();

    localNetRandDev.SetRandomSeed(lnp_par.loc_h);

    DPSNNverboseStart(false,1,0);    
      printf("-700- init()  h=%03d globH=%d  locN=%d\n", 
	   lnp_par.loc_h, lnp_par.globH, lnp_par.locN);
    DPSNNverboseEnd();

    DPSNNverboseStart(true,1,0);    
      //MPI_Barrier(MPI_COMM_WORLD);
      pMessagePassing->barrier();

      if(lnp_par.loc_h == 0){
        //freeMemory=system ("free");
        sysinfo (&si);
        printf ("MEMORY: End of init phase the used memory on %s is %lu kB \n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
    DPSNNverboseEnd();

    #if defined(makeActiveLTD) || defined (makeActiveLTP)
    DPSNNverboseStart(false,1,0);
      checkBackwardSynCountAfterInit();
    DPSNNverboseEnd();
    #endif
  };

//-------END OF INIT METHOD ----------------//

void localNetClass::sendReceiveSynListWithOtherLocalNets()
{//BEGIN send list of synaptic connections to different hosts 
    uint32_t i;
    uint32_t backwardRemoteSynCount;
    struct sysinfo si;

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
	printf("-510- h=%03d sendReceiveSynList: START\n", lnp_par.loc_h);
    DPSNNverboseEnd();

    DPSNNverboseStart(false,1,0);  
      printf("520- h=%03d sendReceiveSynList: BEFORE check BEFORE sort\n",
       lnp_par.loc_h);
    DPSNNverboseEnd();

    DPSNNverboseStart(false,1,0);  
      printf("521- h=%03d localSynCount=%d, forwardSynCount=%d, totSynNum=%d\n",
	     lnp_par.loc_h, localSynCount, forwardSynCount, projectedSynCount);
    DPSNNverboseEnd();

    if(localSynCount + forwardSynCount != projectedSynCount ) {
      printf(
    "ERROR h=%03d localNet_init localSynCount %d + forwardSynCount %d != projectedSynCount %d\n",
	     lnp_par.loc_h, localSynCount, 
	     forwardSynCount, projectedSynCount ); 
	fflush(stdout); exit(0);
    };

    //check list (only thiose toward remote processes) before sort 
    DPSNNverboseStart(true,1,0); 
      checkForwardSynListInitValues(forwardSynList, forwardSynCount);
    DPSNNverboseEnd();

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
      printf("530- h=%03d sendReceiveSynList: AFTER check BEFORE sort\n",
       lnp_par.loc_h);
    DPSNNverboseEnd();

    DPSNNverboseStart(true,1,0); 
      chronoSortTargetHostInForwardSynList.clearAndStartChrono();
    DPSNNverboseEnd();

    //HERE THE SORT of forwardSynList according to targetHost
    //we put q1 in the second part of memPoolB, the first part is used
    //by forwardSynList
    //this is only about synapses toward remote processes
    q1 = (synapseClass *)&memPoolA[sizeof(synapseClass) * localSynCount];
    sortTargetHostInForwardSynList(forwardSynList, q1, forwardSynCount);

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
      printf("539- h=%03d sendReceiveSynList: AFTER sortTargetHostInForwardSynList BEFORE sort\n",
       lnp_par.loc_h);
    DPSNNverboseEnd();

    //this is only about synapses toward remote processes
    checkSortedSynListAndCreateSynapticDistribution(
      forwardSynList, &synTargetHostDistribution, forwardSynCount);
    // CHECK the sum of synHostDistribution.synCount[x] 
    // must be equal to synCount

    //this is only about synapses toward remote processes
    checkTotalSynCountInSynapticDistribution(&synTargetHostDistribution, 
					     forwardSynCount);

    DPSNNverboseStart(true,1,0); 
      chronoSortTargetHostInForwardSynList.stopChrono();
    DPSNNverboseEnd();

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
      //warning this report is only about synapses toward remote processes
      forwardSynListReport(forwardSynList, forwardSynCount,0);
    DPSNNverboseEnd();

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
      printf("540- h=%03d sendReceiveSynList: AFTER sort BEFORE CHECK\n",
       lnp_par.loc_h);
    DPSNNverboseEnd();

    //check list after sort
    DPSNNverboseStart(false,1,0); 
      checkForwardSynListInitValues(forwardSynList, forwardSynCount);
    DPSNNverboseEnd();

    //debug report
    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
	printf("560- h=%03d is going to call pMessagePassing->init()\n", lnp_par.loc_h);
    DPSNNverboseEnd();

    //print distribution of forward synapses 
    //(from each source to each target)
    DPSNNverboseStart(false,16,lnp_par.loc_h);  
    for(i=0;i<lnp_par.globH;i++) {
      uint32_t synCount;
      if(i!=lnp_par.loc_h) {
	synCount = synTargetHostDistribution.synCount[i];
      }else{
	synCount = localSynCount;
      };
      printf("570- S_h=%03d -> T_h=%03d - %d syn, forwOff in remoteTargeDistrib %d\n",
	     lnp_par.loc_h, i, synCount, synTargetHostDistribution.synOffset[i] ); 
	fflush(stdout);
    };
    DPSNNverboseEnd();
    
    //check of TargetHostDistribution before message passing
    DPSNNverboseStart(false,1,0); 
    {
      uint32_t checkGlobalForwardSynCount;
      checkGlobalForwardSynCount=0;
      for(i=0;i<lnp_par.globH;i++) {
	checkGlobalForwardSynCount += synTargetHostDistribution.synCount[i];
      }
      checkGlobalForwardSynCount += localSynCount;

      if((checkGlobalForwardSynCount >= DSD__maxForwardLocSyn) ||
	 (checkGlobalForwardSynCount != (forwardSynCount + localSynCount))) {
	  printf(
       "ERROR h=%03d in localNet::init() error checking globalForwardSynCount=%d\n",
	    lnp_par.loc_h,checkGlobalForwardSynCount); 
	    fflush(stdout);exit(0);
	};
      if((synTargetHostDistribution.synOffset[lnp_par.globH-1] +
	  synTargetHostDistribution.synCount[lnp_par.globH-1] + localSynCount) !=
	 (forwardSynCount + localSynCount)) 
      {
	  printf(
	    "ERROR h=%03d in localNet::init() error checking globalForwardSynOffset o%d c%d l%d f%d \n",
	    lnp_par.loc_h,
	    synTargetHostDistribution.synOffset[lnp_par.globH-1],
	    synTargetHostDistribution.synCount[lnp_par.globH-1],
	    localSynCount,
	    forwardSynCount); 
	    fflush(stdout);
	    exit(0);
      };
    }
    DPSNNverboseEnd();

    //HERE IS THE MESSAGE PASSING OF THE EXPECTED DISTRIBUTION OF SYNAPSES
    //NOTE: synapses directed toward same process are non considered here
    pMessagePassing->sendForwardSynListDimToRemoteHosts(
      &synTargetHostDistribution, 
      &synSourceHostDistribution,0);

    if((synTargetHostDistribution.synCount[lnp_par.loc_h]!=0 )||
       (synSourceHostDistribution.synCount[lnp_par.loc_h]!=0 )){ 
      printf(
    "ERROR h=%03d in localNet::init() erroneous local transfer of %d to %d syn\n",
	      lnp_par.loc_h,
	      synTargetHostDistribution.synCount[lnp_par.loc_h],
	      synSourceHostDistribution.synCount[lnp_par.loc_h]);
      fflush(stdout);exit(0);
    }

    backwardRemoteSynCount=
      synSourceHostDistribution.synOffset[lnp_par.globH-1] + 
      synSourceHostDistribution.synCount[lnp_par.globH-1];

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
      printf("300- h=%03d - MPI_Alltoall told %d total remote syn will arrive\n",
        lnp_par.loc_h, backwardRemoteSynCount);
    DPSNNverboseEnd();

    //EPA PSP this check should be more restricted e.g. to DSD__maxRemoteM
    if(backwardRemoteSynCount >= DSD__maxBackwardLocSyn){
      printf(
      "ERROR on loc_h=%d in localNet::init() backwardSynCount > DSD__maxBackwardLocSyn\n",lnp_par.loc_h); 
      fflush(stdout);exit(0);
    };

    DPSNNverboseStart(false,4,lnp_par.loc_h);
    for(i=0;i<lnp_par.globH;i++) {
      uint32_t synCount;
      if(i!=lnp_par.loc_h) {
	synCount = synSourceHostDistribution.synCount[i];
      }else{
	synCount = localSynCount;
      };
      printf("580- S_h=%03d -> T_h=%03d - %d back syn after dim message passing\n",
	     i, lnp_par.loc_h, synCount);
      fflush(stdout);
    }
    DPSNNverboseEnd();

    {

      //intermSynList is placed in the second part of memPoolA
      //the first part of memPoolA is used by localSynList
      intermSynList=(synapseClass *) &memPoolA[
		    sizeof (synapseClass) * localSynCount];
 
      DPSNNverboseStart(false,1,0); 
      //MPI_Barrier(MPI_COMM_WORLD);
        pMessagePassing->barrier();
 
        if(lnp_par.loc_h==0){
	  //freeMemory=system ("free");
	  sysinfo (&si);
	  printf ("MEMORY: BEFORE send/rec the used memory on %s is %lu kB\n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
        }
      DPSNNverboseEnd();

      pMessagePassing->sendForwardSynListToRemoteHosts(
			       forwardSynList, &synTargetHostDistribution,
			       intermSynList, &synSourceHostDistribution, 0);

      DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
      printf("590- h=%03d sendReceiveSynList: BEFORE check AFTER MESSAGE PASS\n",
	     lnp_par.loc_h);
      DPSNNverboseEnd();

      //before the following line backwardSynCount, was about obly the 
      //synapses transferred by message passing
      backwardSynCount = backwardRemoteSynCount + localSynCount;
      //now backward syn count is about the total pool of incoming synapses

      //forwardSynList (allocated on memPoolB is not anymore necessary
      //from now on memPoolB is used for backwardSynList
      backwardSynList=(synapseClass*)&memPoolB[0];
      
      {//inserting localSynList in the middle of the full list
	//rebuilding a single full backwardSynList
	//EPA - PSP correggere qui almeno il count delle forward syn
	uint32_t h,s;
	uint32_t backwardSynOffset;
	uint32_t intermSynOffset;

        DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
          printf("591- h=%03d localNet_init start inclusion local in middle\n",
	     lnp_par.loc_h);
        DPSNNverboseEnd();

	backwardSynOffset=0; intermSynOffset=0;
	synSourceHostDistribution.synOffset[0]=0;

        for(h=0;h<lnp_par.loc_h;h++) {
	  for(s=0;s<synSourceHostDistribution.synCount[h];s++) {
	    backwardSynList[backwardSynOffset++] = 
	      intermSynList[intermSynOffset++];
	  };
	};

	for(s=0;s<localSynCount;s++) {
	  backwardSynList[backwardSynOffset++] = localSynList[s];
	};

	synSourceHostDistribution.synCount[lnp_par.loc_h]=localSynCount;

	if(lnp_par.loc_h > 0) {
	  synSourceHostDistribution.synOffset[lnp_par.loc_h]=
	    synSourceHostDistribution.synOffset[lnp_par.loc_h-1]+
	    synSourceHostDistribution.synCount[lnp_par.loc_h-1];
	};

        for(h=lnp_par.loc_h+1;h<lnp_par.globH;h++) {
	  for(s=0;s<synSourceHostDistribution.synCount[h];s++) {
	    backwardSynList[backwardSynOffset++] = 
	      intermSynList[intermSynOffset++];
	  };
	  synSourceHostDistribution.synOffset[h]=
	    synSourceHostDistribution.synOffset[h-1]+
	    synSourceHostDistribution.synCount[h-1];
	};

        if(backwardSynCount!=backwardSynOffset) {
	  printf(
          "ERROR h=%d localNet_init backSynCount=%d != backSynOffset=%d\n", 
	 lnp_par.loc_h, backwardSynCount, backwardSynOffset); 
	  fflush(stdout); exit(0); 
	}
        DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);  
          printf("592- h=%03d localNet_init END inclusion local in middle\n",
	     lnp_par.loc_h);
	  DPSNNverboseEnd();
      };

      // memPoolB contains the full backwardSynList (both local and remote) 
      // reordered per source_h

      inflateIntermToBackwardSynList(
	(synapseClass*) &memPoolB[0],
        (synapseClass*) &memPoolA[0], backwardSynCount);

      backwardSynList=(synapseClass*)&memPoolA[0];

      DPSNNverboseStart(true,1,0); 
        //MPI_Barrier(MPI_COMM_WORLD);
        pMessagePassing->barrier();
        if(lnp_par.loc_h==0){
	  //freeMemory=system ("free");
	  sysinfo (&si);
  	  printf ("MEMORY: AFTER inflate the used memory on %s is %lu kB\n",
		lnp_par.hostName, lnp_par.freeMemAtStart-si.freeram/1024);
      }
      DPSNNverboseEnd();
    }

    DPSNNverboseStart(false,1,0);
           printf("593- h=%03d localSynCount=%d, backwardSynCount=%d, totSynNum=%d\n",
	     lnp_par.loc_h, localSynCount, backwardRemoteSynCount, backwardSynCount);
	   fflush(stdout);
    DPSNNverboseEnd();

    DPSNNverboseStart(true,1,0);
      synapticListReport(backwardSynList, backwardSynCount, -1);
    DPSNNverboseEnd();

    DPSNNverboseStart(false,1,0);
      checkSynListInitValues(backwardSynList, backwardSynCount);
    DPSNNverboseEnd();
 
    //The delivery of axonal spikes requires to know if there are synapses directed
    //roward the same process
    //the offset should be not ceseccary, because the fowardSynList
    //has been overwritten by previous memPools management 
    synTargetHostDistribution.synCount[lnp_par.loc_h]=localSynCount;
   
  };

  void localNetClass::inflateIntermToBackwardSynList(
    const synapseClass *intermSynList, 
    synapseClass *backwardSynList, const uint32_t backwardSynCount) {
    uint32_t i,j;
    uint32_t source_h, delay;
    uint32_t delayCntExpected[lnp_par.D-1];
    uint32_t delayCntActual[lnp_par.D-1];
    uint32_t offsetd[lnp_par.D-1];
    uint64_t neuId;
    FILE *fp_output;
    char reportName[80];
    //    uint32_t seedForLIFCAWeightInit;
    //struct sysinfo si;
    uint32_t synSourceCount;

    for (i=0;i<lnp_par.D;i++){
      delayCntExpected[i]=0;
      delayCntActual[i]=0;
    }
 
    for (i=0;i<backwardSynCount;i++)
      delayCntExpected[intermSynList[i].delay]++;

    offsetd[0]=0;
    for (i=1;i<lnp_par.D;i++)
      offsetd[i]=offsetd[i-1]+delayCntExpected[i-1];

    for(i=0;i<backwardSynCount;i++){

      j=offsetd[intermSynList[i].delay]+delayCntActual[intermSynList[i].delay];
      delayCntActual[intermSynList[i].delay]++;
      backwardSynList[j].pre_glob_n = 
	intermSynList[i].pre_glob_n ; 
      backwardSynList[j].post_glob_n = 
	intermSynList[i].post_glob_n ; 
      backwardSynList[j].delay = 
	intermSynList[i].delay;
      backwardSynList[j].preSynNeuralKind = 
	intermSynList[i].preSynNeuralKind;
      backwardSynList[j].weight = 
	intermSynList[i].weight;
 
#if defined(makeActiveLTD) || defined (makeActiveLTP)
      backwardSynList[j].lastActivationTime = 0;
      backwardSynList[j].timeDerivative = (weightType)0.0;
#endif
      delay = backwardSynList[j].delay;
      source_h = backwardSynList[j].pre_glob_n / lnp_par.locN;
      synDelaySourceHostDistribution.synCount[delay][source_h]++;
    }

    for(delay=0;delay<lnp_par.D;delay++){
      synDelaySourceHostDistribution.synOffset[delay][0]=offsetd[delay];
      for(source_h=1;source_h<lnp_par.globH;source_h++)
	synDelaySourceHostDistribution.synOffset[delay][source_h]=
	  synDelaySourceHostDistribution.synOffset[delay][source_h-1]+
	  synDelaySourceHostDistribution.synCount[delay][source_h-1];
    }

    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
    for (i=0;i<lnp_par.D;i++)
      printf("BackwardSynList h=%d: delayCntExpected[%d]=%d, delayCntActual[%d]=%d, offsetd[%d]=%d\n",
	     lnp_par.loc_h,i,delayCntExpected[i],i,delayCntActual[i],i,offsetd[i]);
    DPSNNverboseEnd();

    for (i=0;i<lnp_par.D;i++){
      if(delayCntExpected[i]!=delayCntActual[i])
	printf("ERROR in BackwardSynList: mismatch between expected=%d and actual=%d synapses with delay %d in process h=%d \n",
	       delayCntExpected[i],delayCntActual[i],i,lnp_par.loc_h);
    }

    for (j=0;j<lnp_par.locN*lnp_par.globH;j++){
      synSourceCount = synSourceHostDistribution.synCount[j/lnp_par.locN];
      if(synSourceCount != 0) {
	for (i=0;i<lnp_par.D;i++)
	  backwardSynOffsetInSynList[j].offsetByDelay[i] = -1;
      }
    }
	
    // Calculate synapse offsets in the reordered backward synaptic list
    // for each synapse: 1 offset for each delay
    neuId = (uint64_t)-1;
    for (i=0;i<lnp_par.D;i++){    
      for (j=0;j<delayCntActual[i];j++){
	if(backwardSynList[j+offsetd[i]].pre_glob_n != neuId){
	  neuId = backwardSynList[j+offsetd[i]].pre_glob_n;
	  synSourceCount = synSourceHostDistribution.synCount[neuId/lnp_par.locN];
	  if(synSourceCount != 0)
	    backwardSynOffsetInSynList[neuId].offsetByDelay[i] = offsetd[i] + j;
	}
      }
      neuId = (uint64_t)-1;
      }
    
    DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
    sprintf(reportName,"OffsetSynList_h%d.txt",lnp_par.loc_h);
    fp_output=fopen(reportName,"w");
    fprintf(fp_output,"Offsets by delay for process H=%d\n",lnp_par.loc_h);     
      for(i=0;i<lnp_par.locN*lnp_par.globH;i++){
	fprintf(fp_output,
		"neuId=%d offset[0]=%d offset[1]=%d offset[2]=%d offset[3]=%d offset[4]=%d \n",i,
		backwardSynOffsetInSynList[i].offsetByDelay[0],
		backwardSynOffsetInSynList[i].offsetByDelay[1],
		backwardSynOffsetInSynList[i].offsetByDelay[2],
		backwardSynOffsetInSynList[i].offsetByDelay[3],
		backwardSynOffsetInSynList[i].offsetByDelay[4]);
      }
    DPSNNverboseEnd();
    
  };

  void localNetClass::check_lnp_par_initValues() {
 
    if(lnp_par.D>DSD__maxD) {
      printf("ERROR: localNet::init - maximum delay D=%d out of range\n", lnp_par.D);
      fflush(stdout);exit(0);}
    if(lnp_par.locN>DSD__maxLocN) {
      printf("ERROR: localNet::init - locN=%d wrong value out of range\n", lnp_par.locN);
      fflush(stdout);exit(0);}
    if(lnp_par.loc_h>=DSD__maxGlobH||lnp_par.loc_h>=lnp_par.globH) {
      printf("ERROR: localNet::init - loc_h=%d wrong value out of range\n", lnp_par.loc_h);
      fflush(stdout);exit(0);}
  };

  void localNetClass::forwardSynListReport( 
     synapseClass *synList, 
     uint32_t synCount, uint32_t repNum) {
    FILE *fp_output;
    char reportName[80];
    instrumentedSynapse instrumentedSynapseDummy;
    uint32_t i;
    sprintf(reportName,"localNet_h%d_N%d_SynapticList_FORWARD_R%d.dat",
	    lnp_par.loc_h,lnp_par.locN,repNum);
    fp_output=fopen(reportName,"w");
    fprintf(fp_output,
      "localNet::statusReport - list in localNet::forwardSynList \n");
    for (i=0;i<synCount;i++) {
      instrumentedSynapseDummy.report(synList[i],fp_output,
        lnp_par.locN,lnp_par.factorWeightType_2_Float);
    }
    fclose(fp_output);
  };

 void localNetClass::synapticListReport(synapseClass *synList, 
					uint32_t synCount, uint32_t repNum) {
    FILE *fp_output;
    char reportName[80];
    instrumentedSynapse instrumentedSynapseDummy;
    uint32_t i;
    DPSNNverboseStart(false,1,0);
    if(lnp_par.loc_h==0) {
      printf("loc_h = %d, synapticListReport starting\n", lnp_par.loc_h);
      fflush(stdout); }
    DPSNNverboseEnd();

    if(lnp_par.loc_h==0) {
      sprintf(reportName,"localNet_h%d_N%d_SynapticList_BACKWARD_R%d.dat",
	    lnp_par.loc_h,lnp_par.locN,repNum);
      fp_output=fopen(reportName,"w");
      //fprintf(fp_output,
      //"localNet::statusReport - list in localNet::SynapticList \n");
      for (i=0;i<synCount;i++) {
        instrumentedSynapseDummy.report(synList[i],fp_output,
	  lnp_par.locN,lnp_par.factorWeightType_2_Float);
      }
      fclose(fp_output);
    };
    DPSNNverboseStart(false,1,0);
    if(lnp_par.loc_h==0) {
      printf("loc_h = %d, synapticListReport ending\n", lnp_par.loc_h);
      fflush(stdout); }
    DPSNNverboseEnd();
  };

  void localNetClass::checkTotalSynCountInSynapticDistribution(
	    synapticDistributionClass *pSynapticDistribution, 
	    uint32_t expectedSynCount) {	
    // CHECK the sum of synTargetHostDistribution.synCount[x] 
    // must be equal to synCount 
    uint32_t synCountCheck, i;
    synCountCheck=0;
    for (i=0;i<lnp_par.globH;i++) {
      synCountCheck+=pSynapticDistribution->synCount[i];
      DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
	  printf(
     "checkTotalSynCount on loc_h=%d - synCount=%d syn toward host %d starts at synOffset=%d\n", 
		 lnp_par.loc_h, 
		 pSynapticDistribution->synCount[i], 
		 i, pSynapticDistribution->synOffset[i]);
      DPSNNverboseEnd();
    };
    if(expectedSynCount!=synCountCheck) {
      printf(
	"ERROR total =%d syn toward all host targets does not amount to expected %d, \n",
	synCountCheck, expectedSynCount);
      fflush(stdout);
      exit(0);
    };
  };

void localNetClass::checkSortedSynListAndCreateSynapticDistribution(
	    synapseClass synList[],  
            synapticDistributionClass *pSynHostDistribution, 
            uint32_t synCount)	
{
  //MUST BE EXECUTED AFTER THE SORT OF THE SYN LIST BY TARGET HOST
 
  if(lnp_par.globH==1) {
      pSynHostDistribution->synCount[0]=synCount;
      pSynHostDistribution->synOffset[0]=0;
  }else{
    //counts the number of synapses connected to each host target and verifies the sorting

    //clear the bins in the target host distribution histogram
    uint32_t currentTargetHost, i;
    for (i=0;i<lnp_par.globH;i++) {
      pSynHostDistribution->synCount[i]=0;
      pSynHostDistribution->synOffset[i]=0;
    };

    //the first column in the histogram 
    //under scrutiny is obviously the #0
    currentTargetHost=0;

    //i is the counter of synapses in the list (already sorted by targetHost) 
    for (i=0;i<synCount;i++) {
      //for each synapse there are 3 possibilities:
      //A- the synapse has the same target host of the previous synapse
      //   (increase the count in the current histogramm bin)
      //B- the synapse is in the next bin (NO HOLES)
      //   (increase the histogram bin counter, set to 1 the counter
      //C- there are holes in the histogram
      //   (skip the bins, setting to zero the height and copying the offset)
      
      //case A - same histogram bin of previous synapse
      if((synList[i].post_glob_n / lnp_par.locN) == currentTargetHost)  {
	pSynHostDistribution->synCount[currentTargetHost]++;
	DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
  	  printf(
          "synDistr A - h=%d: another synapse %d added to total %d of synapses directed to host %d\n", 
          lnp_par.loc_h, i, 
          pSynHostDistribution->synCount[currentTargetHost], 
          currentTargetHost);
        DPSNNverboseEnd();
      } else if((synList[i].post_glob_n / lnp_par.locN)==(currentTargetHost+1)) {
	//case B - next histogram bin
	// (in this case no holes in the distribution) 
	currentTargetHost++;
	pSynHostDistribution->synOffset[currentTargetHost] = i;
	pSynHostDistribution->synCount[currentTargetHost]  = 1;
	DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
  	  printf(
          "synDistr B - h=%d: first synapse %d added to total %d of synapses directed to host %d\n", 
          lnp_par.loc_h, i, 
	  pSynHostDistribution->synCount[currentTargetHost], 
	  currentTargetHost);
       DPSNNverboseEnd();
      } else if((synList[i].post_glob_n / lnp_par.locN)>(currentTargetHost+1)) { 
	//case C - set empty columns in the histogram 
	// in this case there are holes in the distribution:
	//i.e. this source host does not project synapses 
	//neither to currentTargetHost (case A), 
	//nor to currentTargetHost+1 (case B)
	do {
	  currentTargetHost++;
	  pSynHostDistribution->synCount[currentTargetHost]=0;
	  //the offset does not move
	  pSynHostDistribution->synOffset[currentTargetHost]=
	    pSynHostDistribution->synOffset[currentTargetHost-1] + 
	    pSynHostDistribution->synCount[currentTargetHost-1];
	  DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
 	      printf("synDistr C - h=%d: skip target host %d\n",
		     lnp_par.loc_h, currentTargetHost);
	      fflush(stdout);
          DPSNNverboseEnd();
  	} while ((currentTargetHost+1)!=(synList[i].post_glob_n / lnp_par.locN));
	//now we should have skipped all the empty bins
	currentTargetHost++; //this should be good... 
	if(currentTargetHost == (synList[i].post_glob_n / lnp_par.locN)) {
	  // ...and now set the offset and count of the first synapse 
	  // directed to this existing target
	    pSynHostDistribution->synOffset[currentTargetHost]=i;
	    pSynHostDistribution->synCount[currentTargetHost]++;
	} else {
	  if (currentTargetHost>=lnp_par.globH) { 
	    printf("ERROR currentTargetHost error - 1\n");fflush(stdout);exit(0);
	  };
	}; 	        
      } else if ((synList[i].post_glob_n / lnp_par.locN)<currentTargetHost) { 
	// we are after the sort so we must not go backwards 
	printf("ERROR synSort not performed\n"); fflush(stdout); exit(0);
      } else if((synList[i].post_glob_n / lnp_par.locN) >= lnp_par.globH) {
	// the target host of the synapse must be in range 
	printf("ERROR: target host too big\n"); fflush(stdout); exit(0);
      };
      if (currentTargetHost>=lnp_par.globH) {
	// the current bin of the histogram must be in range
        printf("ERROR currentTargetHost error -2 \n");fflush(stdout);exit(0);
      };
    };//end of for i over synapses
    //there could be processes not reached by any of the synapses
    //one or more "holes" at the end of h ordering
    { 
      //setting the a valid offset for all final targets processes 
      //not receiving synapses
      uint32_t h;
      for(h = currentTargetHost+1; h < lnp_par.globH ; h ++)
	pSynHostDistribution->synOffset[h]=
	  pSynHostDistribution->synOffset[h-1] +
	  pSynHostDistribution->synCount[h-1];
    }
  };//end of if about globH==1 (single process)
};


  void localNetClass::checkForwardSynListInitValues(
            synapseClass *synList,
            uint32_t synCount) 
  {
    uint32_t i;
    if(synCount >= DSD__maxForwardLocSyn){
	printf("ERROR on loc_h=%d - checkForwardSynListInitValues - too many %d forward syn generated\n",
	       lnp_par.loc_h,synCount);fflush(stdout);exit(0);
    }
    for (i=0;i<synCount;i++) {
      if (synList[i].pre_glob_n>=lnp_par.globN) 
	{printf("ERROR on loc_h=%d - checkForwardSynListInitValues - pre_glob_n\n",
		lnp_par.loc_h);fflush(stdout);exit(0);};
      if (synList[i].post_glob_n>=lnp_par.globN) 
	{printf("ERROR on loc_h=%d - checkForwardSynListInitValues - post_glob_n\n",
		lnp_par.loc_h);fflush(stdout);exit(0);};
      if (synList[i].delay > lnp_par.D) 
	{printf("ERROR on loc_h=%d - checkForwardSynListInitValues - delay (synList[%d].delay=%d)\n",
		lnp_par.loc_h,i,synList[i].delay);fflush(stdout);exit(0);};	
#ifdef LIFCAneuron
      /*PSP 2017-plasticity      if (!((synList[i].preSynNeuralKind == excitatoryLbExc) ||
	    (synList[i].preSynNeuralKind == excitatoryLaExc) ||
            (synList[i].preSynNeuralKind == inhibitoryLaInh)))
	{printf("ERROR on loc_h=%d -chkSynListInitVal (1) - preSynNeuralKind\n",
	lnp_par.loc_h); fflush(stdout); exit(0);};*/
#else
      if (!((synList[i].preSynNeuralKind == excitatoryRS) ||
            (synList[i].preSynNeuralKind == inhibitoryFS)))
	{printf("ERROR on loc_h=%d -chkSynListInitVal - preSynNeuralKind \n",
		lnp_par.loc_h); fflush(stdout); exit(0);};
#endif
    };
  };

 void localNetClass::checkSynListInitValues(
            synapseClass *synList,
            uint32_t synCount) 
{
   uint32_t i;
   for (i=0;i<synCount;i++) {
      if (synList[i].pre_glob_n>=lnp_par.globN) 
	{printf("ERROR on loc_h=%d - checkSynListInitValues - pre_glob_n\n",
		lnp_par.loc_h);fflush(stdout);exit(0);};
      if (synList[i].post_glob_n>=lnp_par.globN) 
	{printf("ERROR on loc_h=%d - checkSynListInitValues - post_glob_n\n",
		lnp_par.loc_h);fflush(stdout);exit(0);};
      if (synList[i].delay > lnp_par.D) 
	{printf("ERROR on loc_h=%d - checkSynListInitValues - delay \n",
		lnp_par.loc_h);fflush(stdout);exit(0);};	

#ifdef LIFCAneuron
      /*PSP 2017-plasticity if (!((synList[i].preSynNeuralKind == excitatoryLbExc) ||
	    (synList[i].preSynNeuralKind == excitatoryLaExc) ||
            (synList[i].preSynNeuralKind == inhibitoryLaInh)))
	{printf("ERROR on loc_h=%d -chkSynListInitVal - preSynNeuralKind \n",
	lnp_par.loc_h); fflush(stdout); exit(0);};*/
#endif
   };
};

#if defined(makeActiveLTD) || defined (makeActiveLTP)
void localNetClass::checkBackwardSynCountAfterInit() {
    uint32_t i, checkBackwardSynCount;
    checkBackwardSynCount=0;

    for(i=0;i<lnp_par.locN;i++)
      {
	checkBackwardSynCount+=n[i].getSizeBackwardSynList();
	DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
	  printf("loc_h=%d, neuron=%d, N_pre=%d\n", 
		 lnp_par.loc_h, i, n[i].getSizeBackwardSynList());
	  fflush(stdout);
	  DPSNNverboseEnd();
      }
    if (checkBackwardSynCount!=backwardSynCount) {
      printf("ERROR on loc_h=%d - check backward syn count failed check=%d back=%d\n",
	     lnp_par.loc_h,checkBackwardSynCount,backwardSynCount); fflush(stdout);exit(0);}
  }
#endif

void localNetClass::reportLocalNetAfterInit() {
  DPSNNverboseStart(false,1,0);
  if(lnp_par.loc_h < 8) {
    printf(
  "490- h=%03d: %d forw syn (%d to itself, %d to remote) %d backw (%d from rem)\n",
    lnp_par.loc_h, 
    forwardSynCount, 
    synTargetHostDistribution.synCount[lnp_par.loc_h],
    (forwardSynCount - synTargetHostDistribution.synCount[lnp_par.loc_h]), 
    backwardSynCount, 
    (backwardSynCount - synSourceHostDistribution.synCount[lnp_par.loc_h]));
  };
  DPSNNverboseEnd();
};

void localNetClass::initChronoWindow() {
      uint32_t startChrono, chronoWindow;
      char *ps_startChrono, *ps_chronoWindow;
      
      ps_startChrono = getenv("env_startPartialChrono_ms");
      ps_chronoWindow = getenv("env_partialChronoWindow_ms");

      if(ps_startChrono!=NULL)
	startChrono = atoi(ps_startChrono);
      else
	startChrono = 0;
      if(ps_chronoWindow!=NULL)
	chronoWindow = atoi(ps_chronoWindow);
      else
	chronoWindow = 0;

      startPartialChrono_ms = startChrono;
      stopPartialChrono_ms = startChrono + chronoWindow;
}

void localNetClass::clearFiringsInChronoWindow() {
  N_firingsInChronoWindow = 0;
  minFiringsInChronoWindow = -1;
  maxFiringsInChronoWindow = 0;
  meanFiringsInChronoWindow = 0;
  varianceFiringsInChronoWindow = 0;
  sigmaFiringsInChronoWindow = 0;
  coeffOfVariationFiringsInChronoWindow = 0;
  M2FiringsInChronoWindow = 0;
  numLapsInChronoWindow = 0;
}

void localNetClass::prepareAxonalSpikeBuffers() {

  uint32_t synTargetCount = 0U;
  uint32_t synSourceCount = 0U;

  //reset the circular buffer for incoming axonal spikes
  for(uint32_t h = 0; h < lnp_par.globH; h++)
    axonalSpikeScheduler[h].reset();

  // REMEMBER: max number of spikes to another process is when all
  // neurons spike, so it must equal how many neurons the spiking
  // process owns (=lnp_par.locN); if the process has no connection
  // with the h-th process's neurons (synCount[h] is null), then no
  // allocation is needed

  // total outgoing synapses per process and total incoming synapses
  // per process (with delay queues)
  for(uint32_t h = 0; h < lnp_par.globH; h++) {
    synTargetCount += ((synTargetHostDistribution.synCount[h] != 0) ? lnp_par.locN : 0);
    synSourceCount += ((synSourceHostDistribution.synCount[h] != 0) ? lnp_par.locN : 0);
  }

  // Reserve contiguous memory pools for outgoing, incoming and
  // delayed incoming spikes arrays
  axonalSpikeDataOnlyClass *pSynTargetPool        = (axonalSpikeDataOnlyClass *)malloc(synTargetCount*sizeof(axonalSpikeDataOnlyClass));
  axonalSpikeDataOnlyClass *pSynSourcePool        = (axonalSpikeDataOnlyClass *)malloc(synSourceCount*sizeof(axonalSpikeDataOnlyClass));
  axonalSpikeDataOnlyClass *pDelayedSynSourcePool = (axonalSpikeDataOnlyClass *)malloc(synSourceCount*sizeof(axonalSpikeDataOnlyClass));

  // assign memory from the pool to...
  for(uint32_t h = 0; h < lnp_par.globH; h++) {

    // ... outgoing spikes array...
    forwardAxonalSpikes[h].list = pSynTargetPool;
    pSynTargetPool += ((synTargetHostDistribution.synCount[h] != 0) ? lnp_par.locN : 0);

    // ... incoming spikes array...
    backwardAxonalSpikes[h].list = pSynSourcePool;
    pSynSourcePool += ((synSourceHostDistribution.synCount[h] != 0) ? lnp_par.locN : 0);

    // ... and delayed spikes array
    delayedBackwardAxonalSpikes[h].list = pDelayedSynSourcePool;
    pDelayedSynSourcePool += ((synSourceHostDistribution.synCount[h] != 0) ? lnp_par.locN : 0);
  }

  // cycle over delays and then processes (so 'h' is the fastest index
  // instead of 'd' and all process synapses are contiguous in memory)
  uint64_t totSchedulerMem = 0U;

  for(uint32_t d = 0; d < DSD__maxD; d++) {

    // assign memory from the pool to the queue of delayed spikes array
    axonalSpikeDataOnlyClass *pCircBufSynSourcePool = (axonalSpikeDataOnlyClass *)malloc(synSourceCount*sizeof(axonalSpikeDataOnlyClass));
    for(uint32_t h = 0; h < lnp_par.globH; h++) {
      axonalSpikeScheduler[h].circBuffer[d].list = pCircBufSynSourcePool;
      pCircBufSynSourcePool +=	((synSourceHostDistribution.synCount[h] != 0) ? lnp_par.locN : 0);
    }

    totSchedulerMem += synSourceCount*sizeof(axonalSpikeDataOnlyClass);
  }

  // this is added to comply with previous code...
  for(uint32_t h = 0; h < lnp_par.globH; h++) {
    if(synSourceHostDistribution.synCount[h] != 0) {
      DPSNNverboseStart(false,1,0);
      printf("localNet_init: From source h=%03u to loc_h=%03u allocating axonalSpikeScheduler memory = %lu byte\n",
  	     h, lnp_par.loc_h, DSD__maxD*lnp_par.locN*sizeof(axonalSpikeDataOnlyClass));
      fflush(stdout);
      DPSNNverboseEnd();
    }
  }

  // fill the displacements arrays to be used in the MPI_Alltoallv
  // spike exchange between processes; their adddresses are pulled
  // from the messagePassingClass they belong to with ad-hoc added
  // methods...
  if (lnp_par.pktLength != 0) {
    pMessagePassing->init_OffsetBaseAddr(forwardAxonalSpikes[0].list, backwardAxonalSpikes[0].list);

    // ... and the offsets are computed with the pointer-to-integer
    // conversion functions of MPI, MPI_Get_address()
    for(uint32_t h = 0; h < lnp_par.globH; h++) {
      pMessagePassing->set_fwdOffsets(forwardAxonalSpikes[h].list, h);
      pMessagePassing->set_bwdOffsets(backwardAxonalSpikes[h].list, h);
    };
  } // pktLength

  DPSNNverboseStart(false,1,0);  
  printf("localNet_init: On loc_h=%03d allocated axonalSpikeScheduler memory = %lu Kb\n",
	 lnp_par.loc_h, totSchedulerMem/1024); 
  fflush(stdout);
  DPSNNverboseEnd();

}
