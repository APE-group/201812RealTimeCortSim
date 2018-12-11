// DPSNN_messagePassing.c
// DPSNN_*.* Distribution/Parallelization 
// performed by Pier Stanislao Paolucci (Roma, Italy project start date 2011),
// starting from the sequral network with axonal conduction delays and STDP
// Created by Eugene M. Izhikevich, May 17, 2004, San Die, CA

// reference paper, named [IzhPol] in the following 
// Eugene M. Izhikevich "Polychronization: Computation with Spikes" 
// Neural Computation 18, 245-282 (2006)
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "DPSNN_debug.h"
#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_messageMeasure.h"
#include "DPSNN_synapse.h"
/* #include "DPSNN_save_and_tell.h" */

#ifdef MPIandDALenvironmentSelected
  #include "DPSNN_messagePassing.h"
  #include "localNetProcess.h"
#else //DALonlyEnvironmentSelected
  #include "localNetProcess.h"
#endif

#ifdef APENET_RDMA
  #include "DPSNN_apenet.h"
#endif

/* here is the definition of the unrolled A2A */
/* BE VERY CAREFUL! */
/* this function doesn't compute the actual size of the derived data type */
/* it is called with, it just uses the MPI_spike size computed elsewhere! */
extern int MPIglobH;
extern int MPIprocessRank;
MPI_Aint sizeof_MPI_spike = 0;

/* Unrolled A2A */
// Perform pairwise exchange - starting from 1 so the local copy
// is last it needs rank (lnp_par.loc_h) and size
// (lnp_par.globH) of communicator which are global
static int unrolled_A2A(const void *sbuf, int scount, MPI_Datatype sdtype,
			      void *rbuf, int rcount, MPI_Datatype rdtype,
			MPI_Comm comm) {
  int sendto, recvfrom;
  void *tmpsend, *tmprecv;
  MPI_Request req;

  for (int32_t step = 1; step < MPIglobH + 1; step++) {

    /* Determine sender and receiver for this step. */
    sendto   =            (MPIprocessRank + step) % MPIglobH;
    recvfrom = (MPIglobH + MPIprocessRank - step) % MPIglobH;

    /* Determine sending and receiving locations */
    tmpsend = (char*)sbuf + sizeof_MPI_spike * (ptrdiff_t)scount * (ptrdiff_t)sendto;
    tmprecv = (char*)rbuf + sizeof_MPI_spike * (ptrdiff_t)rcount * (ptrdiff_t)recvfrom;

    /* post new irecv */
    MPI_Irecv(tmprecv, rcount, rdtype, recvfrom, MPI_ANY_TAG, comm, &req);

    /* send data to children */
    MPI_Send(tmpsend, scount, sdtype, sendto, 0, comm);

    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }

  return MPI_SUCCESS;
}

/* typedef for pointer to function */
/* to overload the collective functions from standard MPI or unrolled Send/Recv loop */
/* actual choice is performed in the init class according to launch parameters */
typedef int (*f_A2A)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, MPI_Comm);

/* put here, the pointer to function is global*/
static f_A2A my_A2A = NULL;

void messagePassingClass::init(const struct DPSNN_parameters lnp_par_initValue, DALProcess *p_initValue, 
			       const int thisTime_ms) {

  lnp_par=lnp_par_initValue;
  p=p_initValue;

  /* here starts a part needed for fixed packet length: */
  if (lnp_par.pktLength != 0) {
    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_spike);
    MPI_Type_commit(&MPI_spike);

    //this is the correct way to know the size of a composite type
    MPI_Type_get_extent(MPI_spike, &lb_MPI_spike, &sizeof_MPI_spike);
  }
  /* here ends a part needed for fixed packet length */

  switch (lnp_par.collectivesAlgo) {

  case standardMPI:
    my_A2A = MPI_Alltoall;
    break;

  case unrolledMPI:
    my_A2A = unrolled_A2A;
    break;

  default:
    {printf(
	    "ERROR unrecognized collectivesAlgo set by the environment\n");
      fflush(stdout);exit(0);
    };
  }
}

/* here starts a part needed for fixed packet length: */
void messagePassingClass::init_OffsetBaseAddr(axonalSpikeDataOnlyClass *fwd, axonalSpikeDataOnlyClass *bwd) {
  MPI_Get_address(fwd,&MPI_fwdOffsetBaseAddress);
  MPI_Get_address(bwd,&MPI_bwdOffsetBaseAddress);
}

void messagePassingClass::set_fwdOffsets(axonalSpikeDataOnlyClass *spikeAddress, int h) {
  MPI_Aint MPI_spikeAddress = MPI_fwdOffsetBaseAddress;

  if (spikeAddress)
    MPI_Get_address(spikeAddress, &MPI_spikeAddress);

  axSpike_forwardOffset[h] = (MPI_spikeAddress - MPI_fwdOffsetBaseAddress)/sizeof_MPI_spike;
}

void messagePassingClass::set_bwdOffsets(axonalSpikeDataOnlyClass *spikeAddress, int h) {
  MPI_Aint MPI_spikeAddress = MPI_bwdOffsetBaseAddress;

  if (spikeAddress)
    MPI_Get_address(spikeAddress, &MPI_spikeAddress);

  axSpike_backwardOffset[h] = (MPI_spikeAddress - MPI_bwdOffsetBaseAddress)/sizeof_MPI_spike;
}
/* here ends a part needed for fixed packet length */

void messagePassingClass::barrier() {
  #ifdef MPIandDALenvironmentSelected
	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Barrier(comm);
  #else
    #warning "MPI barrier not implemented under DAL"
  #endif
}

void messagePassingClass::sendForwardSynListDimToRemoteHosts(
	synapticDistributionClass *synTargetHostDistribution, 
	synapticDistributionClass *synSourceHostDistribution,
	const uint32_t thisTime_ms) 
{
  /*      pSynTargetHostDistribution = synTargetHostDistribution;
      pSynSourceHostDistribution = synSourceHostDistribution;
  */
      // if the probability of reaching different hosts is homogenous, 
      // we expect a mean of M/globH synapses toward each target 
      // Under this simplistic assumption, valid for smal scale neural nets
      // we are going to perform two steps:
      // 1- an ALLTOALL step to inform each target 
      // about how many synaptical connections will be requested 
      // by the following step 
      // ALLTOALLV(...) which actually transfers the list
      // of requested synapses to be established 
      // note: synapses toward the same host could be managed locally, 
      // so a fake request of 0 connections COULD be issued by each host toward itself
    #ifdef MPIandDALenvironmentSelected
      {
	/* MPI_Comm comm; */
	/* comm=MPI_COMM_WORLD; */

  #ifdef APENET_RDMA
      {
	ape_descriptor = (APE_RDMA_desc_t*)malloc(sizeof(APE_RDMA_desc_t));
	uint32_t s,h;

	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	printf("PID %d on %s ready loc_h=%d\n", getpid(), hostname, lnp_par.loc_h);
	DPSNNverboseEnd();
	
	int ret = APE_RDMA_Init(ape_descriptor, lnp_par.loc_h, lnp_par.globH );
	if(ret)
	  {
	    DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms);	 	   
	    printf("391- ERROR while init apenet device h=%03d\n", lnp_par.loc_h);
	    DPSNNverboseEnd();
	    exit(0);
	  }

	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	printf("ape malloc  ms=%d on h=%d\n", thisTime_ms, lnp_par.loc_h);
	DPSNNverboseEnd();

	//usleep(1000000*lnp_par.loc_h);

	for(h = 0; h < lnp_par.globH ;h++)
	  {

	    apeForwardBuffer[h] = (char*)memalign(sysconf(_SC_PAGESIZE), APE_BUFFER_SIZE);
	    apeBackwardBuffer[h]  = (volatile char*)memalign(sysconf(_SC_PAGESIZE), APE_BUFFER_SIZE);
	    /* apeForwardBuffer[h] = (char*)malloc(APE_BUFFER_SIZE); */
	    /* apeBackwardBuffer[h]  = (char*)malloc(APE_BUFFER_SIZE); */

	    if( (apeForwardBuffer[h] == NULL) || (apeBackwardBuffer[h] == NULL) )
	      {
		DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms);
		printf("ape malloc error  ms=%d on h=%d\n", thisTime_ms, lnp_par.loc_h);
		DPSNNverboseEnd();
	      }
	    am_vaddr_t vaddr;
	    am_context_t p_ctx = (am_context_t)0xc1a0c1a0c1a0c1a0ULL;

	    if((ape_descriptor->coords)[lnp_par.loc_h].u == (ape_descriptor->coords)[h].u)
	      {
	    	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	    	printf("ape skip register buf %p for h=%d ms=%d on h=%d\n", vaddr, h, thisTime_ms, lnp_par.loc_h);
	    	DPSNNverboseEnd();

	      } else {
		int retcode = am_register_buf(ape_descriptor->link, &vaddr ,
					      (void*)(apeBackwardBuffer[h]), 
					      APE_BUFFER_SIZE, 
					      p_ctx, 
					      AM_PERSISTENT_BUF);
	
		if(retcode){
		  DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms);
		  printf("error in ape register buf %p size=%d ms=%d on h=%d for h=%d\n", vaddr, APE_BUFFER_SIZE, thisTime_ms, lnp_par.loc_h, h);
		  DPSNNverboseEnd();
		} else {
		  DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms);
		  printf("ape register buf %p size=%d ms=%d on h=%d\n", vaddr, APE_BUFFER_SIZE, thisTime_ms, lnp_par.loc_h);
		  DPSNNverboseEnd();
		}
	      }
	  }
	DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms);
	printf("ape register DONE, let us exchange addresses ms=%d on h=%d\n",  thisTime_ms, lnp_par.loc_h);
	DPSNNverboseEnd();

	//exchange apeBackwardBuffer addresses
	MPI_Alltoall((void *)apeBackwardBuffer, 1, MPI_UNSIGNED_LONG_LONG, 
	    (void *)apeRemoteBackwardBuffer, 1, MPI_UNSIGNED_LONG_LONG,
	    MPI_COMM_WORLD);


      }
  #endif //APENET_RDMA

	if(lnp_par.globH==1) {//single process case
	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	    printf("260- single node sendForwSynListDimToRemHosts - no MPI!\n");
	  DPSNNverboseEnd();
	  synSourceHostDistribution->synCount[0]=
	    synTargetHostDistribution->synCount[0];
	  synSourceHostDistribution->synOffset[0]=
	    synTargetHostDistribution->synOffset[0];
	}else{//more than a single process running
	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	    printf(
	    "270- sendForwSynListDimToRemHosts - loc_h=%d before MPI_Alltoall\n",
		   lnp_par.loc_h);
	  DPSNNverboseEnd();


	  MPI_Alltoall(
	    (void *)&(synTargetHostDistribution->synCount[0]), 1, MPI_INT, 
	    (void *)&(synSourceHostDistribution->synCount[0]), 1, MPI_INT,
	    MPI_COMM_WORLD);
	};
      }
    #else // DALonlyEnvironmentSelected
      { 
	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	printf(
	"265- %d: sendForwSynListDimToRemHosts - start allToAll phase\n",
	lnp_par.loc_h);
	DPSNNverboseEnd();

	CREATEPORTVAR(O);
	CREATEPORTVAR(I);
       	uint32_t h;
	for (h = 0; h < lnp_par.globH;h++) {
	  CREATEPORT(O, PORT_O, 1, h, lnp_par.globH);
	  DAL_write((void*) O, 
		    &(synTargetHostDistribution->synCount[h]), sizeof(int), p);
	};
	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);	 
	printf(
	"275- %d: sendForwSynListDimToRemHosts - launched its DAL_write()s \n",
	       lnp_par.loc_h);
	DPSNNverboseEnd();
	for (h = 0; h < lnp_par.globH;h++) {
	  CREATEPORT(I, PORT_I, 1, h, lnp_par.globH);
	  DAL_read((void*) I, 
		   &(synSourceHostDistribution->synCount[h]), sizeof(int), p);
	};
	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);	 
	printf(
	"285- %d: sendForwSynListDimToRemHosts - completed its DAL_read()s \n",
	       lnp_par.loc_h);
	DPSNNverboseEnd();
      }
    #endif // all to all completed
	{
	  //calculation of the offsets where the incoming synapses 
	  //will be stored in the future
	  uint32_t h, totBackwardSynCount;
	  totBackwardSynCount=synSourceHostDistribution->synCount[0];
	  synSourceHostDistribution->synOffset[0]=0;
	  for(h=1;h<lnp_par.globH;h++){
	    totBackwardSynCount+=synSourceHostDistribution->synCount[h];
	    synSourceHostDistribution->synOffset[h]=
	      synSourceHostDistribution->synOffset[h-1]+
	      synSourceHostDistribution->synCount[h-1];	
	  }
	DPSNNverboseStart(false,1,0);	 
	  printf(
         "285- h=%d, alltoall completed %d incoming syn are expected\n", 
		     lnp_par.loc_h, totBackwardSynCount);
	DPSNNverboseEnd();
	}		

    };

void messagePassingClass::sendForwardSynListToRemoteHosts(
  synapseClass *forwardSynList, 
  synapticDistributionClass *synTargetHostDistribution, 
  synapseClass *backwardSynList, 
  synapticDistributionClass *synSourceHostDistribution,
  const uint32_t thisTime_ms) 
{
  if(lnp_par.globH == 1) {
    //single process case - no MPI nor DAL required
    uint32_t s;
    for(s=0; s < (synTargetHostDistribution->synCount[0]); s++) {
            backwardSynList[s]=forwardSynList[s];
    };
  } else {
    DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);	 
	  printf("360- sendForwardSynListToRemoteHosts - on h=%03d\n",
		lnp_par.loc_h);
    DPSNNverboseEnd();
    #ifdef MPIandDALenvironmentSelected
      {

	/* MPI_Comm comm; */

	int synapse_MPI_datatype_numberOfBlocks;
	int synapse_MPI_datatype_arrayOfBlocklengths[2];
	MPI_Aint synapse_MPI_datatype_arrayOfDisplacements[2];
	MPI_Datatype synapse_MPI_arrayOfDatatype[2];
	MPI_Datatype synapse_MPI_datatype;

	/* comm=MPI_COMM_WORLD; */

  	synapse_MPI_datatype_numberOfBlocks=1;
#if defined(makeActiveLTD) || defined (makeActiveLTP)
	//PSP 2017-09 Plasticity (actually, during transmission it is not needed to transport
	//derivative and activation time. Howover, for easier debugging...)
	synapse_MPI_datatype_arrayOfBlocklengths[0]=5;
#else
	synapse_MPI_datatype_arrayOfBlocklengths[0]=3;
#endif
	synapse_MPI_arrayOfDatatype[0]=MPI_INT;
	synapse_MPI_datatype_arrayOfDisplacements[0]=0;

	MPI_Type_create_struct(
		synapse_MPI_datatype_numberOfBlocks,
		synapse_MPI_datatype_arrayOfBlocklengths,
		synapse_MPI_datatype_arrayOfDisplacements,
		synapse_MPI_arrayOfDatatype,
		&synapse_MPI_datatype);
	MPI_Type_commit(&synapse_MPI_datatype);
     

	//on MPI it seems from experiments on spikes ALLTOALLV that the transmission 
	//from a process to itself in ALLTOALLV is already optimized away
	//from interprocess communication to intraprocess comm.
	//Therefore, no need for further opt here
	MPI_Alltoallv((void *) forwardSynList, 
			(int *) synTargetHostDistribution->synCount, 
			(int *) synTargetHostDistribution->synOffset, 
			synapse_MPI_datatype,
			(void *) backwardSynList, 
			(int *) synSourceHostDistribution->synCount, 
			(int *) synSourceHostDistribution->synOffset, 
			synapse_MPI_datatype,
			MPI_COMM_WORLD);

	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);	 	   
          printf("390- MPI_Alltoallv completed on h=%03d\n", lnp_par.loc_h);
	DPSNNverboseEnd();
      };
    #else // in MPIandDALenvironmentSelected, DAL_used
      {       
	CREATEPORTVAR(O);
	CREATEPORTVAR(I);
        {
	  uint32_t h;
	  int * bufAddress;
	  int bufSize;

          {//intra-process communication
	    uint32_t s;
            for(s=0; 
		s < (synTargetHostDistribution->synCount[lnp_par.loc_h]); 
		s++) {
     backwardSynList[s + synSourceHostDistribution->synOffset[lnp_par.loc_h]] =
     forwardSynList[s + synTargetHostDistribution->synOffset[lnp_par.loc_h]];
	    }; 
          };
	  //interprocess comm. i.e. h != lnp_par.loc_h
	  for (h = 0; h < lnp_par.loc_h; h++) {
	    CREATEPORT(O, PORT_O, 1, h, lnp_par.globH);
	    bufAddress = 
              (int *) &forwardSynList[synTargetHostDistribution->synOffset[h]];
	    bufSize = (3 * sizeof(int)) * 
                      (synTargetHostDistribution->synCount[h]);

	    DAL_write((void*) O, bufAddress, bufSize, p);
	  };

	  for (h = lnp_par.loc_h + 1 ; h < lnp_par.globH; h++) {
	    CREATEPORT(O, PORT_O, 1, h, lnp_par.globH);
	    bufAddress = 
              (int *) &forwardSynList[synTargetHostDistribution->synOffset[h]];
	    bufSize = (3 * sizeof(int)) * 
                      (synTargetHostDistribution->synCount[h]);

	    DAL_write((void*) O, bufAddress, bufSize, p);
	  };
     
	  for (h = 0; h < lnp_par.loc_h ;h++) {
            CREATEPORT(I, PORT_I, 1, h, lnp_par.globH);
	    bufAddress = 
              (int *) &backwardSynList[synSourceHostDistribution->synOffset[h]];
	    bufSize = (3 * sizeof(int) ) * 
	              (synSourceHostDistribution->synCount[h]);
	    DAL_read((void*) I, bufAddress, bufSize, p);
	  };

 	  for (h = lnp_par.loc_h + 1; h < lnp_par.globH;h++) {
            CREATEPORT(I, PORT_I, 1, h, lnp_par.globH);
	    bufAddress = 
              (int *) &backwardSynList[synSourceHostDistribution->synOffset[h]];
	    bufSize = (3 * sizeof(int) ) * 
	              (synSourceHostDistribution->synCount[h]);
	    DAL_read((void*) I, bufAddress, bufSize, p);
	  };
        };
      };
    #endif // MPIandDALenvironmentSelected
  };
};

uint32_t singleHostAxonalSpikeExchEmulCount;

void messagePassingClass::exchangeAxonalSpikesDim(
    const forwardAxonalSpikesClass  *pForwardAxonalSpikes,
    const synapticDistributionClass *pSynTargetHostDistribution, 
    backwardAxonalSpikesClass *pBackwardAxonalSpikes,
    const synapticDistributionClass *pSynSourceHostDistribution,
    const uint32_t thisTime_ms)
{

  uint32_t h;
  uint32_t sendDimCounter;

  for(h=0;h<lnp_par.globH;h++) {
    pBackwardAxonalSpikes[h].count = 0;
    pBackwardAxonalSpikes[h].expectedCount = 0;
  };
 
  for(h=0;h<lnp_par.globH;h++) {
    if(pBackwardAxonalSpikes[h].count!=0) {
      printf("ERROR back count != 0 when starting exchangeAxonalSpikes\n");
      fflush(stdout);exit(0);};
    if(pBackwardAxonalSpikes[h].expectedCount!=0) {
      printf("ERROR back expected count != 0 when starting exchangeAxonalSpikes\n");
      fflush(stdout);exit(0);};
    if(pForwardAxonalSpikes[h].count > DSD__maxAxonalSpike) {
      printf("ERROR:messPass %d ms:pForwAxSpik[h=%d].count=%d>DSD__maxAxonalSpike\n",
	     thisTime_ms,lnp_par.loc_h, pForwardAxonalSpikes[h].count);
      fflush(stdout);exit(0);};
  }

  if(lnp_par.globH==1) {
    h=0;
    if(pForwardAxonalSpikes[h].count > 0) {
 
      pBackwardAxonalSpikes[h].expectedCount = 
        pForwardAxonalSpikes[h].count;

      DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
        printf(
        "exchangeAxonalSpikesDim: forward spikes count=%d, backward exp=%d\n",
	     pForwardAxonalSpikes[h].count,
	     pBackwardAxonalSpikes[h].expectedCount);
      DPSNNverboseEnd();
    };

  } else { // lnp_par.globH > 1

    DPSNNverboseStart(false,thisTime_ms,4);
      for (h = 0; h < lnp_par.globH;h++) {
	
        printf("802- %d ms h=%03d->targ %03d: %d spikes to transm - excAxSpDim\n", 
	     thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count); 
	fflush(stdout);
      };
    DPSNNverboseEnd();

    /* this is where you excise if fixed packet length is chosen */
    if (lnp_par.pktLength == 0) {
      // here it prepares two arrays of flags
      // one sending and one for receiving
      // .Count[h] =0/1 - traffic between me and 'h': no/yes
      // .Offset[h]=sum (of sizes) of all non-0 elements in .Count[h]
      {
	uint32_t sendDimCount = 0, recDimCount = 0;
	uint32_t sendPayloadCount = 0, recPayloadCount = 0;

	axSpike_forwardOffset[0]=0;
	axSpike_backwardOffset[0]=0;
	for(h=0;h<lnp_par.globH;h++) {
	  if(h > 0)
	    axSpike_forwardOffset[h] = 
	      axSpike_forwardOffset[h-1]+
	      axSpike_forwardCount[h-1];
	  if(pSynTargetHostDistribution->synCount[h] != 0 ) {
	    axSpike_forwardCount[h] = 1; // 'h' will receive from me
	    sendDimCount++;
	  }else{
	    axSpike_forwardCount[h] = 0;
	  };
	  if(h>0)
	    axSpike_backwardOffset[h] = 
	      axSpike_backwardOffset[h-1]+
	      axSpike_backwardCount[h-1];
	  if(pSynSourceHostDistribution->synCount[h] != 0 ) {
	    axSpike_backwardCount[h] = 1; // 'h' will send to me
	    recDimCount++;
	  }else{
	    axSpike_backwardCount[h] = 0;
	  };
	}

	// fill up of the sending buffers?
	for(h=0;h<lnp_par.globH;h++) {
	  if(axSpike_forwardCount[h]!=0) {
	    axSpike_forwardPrep[axSpike_forwardOffset[h]] =
	      pForwardAxonalSpikes[h].count;
	  }
	}

	// clean up of the receiving buffers?
	for(int i=0; i < axSpike_backwardOffset[lnp_par.globH-1] + axSpike_backwardCount[lnp_par.globH-1]; i++)
	  axSpike_backwardPrep[i]=0;

	DPSNNverboseStart(true,1,0);
	sendPayloadCount = axSpike_forwardOffset[lnp_par.globH-1]+
	  axSpike_forwardCount[lnp_par.globH-1];
	recPayloadCount = axSpike_backwardOffset[lnp_par.globH-1]+
	  axSpike_backwardCount[lnp_par.globH-1];
	if((sendPayloadCount != sendDimCount) || 
	   (recPayloadCount != recDimCount)) 
	  {
	    printf(
		   "ERROR spikeDim-A t=%d,h=%d, sendPayloadCount=%d, recPayloadCount=%d \n",
		   thisTime_ms, lnp_par.loc_h, sendPayloadCount,recPayloadCount);
	    printf(
		   "ERROR spikeDim-B t=%d,h=%d, sendDimCount=%d, recDimCount=%d \n",
		   thisTime_ms,lnp_par.loc_h,sendDimCount,recDimCount);
	    fflush(stdout);exit(0);
	  }
	DPSNNverboseEnd();
      }
#ifdef MPIandDALenvironmentSelected

      DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
      printf("Dim exchange loc_h=%d ms=%d\n", lnp_par.loc_h, thisTime_ms);
      DPSNNverboseEnd();  

#define axSpikesDimUsesAllToAllv
#ifdef axSpikesDimUsesAllToAllv

      // when here, only the sizes array is ready
      MPI_Alltoallv((void *) axSpike_forwardPrep, 
		    (int *)  axSpike_forwardCount, 
		    (int *)  axSpike_forwardOffset,
		    MPI_INT,
		    (void *) axSpike_backwardPrep,
		    (int *)  axSpike_backwardCount,
		    (int *)  axSpike_backwardOffset,
		    MPI_INT,
		    MPI_COMM_WORLD);  

      DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
      printf("Dim exchange DONE loc_h=%d ms=%d\n", lnp_par.loc_h, thisTime_ms);

      /* char filename[20]; */
      /* sprintf(filename, "%u.stdout", lnp_par.loc_h); */

      /* save_and_tell(filename, */
      /* 		  (void*)axSpike_forwardPrep, (void*)axSpike_backwardPrep, */
      /* 		  axSpike_forwardCount, axSpike_backwardCount, */
      /* 		  axSpike_forwardOffset, axSpike_backwardOffset, */
      /* 		  lnp_par.loc_h, */
      /* 		  lnp_par.globH, */
      /* 		  sizeof(typeof(axSpike_forwardPrep[0]))); */

      DPSNNverboseEnd();  

#else  //axSpikesDimUsesAllToAllv
      {
	MPI_Status status;
	int returnCount;
	for(h=0;h<lnp_par.globH;h++) {
	  if(axSpike_forwardCount[h] != 0) {
	    MPI_Send((void *)&(axSpike_forwardPrep[axSpike_forwardOffset[h]]), 1, MPI_INT,
		     h, thisTime_ms, MPI_COMM_WORLD);
	  }
	}
	for(h=0;h<lnp_par.globH;h++) {
	  if(axSpike_backwardCount[h] != 0) {
	    MPI_Recv((void *)&(axSpike_backwardPrep[axSpike_backwardOffset[h]]), 1, MPI_INT,
		     h, thisTime_ms, MPI_COMM_WORLD, &status);
	    DPSNNverboseStart(false,1,0);
	    MPI_Get_count(&status, MPI_INT, &returnCount);
	    if(returnCount!=1){printf(
				      "ERROR: exAxSpDim: MPI_Get_count != 1\n");
	      fflush(stdout); exit(0);
	    }
	    DPSNNverboseEnd();
	  }
	}
      }
#endif //axSpikesDimUsesAllToAllv
#else // in MPIandDALenvironmentSelected, DAL
      {
	int * bufAddress;
	int bufSize;  
	CREATEPORTVAR(O);
	CREATEPORTVAR(I);
	for (h = 0; h < lnp_par.globH;h++) {
	  if(axSpike_forwardCount[h] != 0) {
	    CREATEPORT(O, PORT_O, 1, h, lnp_par.globH);
	    bufAddress = 
	      (int *)&(
		       axSpike_forwardPrep[axSpike_forwardOffset[h]]);
	    bufSize =  sizeof(int);
	    DAL_write((void*) O, 
		      bufAddress, bufSize, p);
	  };
	};
	for (h = 0; h < lnp_par.globH;h++) {
	  if(axSpike_backwardCount[h] != 0) {
	    CREATEPORT(I, PORT_I, 1, h, lnp_par.globH);
	    bufAddress = 
	      (int *)&(
		       axSpike_backwardPrep[axSpike_backwardOffset[h]]);
	    bufSize = sizeof(int);
	    DAL_read((void*) I,
		     bufAddress, bufSize, p);
	  };
	};
      };
#endif // in MPIandDALenvironmentSelected, end of DAL

      DPSNNverboseStart(false,thisTime_ms,4);
      for (h = 0; h < lnp_par.globH;h++) {
	printf("804- excAxSpDim on h=%d received axSpike: %d ms h=%03d ->targ %03d: %d spike dim to be receiv %d spikes \n",
	       lnp_par.loc_h, thisTime_ms, h, lnp_par.loc_h, 
	       axSpike_backwardCount[h],axSpike_backwardPrep[h]); 
      };
      DPSNNverboseEnd();

      // here the comm buffer is copied back to the working one
      // (not the data array, just the sizes one!)
      for(h=0;h<lnp_par.globH;h++) {
	if(axSpike_backwardCount[h]==0) {
	  pBackwardAxonalSpikes[h].expectedCount=0;
	}else{
	  pBackwardAxonalSpikes[h].expectedCount=
	    axSpike_backwardPrep[axSpike_backwardOffset[h]];
	};
      }
      DPSNNverboseStart(false,thisTime_ms,4);
      for (h = 0; h < lnp_par.globH;h++) {  
	printf("804- excAxSpDim: %d ms h=%03d ->targ %03d: %d bckCount %d bckOff \n", 
	       thisTime_ms, h, lnp_par.loc_h,
	       axSpike_backwardCount[h],axSpike_backwardOffset[h]); 
      };
      DPSNNverboseEnd();

    } else { /* this is where the replacement starts if fixed packet length is chosen */

      // the fake_spike is signed with the max float value
      axonalSpikeDataOnlyClass fake_spike = {FLT_MAX, 0};
      for (uint32_t h=0; h<lnp_par.globH; h++) {

	int lastRowInPacket = (h+1)*lnp_par.pktLength-1;

	// is 'pForwardAxonalSpikes[.].count' the correct count of how
	// many spikes I send to 'h'?
	fake_spike.spikes_in_packet = pForwardAxonalSpikes[h].count;
	axSpike_forwardBuffer[lastRowInPacket] = fake_spike;

	// spill 'pktLength-1' spikes into the packet, starting from the
	// tail and going backwards then adds a 'fake' last one with the
	// total number in the .spikes_in_packet field of the
	// 'pktLength'-th spike
	int count = (pForwardAxonalSpikes[h].count < lnp_par.pktLength
		     ? pForwardAxonalSpikes[h].count
		     : lnp_par.pktLength-1);

	// BEWARE: in non-full packets, we are leaving the unused
	// memory dirty

	// start of tail to pull spikes from (can be 0)
	int startTail = pForwardAxonalSpikes[h].count - count;

	// copy 'pktLength-1' spikes from the tail into
	// axSpike_forwardBuffer; they are contiguous both in the use
	// and the takeoff buffer, so memcpy should be equivalent but
	// a little faster than a 'for' loop
	memcpy(&axSpike_forwardBuffer[h*lnp_par.pktLength],
	       &pForwardAxonalSpikes[h].list[startTail],
	       count*sizeof(axonalSpikeDataOnlyClass));

	// we spilled some spikes, remember how many are left
	axSpike_forwardCount[h] = startTail;

	// if I want padding, here is the place to do it...

      }; // I've scanned each 'h' preparing fixed-length packets

      my_A2A((void *) axSpike_forwardBuffer, lnp_par.pktLength, MPI_spike,
	     (void *) axSpike_backwardBuffer, lnp_par.pktLength, MPI_spike,
	     MPI_COMM_WORLD);

      // here I must decode each received packet's last_spike's to prepare
      // the following (eventually empty) MPI_Alltoallv

      for (uint32_t h=0; h<lnp_par.globH; h++) {
	int lastRowInPacket = (h+1)*lnp_par.pktLength-1;

	// pull the real number of spikes from the last position in
	// the packet and save it in its final destination
	pBackwardAxonalSpikes[h].count =
	  axSpike_backwardBuffer[lastRowInPacket].spikes_in_packet;

	// copy back as many spikes there are in the packet
	int count = (pBackwardAxonalSpikes[h].count < lnp_par.pktLength
		     ? pBackwardAxonalSpikes[h].count
		     : lnp_par.pktLength-1);

	// start of tail to attach spikes to (can be 0)
	int startTail = pBackwardAxonalSpikes[h].count - count;

	// copy 'pktLength-1' spikes from axSpike_backwardBuffer into
	// the tail; spikes are contiguous both in the landing buffer
	// and in their use buffer, so memcpy should be equivalent but
	// a little faster than a 'for' loop
	memcpy(&pBackwardAxonalSpikes[h].list[startTail],
	       &axSpike_backwardBuffer[h*lnp_par.pktLength],
	       count*sizeof(axonalSpikeDataOnlyClass));

	// we filled the tail with spikes, how many yet to fill
	axSpike_backwardCount[h] = startTail;

	// 'pBackwardAxonalSpikes[.].count' is set to the real size while
	// '.expectedCount' is unused, they are checked for equality later
	// on, so the former is dumped into the latter
	pBackwardAxonalSpikes[h].expectedCount = pBackwardAxonalSpikes[h].count;

      }; // I've scanned each 'h'
    } /* this is where the replacement ends if fixed packet length is chosen */

  } //end if globH=1 else if globH!=1

  DPSNNverboseStart(false,thisTime_ms,4);
    for (h = 0; h < lnp_par.globH;h++) {  
      printf(
	"805- excAxSpDim FINAL: %d ms h=%03d ->targ %03d: %d Spikes to be rec \n", 
	 thisTime_ms, h, lnp_par.loc_h, pBackwardAxonalSpikes[h].expectedCount
	); 
    };
  DPSNNverboseEnd();

  //check receiving buffer size before transmitting
  for (h = 0; h < lnp_par.globH;h++) {
    if(pBackwardAxonalSpikes[h].expectedCount > DSD__maxAxonalSpike) {
       printf(
"ERROR:messPass %d ms:pBackAxSpik[h=%d].expectedCount=%d>DSD__maxAxonalSpike\n",
	thisTime_ms,lnp_par.loc_h, pBackwardAxonalSpikes[h].expectedCount);
        fflush(stdout);exit(0);
    };
  };

  sendDimCounter = 1;
  for(h=0; h<lnp_par.globH; h++) {
    if(h != lnp_par.loc_h)
      if(pSynTargetHostDistribution->synCount[h] != 0 )
	spikeDimSize.accumulateMsgSize(sendDimCounter);
  }
};

void messagePassingClass::exchangeAxonalSpikes(
    forwardAxonalSpikesClass *pForwardAxonalSpikes,
    backwardAxonalSpikesClass *pBackwardAxonalSpikes,
    const uint32_t thisTime_ms) 
{
  uint32_t h,i;

  uint32_t a,count,offset; // these var's are used only when packet length is not fixed

  // remove this zeroing when fixed packet length is
  // chosen; offsets arrays must not be erased since
  // they are initialized elsewhere!
  if (lnp_par.pktLength == 0) {

    // clean up all buffers
    for(h=0;h<lnp_par.loc_h;h++) {
      axSpike_forwardCount[h] = 0;
      axSpike_forwardOffset[h] = 0;
      axSpike_backwardCount[h] = 0;
      axSpike_backwardOffset[h] = 0;
      pBackwardAxonalSpikes[h].count=0;
    };
  }

  //debug print range: 852- 898-

  DPSNNverboseStart(false,thisTime_ms,0);
        printf("852- exch Ax Spik:= START on h=%d\n", 
	     lnp_par.loc_h); 
  DPSNNverboseEnd(); 

  DPSNNverboseStart(false,thisTime_ms,0);
    for (h = 0; h < lnp_par.globH;h++) {   
        printf(
       "853- exch Ax Spik:= %d ms sh%03d ->th%03d: %d spikes to be sent\n", 
       thisTime_ms, lnp_par.loc_h, h, 
             pForwardAxonalSpikes[h].count);
        printf(
       "854- exch Ax Spik:= %d ms sh%03d ->th%03d: %d spikes expected\n", 
       thisTime_ms, h, lnp_par.loc_h, 
             pBackwardAxonalSpikes[h].expectedCount); 
    };  
  DPSNNverboseEnd();  

  //check receiving buffer size before transmitting
  for (h = 0; h < lnp_par.globH;h++) {   
    if(pBackwardAxonalSpikes[h].expectedCount > DSD__maxAxonalSpike) {
      printf(
      "ERROR:messPass %u ms:pBackAxSpik[h=%d].expectedCount=%d>DSD__maxAxonalSpike\n",
	thisTime_ms,lnp_par.loc_h, pBackwardAxonalSpikes[h].expectedCount);
        fflush(stdout);exit(0);
    };
  };


#ifndef APENET_RDMA

  // substituting single process comm with self-copy
  if(lnp_par.globH == 1) {
    if(lnp_par.loc_h!=0) {
        printf("ERROR: loc_h !=0 on a system with globH = 1 \n");
        fflush(stdout);exit(0);
    }

    h=0;

    {
      if(pBackwardAxonalSpikes[h].expectedCount > 0) {
	for(i=0; i < pBackwardAxonalSpikes[h].expectedCount; i++) {
	  pBackwardAxonalSpikes[h].list[i] = pForwardAxonalSpikes[h].list[i];
	  pBackwardAxonalSpikes[h].count ++;
  	  if(!((pForwardAxonalSpikes[h].list[i].originalEmissionTime >= thisTime_ms-1) &&
	       (pForwardAxonalSpikes[h].list[i].originalEmissionTime < thisTime_ms))){
	    printf("ERROR at %d ms wrong original emission time %.7f for the neu %d preparing ax spike MPI\n",
		   thisTime_ms,pForwardAxonalSpikes[h].list[i].originalEmissionTime,
		   pForwardAxonalSpikes[h].list[i].pre_glob_n);
	    fflush(stdout);exit(0);
	  };

	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("exch Ax Spik: after increm of the ax spike counter %d\n",
		 pBackwardAxonalSpikes[h].count);
	  DPSNNverboseEnd();  
	};
      }
    }
  
  } else { //not a single h

    if (lnp_par.pktLength == 0) { // replace this part of function if fixed packet length is chosen

      { //prepare forward offset and count 
	offset = 0;

	for(h=0;h<lnp_par.globH;h++) {
	  axSpike_forwardOffset[h] = offset;
	  count = pForwardAxonalSpikes[h].count;
	  axSpike_forwardCount[h] = count;
	  offset += count;
	};
      }

      { //copy data into the transmission buffer
	for(h=0;h<lnp_par.globH;h++) {
	  count = axSpike_forwardCount[h];
	  offset = axSpike_forwardOffset[h];
	  for(a=0; a< count; a++) {
	    axSpike_forwardBuffer[offset + a] = 
	      pForwardAxonalSpikes[h].list[a];
#ifdef LIFCAneuron
	    if(!((axSpike_forwardBuffer[offset + a].originalEmissionTime >= thisTime_ms-1) &&
		 (axSpike_forwardBuffer[offset + a].originalEmissionTime < thisTime_ms))){
	      printf("ERROR from loc_h=%d to target h=%d at %d ms wrong original emission time %.7f preparing ax spike MPI\n",
		     lnp_par.loc_h,h,thisTime_ms, axSpike_forwardBuffer[offset + a].originalEmissionTime);
	      fflush(stdout);exit(0);
	    };
#else // LIFCAneuron
	    if(axSpike_forwardBuffer[offset + a].originalEmissionTime
	       != thisTime_ms) {
	      printf("ERROR wrong original emission time preparing ax spike MPI\n");
	      fflush(stdout);exit(0);
	    };
#endif // LIFCAneuron
	  }
	};
      };

      { //debug section - print forward data 
	DPSNNverboseStart(false,thisTime_ms,0);        
	for(h=0;h<lnp_par.globH;h++) {
	  count = axSpike_forwardCount[h];
	  offset =  axSpike_forwardOffset[h];
	  for(a=0; a< count; a++) {
	    axonalSpikeDataOnlyClass axSp;
	    axSp = axSpike_forwardBuffer[offset + a];
	    printf(
		   "858- exch Ax Spik- on h=%d prep2 ALLTOALLV: %d ms AxSp(glob_n=%d,off=%d,orig %.7f ms)\n",
		   lnp_par.loc_h, thisTime_ms, 
		   axSp.pre_glob_n, offset + a, axSp.originalEmissionTime);
	  }
	};
	DPSNNverboseEnd();
      };

      { //prepare backward offset and count 
	offset = 0;
	for(h=0;h<lnp_par.globH;h++) {
	  axSpike_backwardOffset[h] = offset;
	  count = pBackwardAxonalSpikes[h].expectedCount;
	  axSpike_backwardCount[h] = count; 
	  offset += count;
	};
      };

      {//debug section prints backward offset and count before transmission
	DPSNNverboseStart(false,thisTime_ms,0);
	for(h=0;h<lnp_par.globH;h++) {
	  printf(
		 "860- exch Ax Spik - %d ms from h=%d: ->targ h=%d: backCount=%d,backOff=%d\n", 
		 thisTime_ms, h, lnp_par.loc_h,
		 axSpike_backwardCount[h],
		 axSpike_backwardOffset[h]);
	  fflush(stdout);
	};
	DPSNNverboseEnd();
      };
      {//debug section prints forward offset and count before transmission
	DPSNNverboseStart(false,thisTime_ms,0);
	for(h=0;h<lnp_par.globH;h++) {
	  printf(
		 "864- exch Ax Spik - %d ms from h=%d: ->targ h=%d: forwCount=%d,forwOff=%d\n", 
		 thisTime_ms, lnp_par.loc_h, h, 
		 axSpike_forwardCount[h],
		 axSpike_forwardOffset[h]);
	  fflush(stdout);
	};
	DPSNNverboseEnd();
      };

      DPSNNverboseStart(false,thisTime_ms,0);
      {
	uint32_t sendPayloadCount,recPayloadCount;
	sendPayloadCount = axSpike_forwardOffset[lnp_par.globH-1]+
	  axSpike_forwardCount[lnp_par.globH-1];
	recPayloadCount = axSpike_backwardOffset[lnp_par.globH-1]+
	  axSpike_backwardCount[lnp_par.globH-1];
	if(lnp_par.loc_h==0) {
	  printf("863- exch Ax Spikspike - Payload to/from h (spikes unit) t=%d, h=%d, sendPayloadCount=%d, recPayloadCount=%d\n",
		 thisTime_ms, lnp_par.loc_h, sendPayloadCount,recPayloadCount);
	  fflush(stdout);
	}
      }
      DPSNNverboseEnd();
 
#ifdef MPIandDALenvironmentSelected //not DAL
      { 
	/* MPI_Comm comm; */
     
	int axSpike_datatype_numberOfBlocks;
	int axSpike_datatype_arrayOfBlocklengths[2];
	MPI_Aint axSpike_datatype_arrayOfDisplacements[2];
	MPI_Datatype axSpike_arrayOfDatatype[2];
	MPI_Datatype axSpike_datatype;

      /* comm=MPI_COMM_WORLD; */
      /*
      axSpike_datatype_numberOfBlocks=1;
      //2 int in the axonalSpikeDataOnlyClass
      axSpike_datatype_arrayOfBlocklengths[0]=4;
      axSpike_datatype_arrayOfDisplacements[0]=0;
      axSpike_arrayOfDatatype[0]=MPI_INT;
      */


	axSpike_datatype_numberOfBlocks=2;
	//1 double in the axonalSpikeDataOnlyClass (originalEmissionTime)
	axSpike_datatype_arrayOfBlocklengths[1]=1;
	axSpike_datatype_arrayOfDisplacements[1]=0;
	axSpike_arrayOfDatatype[1]=MPI_DOUBLE;
	//1 int in the axonalSpikeDataOnlyClass (neuron id)
	axSpike_datatype_arrayOfBlocklengths[0]=1;
	axSpike_datatype_arrayOfDisplacements[0]=8;
	axSpike_arrayOfDatatype[0]=MPI_INT;

	MPI_Type_create_struct(
			       axSpike_datatype_numberOfBlocks,
			       axSpike_datatype_arrayOfBlocklengths,
			       axSpike_datatype_arrayOfDisplacements,
			       axSpike_arrayOfDatatype,
			       &axSpike_datatype);
	MPI_Type_commit(&axSpike_datatype);

#define axSpikesPayloadUsesAllToAllv
#ifdef axSpikesPayloadUsesAllToAllv
	MPI_Alltoallv((void *) axSpike_forwardBuffer, 
		      (int *)  axSpike_forwardCount, 
		      (int *)  axSpike_forwardOffset,
		      axSpike_datatype,
		      (void *) axSpike_backwardBuffer,
		      (int *)  axSpike_backwardCount,
		      (int *)  axSpike_backwardOffset,
		      axSpike_datatype,
		      MPI_COMM_WORLD);

	/* char filename1[20]; */
	/* sprintf(filename1, "data_%u.stdout", lnp_par.loc_h); */

	/* save_and_tell(filename1, */
	/* 		    (void*)axSpike_forwardBuffer, (void*)axSpike_backwardBuffer, */
	/* 		    axSpike_forwardCount, axSpike_backwardCount, */
	/* 		    axSpike_forwardOffset, axSpike_backwardOffset, */
	/* 		    lnp_par.loc_h, */
	/* 		    lnp_par.globH, */
	/* 		    sizeof(typeof(axSpike_forwardBuffer[0]))); */

#else //axSpikesPayloadUsesAllToAllv
      // WARNING
      // The following code works for large configurations only if the
      // option "-mca btl_sm_use_cma 1" is used
	{
	  MPI_Status status;
	  int returnCount;
	  MPI_Request requestS[lnp_par.globH];
	  for(h=0;h<lnp_par.globH;h++) {
	    if(axSpike_forwardCount[h] != 0) {
	      MPI_Send((void *)&(axSpike_forwardBuffer[axSpike_forwardOffset[h]]), 
		       axSpike_forwardCount[h], axSpike_datatype, h, 
		       thisTime_ms, MPI_COMM_WORLD);
	      DPSNNverboseStart(false,1,0);
	      {printf("MPI Axonal Spikes Exchange - At t=%dms between procs %d ==> %d MPI_Send %d spikes\n",
		      thisTime_ms,lnp_par.loc_h,h,axSpike_forwardCount[h]);fflush(stdout);}
	      DPSNNverboseEnd();
	    }
	  }
	  //barrier();
	  for(h=0;h<lnp_par.globH;h++) {
	    if(axSpike_backwardCount[h] != 0) {
	      MPI_Recv((void *)&(axSpike_backwardBuffer[axSpike_backwardOffset[h]]), 
		       axSpike_backwardCount[h], axSpike_datatype, h, 
		       thisTime_ms, MPI_COMM_WORLD, &status);
	      DPSNNverboseStart(false,1,0);
	      MPI_Get_count(&status, axSpike_datatype, &returnCount);
	      {printf("MPI Axonal Spikes Exchange - At t=%dms between procs %d <== %d MPI_Recv %d spikes\n",
		      thisTime_ms,lnp_par.loc_h,h,returnCount);fflush(stdout);}
	      DPSNNverboseEnd();
	    }
	  }
	}
#endif //axSpikesPayloadUsesAllToAllv     

      };
#else //not MPIandDALenvironmentSelected is DAL
      {
	CREATEPORTVAR(O);
	CREATEPORTVAR(I);
	for (h = 0; h < lnp_par.globH;h++) {
	  if(axSpike_forwardCount[h] != 0) {
	    CREATEPORT(O, PORT_O, 1, h, lnp_par.globH);
	    DAL_write((void*) O, 
		      &(axSpike_forwardBuffer[axSpike_forwardOffset[h]]), 
		      axSpike_forwardCount[h] * sizeof(int) * 2, p);
	  };
	};
	for (h = 0; h < lnp_par.globH;h++) {
	  if(axSpike_backwardCount[h] != 0) {
	    CREATEPORT(I, PORT_I, 1, h, lnp_par.globH);
	    DAL_read((void*) I,
		     &(axSpike_backwardBuffer[axSpike_backwardOffset[h]]), 
		     axSpike_backwardCount[h] * sizeof(int) * 2, p);
	  };
	};
      }
#endif //not MPIandDALenvironmentSelected is DAL

      {//debug section prints backward offset and count after transmission
	DPSNNverboseStart(false,thisTime_ms,0);
	for(h=0;h<lnp_par.globH;h++) {
	  printf(
		 "870- exch Ax Spik - on h=%d rec4 ALLTOALLV: %d ms h=%d->targ h=%d: backCount=%d,backOff=%d\n", 
		 lnp_par.loc_h, thisTime_ms,  h, lnp_par.loc_h,
		 axSpike_backwardCount[h],
		 axSpike_backwardOffset[h]);
	  fflush(stdout);
	};
	DPSNNverboseEnd();
      };

      {//copy received data into final structures
	for(h=0;h<lnp_par.globH;h++) {
	  offset = axSpike_backwardOffset[h];
	  count = axSpike_backwardCount[h];

	  for(a=0;a< count; a++) {
	    pBackwardAxonalSpikes[h].list[a] = 
	      axSpike_backwardBuffer[offset + a];
	    pBackwardAxonalSpikes[h].count++;	   
	  }
	};
      };

      {//debug section print received data
	DPSNNverboseStart(false,thisTime_ms,0);
	for(h=0;h<lnp_par.globH;h++) {
	  offset = axSpike_backwardOffset[h];
	  count = axSpike_backwardCount[h];
	  for(a=0;a<count;a++) {
	    axonalSpikeDataOnlyClass axSp;
	    axSp = axSpike_backwardBuffer[offset + a];
	    printf(
		   "875- exch Ax Spik- rec5 ALLTOALLV: %d ms r_h=%d->h=%d: rec(glob_n=%d,orig %.7f ms)\n",
		   thisTime_ms, h, lnp_par.loc_h, 
		   axSp.pre_glob_n, axSp.originalEmissionTime);
	    if(!((axSpike_backwardBuffer[offset + a].originalEmissionTime >= thisTime_ms-1) &&
		 (axSpike_backwardBuffer[offset + a].originalEmissionTime < thisTime_ms))){
	      printf("ERROR wrong original emiss. time receiving ax spike MPI\n");
	      fflush(stdout);exit(0);
	    };
	  };
	};
	DPSNNverboseEnd();
      };

      DPSNNverboseStart(false,thisTime_ms,4);	 	   
      printf(
	     "895- MPC::exchAxSpikes MPI_Alltoallv completed on h=%03d\n", lnp_par.loc_h);
      fflush(stdout);
      DPSNNverboseEnd();

    } else { // pktLength true

      MPI_Alltoallv((void *) pForwardAxonalSpikes[0].list,
		    (int *)  axSpike_forwardCount,
		    (int *)  axSpike_forwardOffset, MPI_spike,
		    (void *) pBackwardAxonalSpikes[0].list,
		    (int *)  axSpike_backwardCount,
		    (int *)  axSpike_backwardOffset, MPI_spike,
		    MPI_COMM_WORLD);

      /* char filename6[20]; */
      /* sprintf(filename6, "dump_%u.post_ATAV", lnp_par.loc_h); */

      /* save_and_tell(filename6, */
      /* 		  (void*)pForwardAxonalSpikes[0].list, (void*)pBackwardAxonalSpikes[0].list, */
      /* 		  FC, BC, */
      /* 		  axSpike_forwardOffset, axSpike_backwardOffset, */
      /* 		  lnp_par.loc_h, */
      /* 		  lnp_par.globH, */
      /* 		  sizeof(typeof(pForwardAxonalSpikes[0].list[0]))); */

    } // pktLength

  } //endif not a single h

#else // APENET_RDMA defined
  {
    int n_expected_rEvents = 0 ;
    int n_expected_sEvents = 0;
    int n_rEvents = 0;
    int n_sEvents = 0;
    int retcode = 0;
    MPI_Status status;
    int returnCount;
    int sreq = 0;
    int rreq = 0;
    int nreqs = 0;
    MPI_Request requests[2*lnp_par.globH];



    for (h = 0; h < lnp_par.globH;h++) {
      if( h == lnp_par.loc_h) {
	DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	printf("ape this is me %d, just copy ms=%d from loc_h=%d to h=%d bytecount=%d bytes\n",(ape_descriptor->coords)[h].u, thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass));
	DPSNNverboseEnd();
	for(i=0; i < pBackwardAxonalSpikes[h].expectedCount; i++) {
	  pBackwardAxonalSpikes[h].list[i] = pForwardAxonalSpikes[h].list[i]; 
	  pBackwardAxonalSpikes[h].count = pBackwardAxonalSpikes[h].expectedCount;
	}
      }
    }
 
    for (nreqs = 0, h = (lnp_par.loc_h + 1)% lnp_par.globH; h != lnp_par.loc_h; h = (h+1) % lnp_par.globH, ++rreq) {
      if(pBackwardAxonalSpikes[h].expectedCount > 0 )
	{
	  if((ape_descriptor->coords)[h].u == (ape_descriptor->coords)[lnp_par.loc_h].u)
	    {
	      DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms); 
	      printf("ape MPI Recv ms=%d on loc_h=%d from h=%d bytecount=%d bytes\n",thisTime_ms, 
		     lnp_par.loc_h, h, pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass)); 
	      DPSNNverboseEnd(); 

	      MPI_Irecv((void *)& (pBackwardAxonalSpikes[h].list[0]),
			pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass),
			MPI_CHAR,
			h,
			thisTime_ms,
			MPI_COMM_WORLD,
			&(requests[nreqs]));

	      DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms); 
	      printf("ape MPI done Recv ms=%d on loc_h=%d from h=%d bytecount=%d bytes\n",thisTime_ms, 
		     lnp_par.loc_h, h, pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass)); 
	      DPSNNverboseEnd(); 

	      pBackwardAxonalSpikes[h].count = pBackwardAxonalSpikes[h].expectedCount;
	      nreqs++;
	    }
	  else
	    {
	      n_expected_rEvents++; 

	    }
	}
    }


    for (h = (lnp_par.loc_h + lnp_par.globH - 1)% lnp_par.globH; h != lnp_par.loc_h; h = (h+lnp_par.globH-1)% lnp_par.globH, ++sreq) {
      if(pForwardAxonalSpikes[h].count > 0 && ((ape_descriptor->coords)[h].u == (ape_descriptor->coords)[lnp_par.loc_h].u))
	{

	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("ape MPI Send ms=%d from loc_h=%d to h=%d bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass));
	  DPSNNverboseEnd();

	  MPI_Isend( (void *)& (pForwardAxonalSpikes[h].list[0]),
		    pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass),
		     MPI_CHAR,  h, thisTime_ms, MPI_COMM_WORLD, &(requests[nreqs]));
	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("ape MPI Send done ms=%d from loc_h=%d to h=%d bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass));
	  DPSNNverboseEnd();
	  nreqs++;
	}
    }

    if(nreqs)
      {
	MPI_Waitall(nreqs, requests, MPI_STATUSES_IGNORE);

      }

    DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
    printf("ape MPI done ms=%d sent=%d recv'd=%d loc_h=%d\n",thisTime_ms, sreq,rreq, lnp_par.loc_h);
    DPSNNverboseEnd();

/* 	DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms); */
/* 	printf("ape this is not me %d != %d  ms=%d loc_h=%d  h=%d bytecount=%d bytes\n", (ape_descriptor->coords)[lnp_par.loc_h].u, (ape_descriptor->coords)[h].u, thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass)); */
/* 	  DPSNNverboseEnd(); */
/* 	if(pForwardAxonalSpikes[h].count != 0  */
/* 	   && ((ape_descriptor->coords)[h].u == (ape_descriptor->coords)[lnp_par.loc_h].u)) { */
	 
/* 	  DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms); */
/* 	  printf("ape MPI Send  ms=%d from loc_h=%d to h=%d bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass)); */
/* 	  DPSNNverboseEnd(); */

/* 	  MPI_Send( (void *)& (pForwardAxonalSpikes[h].list[0]), */
/* 		    pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass),  */
/* 		    MPI_CHAR,  h, thisTime_ms, MPI_COMM_WORLD);  */
/* 	  DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms); */
/*     			printf("ape MPI Send done ms=%d from loc_h=%d to h=%d bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h, pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass)); */
/*     			DPSNNverboseEnd(); */
/* 	} */
	
/* 	if(pBackwardAxonalSpikes[h].expectedCount > 0) { */
/* 	  if(((ape_descriptor->coords)[h].u == (ape_descriptor->coords)[lnp_par.loc_h].u)) { */

/* 	    DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms); */
/* 	    printf("ape MPI Recv ms=%d on loc_h=%d from h=%d bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h, pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass)); */
/* 	    DPSNNverboseEnd(); */
	    
	    
/* 	    MPI_Recv((void *)& (pBackwardAxonalSpikes[h].list[0]),  */
/* 		     pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass),  */
/* 		     MPI_CHAR,  h, thisTime_ms, MPI_COMM_WORLD, &status); */

/* 	    DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms); */
/* 	    printf("ape MPI Recv done ms=%d on loc_h=%d from h=%d bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h, pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass)); */
/* 	    DPSNNverboseEnd(); */

/* 	    DPSNNverboseStart(true,1,0); */
/* 	    MPI_Get_count(&status, MPI_CHAR, &returnCount); */
/* 	    if(returnCount!=pBackwardAxonalSpikes[h].expectedCount*sizeof(axonalSpikeDataOnlyClass)){ */
/* 	      printf("ERROR: Axon Spikes ex: MPI_Get_count != 1\n"); */
/* 	      fflush(stdout); exit(0); */
/* 	    } */

/* 	    pBackwardAxonalSpikes[h].count = pBackwardAxonalSpikes[h].expectedCount; */
/* 	    DPSNNverboseEnd(); */
/* 	  } else { */
/* 	    n_expected_rEvents++; */
/* 	    DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms); */
/* 	    printf("ape REVENTS = %d h=%d from = %d\n",n_expected_rEvents,lnp_par.loc_h, h  ); */
/* 	    DPSNNverboseEnd(); */
/* 	  } */
/* 	} */


  

    
  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
  printf("ape before put cycle on h=%d from = %d ms=%d\n",lnp_par.loc_h, h, thisTime_ms );
  DPSNNverboseEnd();
    for (h = 0; h < lnp_par.globH;h++) {
      if(pForwardAxonalSpikes[h].count > 0) {
	if( h != lnp_par.loc_h && ((ape_descriptor->coords)[h].u != (ape_descriptor->coords)[lnp_par.loc_h].u)) {
	  for(i=0; i < pForwardAxonalSpikes[h].count; i++) {
	    ((axonalSpikeDataOnlyClass*)(apeForwardBuffer[h]))[i] = pForwardAxonalSpikes[h].list[i];
	  }
	
	  size_t bytecount = pForwardAxonalSpikes[h].count*sizeof(axonalSpikeDataOnlyClass);
	  if(bytecount%16) {
	    bytecount= ((bytecount >> 4) + 1)<<4;
	  }
	 
	  if(bytecount > APE_BUFFER_SIZE)
	    {
	    
	      printf("ape ERROR bytecount exceeds buffer size!! ms=%d loc_h=%d"
		     " bytecount=%d\n",thisTime_ms, lnp_par.loc_h, bytecount);
	      fflush(stdout); exit(0);
	    }
	  am_context_t p_send_ctx = (am_context_t)0xc1c1c1c100000000ULL;

	  retcode = am_put(ape_descriptor->link, 
			   (ape_descriptor->coords)[h],
			   ape_descriptor->ports[h],
			   thisTime_ms,  
			   AM_BLOCKING,
			   p_send_ctx,
			   (void*)(apeForwardBuffer[h]),
			   bytecount,
			   (am_vaddr_t)(apeRemoteBackwardBuffer[h]));

	  if(retcode != APELINK_NOERROR) {

	    printf("[%d] error %d in am_put to process %d, exiting...\n", 
		   lnp_par.loc_h, retcode, h);
	    fflush(stdout); exit(0);	     
	  }
	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("ape after am_put ms=%d from loc_h=%d to h=%d dst= %d %d %d"
		 " bytecount=%d bytes\n",thisTime_ms, lnp_par.loc_h, h,  
		 (ape_descriptor->coords)[h].s.x, (ape_descriptor->coords)[h].s.y, 
		 (ape_descriptor->coords)[h].s.z,  bytecount);
	  DPSNNverboseEnd();

	  n_expected_sEvents++;

	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("ape SEVENTS = %d h=%d from = %d\n",n_expected_sEvents,lnp_par.loc_h, h  );
	  DPSNNverboseEnd();
	}
      }
    };

   DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
    printf("ape entering in wait event cycle ms=%d on loc_h=%d n_expected_rEvents=%d n_expected_sEvents=%d\n",thisTime_ms, lnp_par.loc_h, n_expected_rEvents, n_expected_sEvents);
    DPSNNverboseEnd();  

  
    while((n_rEvents < n_expected_rEvents) || (n_sEvents < n_expected_sEvents)) {
      volatile am_event_t ev;
      retcode = am_wait_event(ape_descriptor->link,(am_event_t*) (&ev), 0 );
      
      if(retcode < 0) {
	if(retcode == APELINK_WOULDBLOCK) {
	  continue; 
	} else {
	  printf("ERROR in am_wait_event got retcode=%d, exiting\n", retcode);
	  fflush(stdout); exit(0);
	  break;
	}
      }
      
      switch(ev.ev_type){
      case AM_TYPE_RECV:
	{
	  volatile am_event_recv_t* recv = &ev.ev.ev_recv;
		  (ape_descriptor->rx_pkts)++;
	  n_rEvents++;	  
	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("ape received packet ms=%d on loc_h=%d data=%d coord=%d %d %d"
		 " port=%d nrecv=%d total_recvd=%d FirstElementEmissionTime=%d"
		 " ms neuron=%d buf=%p buflen=%d size=%d\n",
		 thisTime_ms, lnp_par.loc_h, recv->data, recv->src_addr.s.x, 
		 recv->src_addr.s.y, recv->src_addr.s.z, 
		 recv->src_port, n_rEvents, ape_descriptor->rx_pkts,  
		 ((axonalSpikeDataOnlyClass*)recv->buf)[0].originalEmissionTime, 
		 ((axonalSpikeDataOnlyClass*)recv->buf)[0].pre_glob_n,
		 (axonalSpikeDataOnlyClass*)recv->buf, recv->buflen, recv->len);
	  DPSNNverboseEnd();

	  
	  if(thisTime_ms != ((axonalSpikeDataOnlyClass*)recv->buf)[0].originalEmissionTime) 
	    {
	      DPSNNverboseStart(true,thisTime_ms,lnp_par.debugPrintEnable_ms);
	      printf("ape AAAAHHHHHHH ms=%d on loc_h=%d FirstElementEmissionTime=%d ms"
		     " neuron=%d vaddr=%p buf=%p\n",thisTime_ms, lnp_par.loc_h,   
		     ((axonalSpikeDataOnlyClass*)recv->buf)[0].originalEmissionTime,  
		     ((axonalSpikeDataOnlyClass*)recv->buf)[0].pre_glob_n, 
		     (char*)recv->vaddr, (char*)recv->buf);
	      
	      DPSNNverboseEnd();

	      printf("ape dump buffer ms=%d\n",thisTime_ms );
	      for(int g = 0; g < AM_MAX_MSG_SIZE/sizeof(axonalSpikeDataOnlyClass); g++)
		{
		  
		  printf("%d %d - ",  ((axonalSpikeDataOnlyClass*)recv->buf)[g].originalEmissionTime, ((axonalSpikeDataOnlyClass*)recv->buf)[g].pre_glob_n);
		}
	      printf("ape end buffer ms=%d\n",thisTime_ms );

	      fflush(stdout);
	      sleep(100);
	      exit(0);
	      
	    }

	  break;
	}
      case AM_TYPE_SENT:
	{
	  volatile am_event_sent_t* sent = &ev.ev.ev_sent;
	  (ape_descriptor->tx_pkts)++;
	  n_sEvents++;
	  DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
	  printf("ape sent packet ms=%d on h_loc=%d data=%d nsent=%d"
		 " total_sent=%d dest_addr=%d %d %d port=%d size=%d vaddr=%p\n",
		 thisTime_ms, lnp_par.loc_h, sent->data, n_sEvents, 
		 ape_descriptor->tx_pkts, sent->dst_addr.s.x, sent->dst_addr.s.y, 
		 sent->dst_addr.s.z, sent->dst_port, sent->len, sent->rem_vaddr );
	  DPSNNverboseEnd();  

	  break;
	}
      default:
	{
	  printf("ape ERROR unknown event type ms=%d on loc_h=%d\n",
		 thisTime_ms, lnp_par.loc_h);
	}
      }
    }

    DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
    printf("ape after of wait event cycle ms=%d on loc_h=%d\n",thisTime_ms, lnp_par.loc_h);
    DPSNNverboseEnd(); 
 
    for (h = 0; h < lnp_par.globH;h++) {
      if(pBackwardAxonalSpikes[h].expectedCount > 0)
	{
	  if(h != lnp_par.loc_h) {
	    if( ((ape_descriptor->coords)[h].u != (ape_descriptor->coords)[lnp_par.loc_h].u)){
	       for(i=0; i < pBackwardAxonalSpikes[h].expectedCount; i++) {
		 pBackwardAxonalSpikes[h].list[i] = ((axonalSpikeDataOnlyClass*)(apeBackwardBuffer[h]))[i];
		 pBackwardAxonalSpikes[h].count = pBackwardAxonalSpikes[h].expectedCount;
	       }
	     }
	  }
	}
    }

   if(thisTime_ms == (lnp_par.totalSimTime_ms-1))
    {
    	for(h = 0; h < lnp_par.globH ;h++)
    	{
	  if((ape_descriptor->coords)[lnp_par.loc_h].u == (ape_descriptor->coords)[h].u)
	      {

		DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
    		printf("ape unregister buf %p for h=%d ms=%d on h=%d\n", 
		       apeBackwardBuffer[h],h, thisTime_ms, lnp_par.loc_h);
    		DPSNNverboseEnd();

	      } else {
    		am_vaddr_t vaddr;
    		am_unregister_buf(ape_descriptor->link,
    				(am_vaddr_t)(apeBackwardBuffer[h]),
    				(void*)(apeBackwardBuffer[h]),
    				APE_BUFFER_SIZE,
    				AM_PERSISTENT_BUF);
    		DPSNNverboseStart(false,thisTime_ms,lnp_par.debugPrintEnable_ms);
    		printf("ape unregister buf %p ms=%d on h=%d\n", apeBackwardBuffer[h],
		       thisTime_ms, lnp_par.loc_h);
    		DPSNNverboseEnd();

    		free((void*)apeBackwardBuffer[h]);
    		free((void*)apeForwardBuffer[h]);
	  }
    	}
    	APE_RDMA_Finalize(ape_descriptor);
    }
  }
#endif // APENET_RDMA defined

  for (h = 0; h < lnp_par.globH;h++) {
    if(pBackwardAxonalSpikes[h].count != pBackwardAxonalSpikes[h].expectedCount) {
      printf("ERROR in exAxSpikes, expected %d spikes, received %d \n",
	   pBackwardAxonalSpikes[h].expectedCount, pBackwardAxonalSpikes[h].count);
      fflush(stdout);exit(0);
    }
  };

  DPSNNverboseStart(false,thisTime_ms,4);
    for (h = 0; h < lnp_par.globH;h++) {
      printf(
      "898- %d ms h=%03d ->targ %03d: %d spikes received - exchAxSpikes\n", 
	     thisTime_ms, h, lnp_par.loc_h, pBackwardAxonalSpikes[h].count);
      fflush(stdout); 
    };
  DPSNNverboseEnd();

  for(h=0; h<lnp_par.globH; h++) {
    if(h != lnp_par.loc_h)
      if(pForwardAxonalSpikes[h].count != 0)
	spikePayloadSize.accumulateMsgSize(pForwardAxonalSpikes[h].count);
  }
};





