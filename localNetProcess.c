// localNetProcess.c
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
//#include <time.h>
//#include <unistd.h>

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_getParameters.h"

#include "DPSNN_localNet.h"
#include "DPSNN_stopWatch.h"
#include "DPSNN_stat.h"
#include "DPSNN_debug.h"

#include "localNetProcess.h"
    
#ifdef MPIandDALenvironmentSelected
    extern int MPIglobH;
    extern int MPIprocessRank;
    extern localNetProcess_localState lnp_ls;
//acronym lnp_ls_par for localNetProcess_localState_parameters
#else
  #ifdef DALonlyEnvironmentSelected
    #define lnp_ls p->local 
  #endif
#endif

void localNetProcess_init(DALProcess *p)
{
   p->local->id=GETINDEX(0);
   sprintf(p->local->name, "C0- localNetProcess_%d", p->local->id);

   DPSNNverboseStart(false,1,0);
      printf("010- %s: init() - started\n",p->local->name);
      fflush(stdout);
   DPSNNverboseEnd();

   //get from the environment some parameters for this process 
   lnp_ls_par = getParameters(p->local->id);

   //program the total simulated time
     //note - 1 simulated day = 
     // 60*60*24 simulated seconds == 86 400 000 ms
   lnp_ls_stopWatch.programSimTime_ms(lnp_ls_par.totalSimTime_ms);
 
   DPSNNverboseStart(true,1,0);
   if(lnp_ls_par.loc_h == 0) {
      printf("localNetProcess-%s: init() - globH = %d\n",
	     p->local->name,lnp_ls_par.globH);
      fflush(stdout);
   };
   DPSNNverboseEnd();

   //printing the parameters of this run
   lnp_ls_stat.prep(lnp_ls_par);

   if(lnp_ls_par.globH<=4) {
     if(lnp_ls_par.loc_h==0) {
        lnp_ls_stat.parameters.openFile(0, lnp_ls_par.moduloSec);
        lnp_ls_stat.parameters.write(0, lnp_ls_par.moduloSec);
        lnp_ls_stat.parameters.closeFile(0, lnp_ls_par.moduloSec);
      }
   };
   
   if((lnp_ls_par.loc_h<8) || (lnp_ls_par.loc_h == (lnp_ls_par.globH-1))) {
     lnp_ls_stat.spikingRates.openFile(0, lnp_ls_par.moduloSec);
   };

   lnp_ls_messagePassing.init(lnp_ls_par, p, 0);

   lnp_ls_neuralNet.clearAllChronometers();

   lnp_ls_neuralNet.init(lnp_ls_par, 
			 &(lnp_ls_messagePassing), &(lnp_ls_stopWatch),
			 &(lnp_ls_stat));
   lnp_ls_neuralNet.printAllInitChronoResults();

   p->local->fireCount=0;

   DPSNNverboseStart(true,1,0);
   if(lnp_ls_par.loc_h == 0) {

      printf("500- %s: init() - completed - from now it will fire()\n",
	   p->local->name);
      fflush(stdout);
   };
   DPSNNverboseEnd();
}

int localNetProcess_fire(DALProcess *p)
{ 
  stopWatchStatusEnum stopWatchStatus;
  int thisSimTime_ms;

  thisSimTime_ms = lnp_ls_stopWatch.getSimTime_ms();
  stopWatchStatus = lnp_ls_stopWatch.getStatus();
 
  DPSNNverboseStart(true,1,0);
  if(lnp_ls_par.loc_h == 0) {
    if(thisSimTime_ms == 0) {
      printf("501- h=000 - sim step %d ms - START\n",
	   thisSimTime_ms);
    }
  };
  DPSNNverboseEnd();
 
  //next line tests if end of simulation reached
  if(stopWatchStatus != maxSimTime_ms_reached) {

    //simulation time step of this localNet
    lnp_ls_neuralNet.completeTimeStep_ms();

    //advance the stopWatch of this localNet
    lnp_ls_stopWatch.advanceSimTime_ms();

   DPSNNverboseStart(false,1,0);
    if(lnp_ls_par.loc_h==0) { 
      printf("950- h=000 - sim step %d ms - END\n",
	   thisSimTime_ms);
      fflush(stdout);
    };
   DPSNNverboseEnd();

    //normal return
    return(0);  
  } else {
   if((lnp_ls_par.loc_h<8) || (lnp_ls_par.loc_h == (lnp_ls_par.globH-1))) {
    lnp_ls_stat.spikingRates.closeFile(thisSimTime_ms/1000, lnp_ls_par.moduloSec);
   };

   DPSNNverboseStart(false,1,0);
    if(lnp_ls_par.loc_h==0) {
      printf("530- h=000 - maxSimTime_ms reached at %d ms\n",
	   thisSimTime_ms);
      fflush(stdout); 
    };
   DPSNNverboseEnd();
  
    //sim closing actions
    if(lnp_ls_par.loc_h == 0) {
       DAL_send_event(EVENT_stop_fsm, p);
    };
    #warning "check DAl fire return convention"
    return(1);
  };
};

void localNetProcess_finish(DALProcess *p) {
  DPSNNverboseStart(false,1,0);
	printf("%s: finish() - started\n",p->local->name);
	fflush(stdout);
	fflush(stdout);
  DPSNNverboseEnd();
  if(lnp_ls_par.chrono==1) {
    if(lnp_ls_par.loc_h <2 || (lnp_ls_par.loc_h >= (lnp_ls_par.globH-2))) { 
      lnp_ls_neuralNet.printAllSimulChronoResults();
    }
    lnp_ls_neuralNet.printAllStatChronoResults();
    lnp_ls_neuralNet.printStatFiringsInChronoWindow(lnp_ls_par.loc_h);
  };
}
