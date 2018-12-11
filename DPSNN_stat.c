// DPSNN_stat.c
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

#include "DPSNN_parameters.h"
#include "DPSNN_stat.h"
#include "DPSNN_debug.h"
#include "DPSNN_chrono.h"
#include "DPSNN_dataStructDims.h"

int get_syn_pre_glob_n(const int i, const int j);
float get_syn_value(const int i, const int j);

int totalStat_fopen;
int max_totalStat_fopen;

void statBasicClass::test_prepCalled() {
  if(lnp_par.loc_h==0) { 
    if(prepCalled != true) { 
      printf("ERROR stat_prep... not called\n");
      fflush(stdout); exit(0);} 
  } 
};

void statBasicClass::prep(
  const struct DPSNN_parameters lnp_par_initValue, 
  const char prefixName_value[]) 
{
  if(writeDisable==false) {
    prepCalled = true;
    lnp_par=lnp_par_initValue;
    strcpy(prefixName, prefixName_value);
  }
};

void statBasicClass::test_fopenStatus(const char methodName[],
const bool stat_booleanValue, const int sec, const int moduloSec) {
    if(fopenStatus != stat_booleanValue) {
      printf("ERROR on h=%d stat fopenStatus wrong, called by %s at sec=%d moduloSec=%d\n", 
	     lnp_par.loc_h, methodName, sec, moduloSec);
      fflush(stdout);
      exit(0); }
};

void statBasicClass::openFile(const uint sec, const int moduloSec) 
{
  double totalSimTime;
  totalSimTime = atoi(getenv("env_totalSimTime_ms"));
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
      if(sec % moduloSec == 0) {
        test_prepCalled();
        test_fopenStatus("openFile", false,sec, moduloSec);
        if(lnp_par.howManyOutputFiles==longOutput_1){
          sprintf(fName, "%s-sec%d-C%d-H%d-h%03d.dat",
	      prefixName,sec,
	      lnp_par.globCFT,lnp_par.globH,lnp_par.loc_h);
        } else {
	  if(sec==0) {
	    sprintf(fName, "%s-start-C%d-H%d-h%d.dat",
		  prefixName,
		  lnp_par.globCFT,lnp_par.globH,lnp_par.loc_h);
	  //} else if(sec==((totalSimTime/1000)-2)){
	  } else if(sec % 2 == 0){
	    sprintf(fName, "%s-lasteven-C%d-H%d-h%d.dat",
		  prefixName,
		  lnp_par.globCFT,lnp_par.globH,lnp_par.loc_h);
	  } else {
	    sprintf(fName, "%s-lastodd-C%d-H%d-h%d.dat",
		  prefixName,
		  lnp_par.globCFT,lnp_par.globH,lnp_par.loc_h);
	  };
	};
      }

      DPSNNverboseStart(false,1,0);
        printf("openFile: %s at sec %d on loc_h %d\n",
	       fName, sec, lnp_par.loc_h);
        fflush(stdout);
      DPSNNverboseEnd();
	//3 because stdin, stdout, stderr are there
	if( (totalStat_fopen + 3) <= DSD__FOPEN_MAX) {
        fp = fopen(fName,"w");
        if(fp == NULL) {printf("ERROR NULL fp\n"); fflush(stdout); exit(0);};
	fopenStatus=true;
        totalStat_fopen++;
        if(totalStat_fopen > max_totalStat_fopen) {
          max_totalStat_fopen = totalStat_fopen; 
        } 
      }else{
        printf("ERROR too much fopen while opening %s DSD__FOPEN_MAX=%d\n",
	       fName,DSD__FOPEN_MAX);fflush(stdout);exit(0);
      }    
    };
  };
};  

void statBasicClass::closeFile(const uint sec, const int moduloSec) 
{
  int status = 1;
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
      if(sec % moduloSec == 0) {
        debugPrintStatStatus(sec, moduloSec);
        test_prepCalled();
        test_fopenStatus("closeFile", true, sec, moduloSec);
        totalStat_fopen--; 
#ifndef DALonlyEnvironmentSelected
        if(totalStat_fopen<3) { 
          printf("ERROR totalStat_fopen < 3 in %s\n",prefixName);
          fflush(stdout); exit(0);
        };
#endif
        status = fclose(fp);
        DPSNNverboseStart(false,1,0);
          printf("closeFile: %s \n",fName);
	  fflush(stdout);
        DPSNNverboseEnd();
        if(status!=0) { \
          printf("ERROR fclose!=0 in stat_sfclose in %s\n",prefixName);
          fflush(stdout); exit(0);
        };
        fopenStatus=false;
      };
    };
  };
};

void statBasicClass::debugPrintStatStatus(const int sec, const int moduloSec) 
{
  DPSNNverboseStart(false,1,0);
    printf("\n begin status report for statistic function\n");
    printf("at sec=%d, moduloSec=%d writeDisable=%d, prepCalled=%d, fopenStatus=%d, ", 
	   sec, moduloSec, writeDisable, prepCalled, fopenStatus);
    if(prepCalled==true) {printf("prefixName=%s, ",prefixName);};
    if(fopenStatus==true) {printf("fName=%s, ",fName);};
    printf("\n end\n");
    fflush(stdout);
  DPSNNverboseEnd();
}

void statParametersClass::write(const uint32_t sec, const int moduloSec) {
  int i;
  FILE * locFp;

  printf("statParameters.write called on lnp_par.loc_h=%03d\n", lnp_par.loc_h);

    if(writeDisable==false) {
        if(sec == 0) {
          test_prepCalled();
          test_fopenStatus("write", true, sec, moduloSec);    
        for(i=0;i<2;i++) {
          if(i==0) { locFp = stdout;
          } else   { locFp = fp;}
          fprintf(locFp, "\n011- code rev        : %s\n",lnp_par.codeRev); fflush(locFp);

          fprintf(locFp, "012- neurons            : %4d=globN %4d=locN %4d=globNe  %4d=globNi\n", 
	      lnp_par.globN, lnp_par.locN, lnp_par.globNe, lnp_par.globNi); fflush(locFp);
          fprintf(locFp, "013- synapses           :      %4d=D    %2.2f=M/D\n", 
	       lnp_par.D, float(lnp_par.M)/float(lnp_par.D));fflush(locFp);
          fprintf(locFp, "014- distribution       : %4d=loc_h %4d=globH\n", 
	      lnp_par.loc_h, lnp_par.globH);fflush(locFp);
          fprintf(locFp, "016- programmed sim time:  tot=%8d ms  moduloSec=%3d\n",
		lnp_par.totalSimTime_ms,lnp_par.moduloSec);fflush(locFp);
	  fprintf(locFp, "016a- start and stop SynPeriodicProbe_s: %d %d\n",
		  lnp_par.startSynPeriodicProbe_s, lnp_par.stopSynPeriodicProbe_s);fflush(locFp);
          if(lnp_par.chrono==0) {
	    fprintf(locFp, "017- chrono performance :  not active\n"); fflush(locFp);
          } else if (lnp_par.chrono==1) { 
	    fprintf(locFp, "017- chrono performance:  ACTIVATED\n"); fflush(locFp);
          } else { fprintf(locFp, "ERROR chrono env value unrecognized\n");
	  fflush(locFp);exit(0);};
      };
    };
  };
};

double presentTime;

void statSpikingRatesClass::write(const uint32_t sec, 
  const int moduloSec, 
  const int N_firings,
  const int N_thalamicInputs) 
{
  int writeCount;
  double previousTime,deltaTime;
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
  
      previousTime = presentTime;
      presentTime = MPI_Wtime();
      if (sec > 0)
        deltaTime = presentTime - previousTime;
      else
        deltaTime = 0;

      test_prepCalled();
      test_fopenStatus("write",true, sec, moduloSec);
      writeCount =
	fprintf(fp,"sec=%d, firings=%06d, rates: tot %03.3f elapsedTime=%.3f sec\n",
	  sec, N_firings,float(N_firings)/float(lnp_par.locN),deltaTime);
      fflush(fp);
      if(writeCount<0) { 
          printf("ERROR fprintf returns negative in SpikingRates\n"); 
          fflush(stdout); 
          exit(0);
      };     
    };
  };
};

void statSpikingRatesPerPopClass::write(const uint32_t sec,
					const uint32_t milliSec, 
					const uint32_t moduloSec, 
					const uint32_t ratesSampling, 
					uint32_t *N_firingsPerPop,
					const uint32_t *neuSubPopCount,
					const uint32_t neuSubPopTotal) 
{
  int writeCount;
  uint32_t subPop,subPopTot;
    
  if(writeDisable==false) {
    if( (sec >= lnp_par.startStatFiles_s) &&
        (sec <= lnp_par.stopStatFiles_s ) ) {

      if((sec % moduloSec) == 0) {
	if ((sec >= lnp_par.startRatesPerPop_s)&&
             (sec <= lnp_par.stopRatesPerPop_s)   ) {
          if((milliSec % ratesSampling) == 0) { 
	    test_prepCalled();
	    test_fopenStatus("write",true,sec,moduloSec);
	    fprintf(fp,"%d ",milliSec-ratesSampling);
      
	    for(subPopTot=0; subPopTot<lnp_par.locCFT*neuSubPopTotal;
	      subPopTot+=neuSubPopTotal){
	      for(subPop=0;subPop<neuSubPopTotal;subPop++){
	        writeCount =
	          fprintf(fp,"%.7g ",((float)(N_firingsPerPop[subPopTot+subPop])*1000.0) /
			 ((float)(neuSubPopCount[subPop])*(float)(ratesSampling)));
	      }
	    }

	    fprintf(fp," \n");
	    fflush(fp);
	    if(writeCount<0) { 
	      printf("ERROR fprintf returns negative in SpikingRatesPerPop\n"); 
	      fflush(stdout); 
	      exit(0);
	    }

	    //Clear array of firings per population for the next usage
	    for(uint32_t i=0;i<lnp_par.globCFT*neuSubPopTotal;i++)
	      N_firingsPerPop[i] = 0;
	  }
        }
      }
    };
  };
};

void statSpikesClass::write(const uint32_t sec, const int moduloSec,
			    const int N_firings,
			    const int thisMs,
			    const int *data,
			    const double *emissionTime)
{
  int i;
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
      if((sec % moduloSec) == 0) { 
        test_prepCalled();
        test_fopenStatus("write", true, sec, moduloSec);
        for (i=0;i < N_firings;i++) 
        {
  	  fprintf(fp, "%d %.7g\n", data[i], emissionTime[i]);
          fflush(fp);
        };
      };
    };
  };		  
};


void statSpikeCountPerMsClass::write(const uint32_t sec, const int moduloSec,
				     const int thisMs,
				     const int numberOfSpikesInMs)
{
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
      if((sec % moduloSec) == 0) { 
        test_prepCalled();
        test_fopenStatus("write", true, sec, moduloSec);
        fprintf(fp, "%07d  %6d\n", thisMs, numberOfSpikesInMs);
        fflush(fp);
      };
    };
  };		  
};

void statThalamicInputClass::write(const uint32_t sec, const int moduloSec,
				   const int thisMs,
				   const int numberOfTargetNeu)
{
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
      if((sec % moduloSec) == 0) { 
        test_prepCalled();
        test_fopenStatus("write", true, sec, moduloSec);
        fprintf(fp, "%07d  %6d\n", thisMs, numberOfTargetNeu);
        fflush(fp);
      };
    };
  };		  
};

void statFloatEvery_msClass::write(const uint32_t sec, const int moduloSec,
				   const uint32_t t_ms,
				   const int neuron_id,
				   const float value) 
{
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
       (sec <= lnp_par.stopStatFiles_s)){
      if((sec % moduloSec) == 0) {
        test_prepCalled();
        test_fopenStatus("write", true, sec, moduloSec);
        if( lnp_par.debugPrintEnable_ms == 0 ||
	    ((t_ms >= lnp_par.debugPrintEnable_ms) &&
             (t_ms <= (lnp_par.debugPrintEnable_ms + 20) ))) 
        {
            fprintf(fp, "%05d %03d %04d %03.2f\n",sec, t_ms, neuron_id, value);
            fflush(fp);
	};
      };   
    };
  };
};

void stat3FloatsEvery_msClass::write(const uint32_t sec, const int moduloSec,
				     const uint32_t t_ms,
				     const int neuron_id,
				     const float value0,
				     const float value1,
				     const float value2)
{
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
      (sec <= lnp_par.stopStatFiles_s)){
      if((sec % moduloSec) == 0) {
        test_prepCalled();
        test_fopenStatus("write", true, sec, moduloSec);
        if( lnp_par.debugPrintEnable_ms == 0 ||
	    ((t_ms >= lnp_par.debugPrintEnable_ms) &&
             (t_ms <= (lnp_par.debugPrintEnable_ms + 20) ))) 
        {
          fprintf(fp, "%05d %03u %04d %03.2f %03.2f %03.2f\n",
          sec, t_ms%1000, neuron_id, value0, value1, value2);
          fflush(fp);
        };  
      };
    };
  };
};

void stat_LIFCAivu_Class::write(const uint32_t sec, const int moduloSec,
				const uint32_t t_ms,
				const double currentTime,
				const int neuron_id,
				const double value0,
				const double value1,
				const double value2)
{
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
        (sec <= lnp_par.stopStatFiles_s)) {
      test_prepCalled();
      test_fopenStatus("write", true, sec, moduloSec);
      if( lnp_par.debugPrintEnable_ms == 0 ||
	    ((t_ms >= lnp_par.debugPrintEnable_ms) &&
             (t_ms <= (lnp_par.debugPrintEnable_ms + 20) ))) 
      {
          fprintf(fp, "%05u %03.9f %04d %03.7f %03.7f %03.7f \n",
		  sec, currentTime, neuron_id, value0, value1, value2);
          fflush(fp);
      };  
    };
  };
};

void statSynClass::write(const uint32_t sec, const int moduloSec,
			 const int this_ms, const int synIndex,
			 const int preSynNeu, const int postSynNeu,
			 const int synDelay,
			 const float synWeight, const float synDeriv)
{
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
        (sec <= lnp_par.stopStatFiles_s)) {
      if((sec % moduloSec) == 0) {
          test_prepCalled();
          test_fopenStatus("write", true, sec, moduloSec);
          fprintf(fp,"%04d %04d %02d %2.4f %2.4f %05d %03d\n", 
  	      preSynNeu, postSynNeu,
	      synDelay, synWeight, synDeriv,  
	      sec, this_ms % 1000);
          fflush(fp);
      };
    };
  };
};

void statSynPeriodicProbeClass::write(const uint32_t sec,
				      const int moduloSec,
				      const synapseClass synToBeReported)
{
  instrumentedSynapse instrumentedSynapseDummy;
  if(writeDisable==false) {
    if ((sec >= lnp_par.startStatFiles_s)&&
        (sec <= lnp_par.stopStatFiles_s)) {
      if((sec % moduloSec) == 0) {
        if ((sec >= lnp_par.startSynPeriodicProbe_s)&&
            (sec <= lnp_par.stopSynPeriodicProbe_s)   ){
          test_prepCalled();
          test_fopenStatus("write", true, sec, moduloSec);
          instrumentedSynapseDummy.report(synToBeReported,fp,
	    lnp_par.locN,lnp_par.factorWeightType_2_Float);
          fflush(fp);
        };
      };
    };
  };
};

void statSTDPeventClass::write(const uint32_t sec, const int moduloSec,
			       const uint32_t t_ms,
			       const int pre_glob_n, const int post_glob_n,
			       const float contribute)
{
  if(writeDisable==false) {
  if ((sec >= lnp_par.startStatFiles_s)&&
      (sec <= lnp_par.stopStatFiles_s)) {
     if((sec % moduloSec) == 0) {
       if ((sec >= lnp_par.startSTDPevent_s)&&
          (sec <= lnp_par.stopSTDPevent_s)   ) {

         test_prepCalled();
         test_fopenStatus("write", true, sec, moduloSec);
         if( lnp_par.debugPrintEnable_ms == 0 ||
	    ((t_ms >= lnp_par.debugPrintEnable_ms) &&
             (t_ms <= (lnp_par.debugPrintEnable_ms + 20) ))) 
         {
          fprintf(fp,"%8u %7d %7d %2.7f \n",
		  t_ms, pre_glob_n, post_glob_n, contribute);
          fflush(fp);
         };
       };
      };
    };
  };
};


void statMessageTrafficClass::write(const uint32_t sec, const int moduloSec,
				    const uint32_t t_ms, uint32_t *forwardAxonalSpikesCount)
{
  uint32_t target_h;
  if(writeDisable==false) {
  if ((sec >= lnp_par.startStatFiles_s)&&
      (sec <= lnp_par.stopStatFiles_s)) {
      if((sec % moduloSec) == 0) {
        test_prepCalled();
        test_fopenStatus("write", true, sec, moduloSec);
        if( lnp_par.debugPrintEnable_ms == 0 ||
	  ((t_ms >= lnp_par.debugPrintEnable_ms) &&
	   (t_ms <= (lnp_par.debugPrintEnable_ms + 20) ))) {
	  fprintf(fp,"%04u",t_ms);
	  for(target_h=0;target_h<lnp_par.globH;target_h++) 
	    fprintf(fp," %6d",forwardAxonalSpikesCount[target_h]);
	  fprintf(fp,"\n");
	  fflush(fp);
	};
      };
    };
  };
};

void statClass::prep(struct DPSNN_parameters lnp_par_initValue) {
  lnp_par=lnp_par_initValue;

  //printf("on this system DSD__FOPEN_MAX value is %d\n",DSD__FOPEN_MAX);
  totalStat_fopen = 3; //stdin, stdout, stderr included
  max_totalStat_fopen = totalStat_fopen;
  fflush(stdout);
  
  //set the set of active statistic collectors
  spikes.writeDisable = lnp_par.spikes_writeDisable;
  spikePerMs.writeDisable = true;

  parameters.writeDisable = true;
  DPSNNverboseStart(false,1,0);
    if(lnp_par.globH <= 4) {
      parameters.writeDisable = false;
    };
  DPSNNverboseEnd();

  syn.writeDisable = true;
  DPSNNverboseStart(false,1,0);
    if(lnp_par.loc_h == 0) {
      syn.writeDisable = false;
    };
  DPSNNverboseEnd();

  spikingRates.writeDisable = true;
  DPSNNverboseStart(false,1,0);
  if(lnp_par.loc_h < 8  || (lnp_par.loc_h == (lnp_par.globH-1))) {
    spikingRates.writeDisable = false;
  };
  DPSNNverboseEnd();
  
  spikingRatesPerPop.writeDisable = lnp_par.statRatesPerPop_writeDisable;

  thalamicInput.writeDisable = true;
  inputCurrent.writeDisable = true;
  ivu0.writeDisable = true;
  ivu1.writeDisable = true;
  ivu2.writeDisable = true;
  ivu999.writeDisable = true;
  
  statSTDPevent.writeDisable = true;  
  DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h == 0) {
        statSTDPevent.writeDisable =lnp_par.statSTDPevent_writeDisable;   
    };
  DPSNNverboseEnd();

  statSynPeriodicProbeA.writeDisable = true;  
  statSynPeriodicProbeB.writeDisable = true;  
  DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h == 0) {
        statSynPeriodicProbeA.writeDisable =
	  lnp_par.statSynPeriodicProbe_writeDisable;   
        statSynPeriodicProbeB.writeDisable =
	  lnp_par.statSynPeriodicProbe_writeDisable;   
    };
  DPSNNverboseEnd();

  LIFCAivu0.writeDisable = true;
  LIFCAivu1.writeDisable = true;
  LIFCAivu2.writeDisable = true;
  LIFCAivu999.writeDisable = true;
  messageTraffic.writeDisable = true;

  //and calls the statistic preparation functions
  //if writeDisable==true nothing is done
  parameters.prep(lnp_par,"params");
  spikingRates.prep(lnp_par,"rates");
  spikingRatesPerPop.prep(lnp_par,"ratesPerPop");
  spikes.prep(lnp_par,"spikes");
  spikePerMs.prep(lnp_par,"spikePerMs");
  thalamicInput.prep(lnp_par,"thalInput");
  inputCurrent.prep(lnp_par,"inputCurrent");
  ivu0.prep(lnp_par,"ivu0");
  ivu1.prep(lnp_par,"ivu1");
  ivu2.prep(lnp_par,"ivu2");
  ivu999.prep(lnp_par,"ivu999");
  syn.prep(lnp_par,"syn");
  statSTDPevent.prep(lnp_par,"STDPevent");
  statSynPeriodicProbeA.prep(lnp_par, "synPeriodicProbe_A");
  statSynPeriodicProbeB.prep(lnp_par, "synPeriodicProbe_B");
  LIFCAivu0.prep(lnp_par,"LIFCAivu0");
  LIFCAivu1.prep(lnp_par,"LIFCAivu1");
  LIFCAivu2.prep(lnp_par,"LIFCAivu2");
  LIFCAivu999.prep(lnp_par,"LIFCAivu999");
  messageTraffic.prep(lnp_par,"messageTraffic");
};

void statClass::openFiles_forThisSecond(const uint32_t sec, const int moduloSec) {
  currentStat_sec=sec;


  //where writeDisable is set in statClass:prep()
  //nothing is done 
  if ((sec >= lnp_par.startStatFiles_s)&&
      (sec <= lnp_par.stopStatFiles_s)){

    if((sec % moduloSec)==0) {

      DPSNNverboseStart(false,1,0);
        printf(
        "openFiles_forThisSecond at sec=%d moduloSec=%d START on loc_h=%d\n",
        sec,moduloSec,lnp_par.loc_h);
        fflush(stdout);
      DPSNNverboseEnd();
      if(spikes.writeDisable == false)
        spikes.openFile(sec,moduloSec);
      if(spikePerMs.writeDisable == false)
        spikePerMs.openFile(sec,moduloSec);

      if ((sec >= lnp_par.startRatesPerPop_s)&&
          (sec <= lnp_par.stopRatesPerPop_s)   ) {
        if(spikingRatesPerPop.writeDisable == false)
           spikingRatesPerPop.openFile(sec,moduloSec);
      };
      
      if(thalamicInput.writeDisable == false)
        thalamicInput.openFile(sec,moduloSec);
      if(inputCurrent.writeDisable == false)
        inputCurrent.openFile(sec,moduloSec);
      if(ivu0.writeDisable == false)
        ivu0.openFile(sec,moduloSec);
      if(ivu1.writeDisable == false)
        ivu1.openFile(sec,moduloSec);
      if(ivu2.writeDisable == false)
        ivu2.openFile(sec,moduloSec);
      if(ivu999.writeDisable == false)
        ivu999.openFile(sec,moduloSec);
      if(syn.writeDisable == false)
        syn.openFile(sec,moduloSec);
      
      if ((sec >= lnp_par.startSTDPevent_s)&&
          (sec <= lnp_par.stopSTDPevent_s)   ) {
        if(statSTDPevent.writeDisable == false)
          statSTDPevent.openFile(sec,moduloSec);
      };
      
      if ((sec >= lnp_par.startSynPeriodicProbe_s)&&
          (sec <= lnp_par.stopSynPeriodicProbe_s)   ) {
        if(statSynPeriodicProbeA.writeDisable == false)
            statSynPeriodicProbeA.openFile(sec,moduloSec);
	if(statSynPeriodicProbeB.writeDisable == false)
            statSynPeriodicProbeB.openFile(sec,moduloSec);
      };
      if(LIFCAivu0.writeDisable == false)
        LIFCAivu0.openFile(sec,moduloSec);
      if(LIFCAivu1.writeDisable == false)
        LIFCAivu1.openFile(sec,moduloSec);
      if(LIFCAivu2.writeDisable == false)
        LIFCAivu2.openFile(sec,moduloSec);
      if(LIFCAivu999.writeDisable == false)
        LIFCAivu999.openFile(sec,moduloSec);
      if(messageTraffic.writeDisable == false)
        messageTraffic.openFile(sec,moduloSec);

      DPSNNverboseStart(false,1,0);
        printf(
          "openFiles_forThisSecond at sec=%d moduloSec=%d END on loc_h=%d\n",
          sec,moduloSec,lnp_par.loc_h);
        fflush(stdout);
      DPSNNverboseEnd();
    };
  };
};

void statClass::closeFiles_forThisSecond(const uint32_t sec, const int moduloSec) {

  if(currentStat_sec != sec) {
    printf("ERROR statClass::closeFiles_forThisSecond sec%d!=currentStat_sec%d\n",
	   currentStat_sec,sec); fflush(stdout);exit(0);}

  DPSNNverboseStart(false,1,0);
    printf("closeFiles_forThisSecond at sec=%d START on loc_h=%d\n",
      sec,lnp_par.loc_h);
    fflush(stdout);
  DPSNNverboseEnd();

  //where writeDisable is set in statClass:prep() 
  //nothing is done
  if ((sec >= lnp_par.startStatFiles_s)&&
      (sec <= lnp_par.stopStatFiles_s)){

    if((sec % moduloSec)==0) {
      if(spikes.writeDisable == false)
        spikes.closeFile(sec,moduloSec);
      if(spikePerMs.writeDisable == false)
        spikePerMs.closeFile(sec,moduloSec);

      if ((sec >= lnp_par.startRatesPerPop_s)&&
          (sec <= lnp_par.stopRatesPerPop_s)   ) {
         if(spikingRatesPerPop.writeDisable == false)
           spikingRatesPerPop.closeFile(sec,moduloSec);
      };
      if(thalamicInput.writeDisable == false)
        thalamicInput.closeFile(sec,moduloSec);
      if(inputCurrent.writeDisable == false)
        inputCurrent.closeFile(sec,moduloSec);
      if(ivu0.writeDisable == false)
        ivu0.closeFile(sec,moduloSec);
      if(ivu1.writeDisable == false)
        ivu1.closeFile(sec,moduloSec);
      if(ivu2.writeDisable == false)
        ivu2.closeFile(sec,moduloSec);
      if(ivu999.writeDisable == false)
        ivu999.closeFile(sec,moduloSec);
      if(syn.writeDisable == false)
        syn.closeFile(sec,moduloSec);
      
      if ((sec >= lnp_par.startSTDPevent_s)&&
          (sec <= lnp_par.stopSTDPevent_s)   ) {
        if(statSTDPevent.writeDisable == false)
          statSTDPevent.closeFile(sec,moduloSec);
      };
      
      if ((sec >= lnp_par.startSynPeriodicProbe_s)&&
          (sec <= lnp_par.stopSynPeriodicProbe_s)   ) {
        if(statSynPeriodicProbeA.writeDisable == false)
          statSynPeriodicProbeA.closeFile(sec,moduloSec);
        if(statSynPeriodicProbeB.writeDisable == false)
          statSynPeriodicProbeB.closeFile(sec,moduloSec);
      };
      
      if(LIFCAivu0.writeDisable == false)
        LIFCAivu0.closeFile(sec,moduloSec);
      if(LIFCAivu1.writeDisable == false)
        LIFCAivu1.closeFile(sec,moduloSec);
      if(LIFCAivu2.writeDisable == false)
        LIFCAivu2.closeFile(sec,moduloSec);
      if(LIFCAivu999.writeDisable == false)
        LIFCAivu999.closeFile(sec,moduloSec);
      if(messageTraffic.writeDisable == false)
        messageTraffic.closeFile(sec,moduloSec);
    }
  }
  DPSNNverboseStart(false,1,0);
    printf("closeFiles_forThisSecond at sec=%d DONE on loc_h=%d\n",
	   sec, lnp_par.loc_h);
    fflush(stdout);
  DPSNNverboseEnd();
};

void statClass::check_fileClosure(const int sec, const int moduloSec) {
  DPSNNverboseStart(false,1,0);  
  printf("checking file closure max_totalStat_fopen = %d should be 3 on loc_h=%d\n", 
	 max_totalStat_fopen,lnp_par.loc_h);
  fflush(stdout);
  DPSNNverboseEnd();
  DPSNNverboseStart(true,1,0);  
  if(totalStat_fopen != 3) {
    printf("ERROR in stat check_fclose totalStat_fopen!=3 at %d sec %d moduloSec\n",sec, moduloSec); 
    fflush(stdout); exit(0); };
  DPSNNverboseEnd();
};

