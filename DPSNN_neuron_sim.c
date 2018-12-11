// DPSNN_neuron_sim.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

// Saves spiking data each second in file spikes.dat
// To plot spikes, use MATLAB code: load spikes.dat;plot(spikes(:,1),spikes(:,2),'.');
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "DPSNN_debug.h"
#include "DPSNN_random.h"

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_synapse.h"
#include "DPSNN_neuron.h"
#include "DPSNN_stat.h"


#if defined(makeActiveLTD) || defined (makeActiveLTP)
  // exponential decay of Long Term Potentiation/Depression
  //only for synapses of excitatory neurons 
  // ...  Tau=20ms
  //LTD[0]=0.12;
  //LTD[t+1]=LTD[t]*0.95;
  //LTP[1]=0.1;
  //LTP[t+1]=LTP[t]*0.95;

//#define DSD__maxLongTermPlasticity_ms 128
//float LTPtable[DSD__maxLongTermPlasticity_ms];
//float LTDtable[DSD__maxLongTermPlasticity_ms];

#define DSD__maxSTDPspan_ms 128
float causalSTDPtable[DSD__maxSTDPspan_ms];
float antiCausalSTDPtable[DSD__maxSTDPspan_ms];


void neuronClass::initLongTermPlasticity() {
  int i;

  if(lnp_par.plasticityAlgo!=SongSTDPenum) {
    if(lnp_par.loc_h==0 && loc_n==0) {
      printf("ERROR in initLongTermPlasticity, Unknown plasticityAlgo requested\n");
      fflush(stdout);
      exit(0);
    }
  }

  if(lnp_par.loc_h==0 && loc_n==0) {
    if(lnp_par.SongSTDP_param.tauSTDP_ms*3>DSD__maxSTDPspan_ms) {
      printf("initLongTermPlasticity POSSIBLE ERROR: 3*tauSTDP>DSD__maxSTDPspan_ms\n");    fflush(stdout); exit(0); };
    if(lnp_par.SongSTDP_param.tauSTDP_ms <=0) {
      printf("initLongTermPlasticity ERROR: tauSTDP <= 0\n");
      fflush(stdout); exit(0); };
    if(lnp_par.SongSTDP_param.STDPmultiplier < 0) {
      printf("initLongTermPlasticity ERROR: STDPmultiplier < 0\n");
      fflush(stdout); exit(0);};
  }
  
  DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h==0 && loc_n==0) {
      printf("initLongTermPlasticity - SongSTDP plasticity selected - loc_h=0, loc_n=0\n");
      fflush(stdout);
    }
  DPSNNverboseEnd();
  
  causalSTDPtable[0] = lnp_par.SongSTDP_param.causalSTDPmaxAbsValue *
    lnp_par.SongSTDP_param.STDPmultiplier;
  causalSTDPtable[1] = causalSTDPtable[0] ;
  //antiCausalSTDPtable[0] = antiCausalSTDPsign * antiCausalSTDPmaxAbsValue;
  antiCausalSTDPtable[0] = lnp_par.SongSTDP_param.antiCausalCoeff *
                           lnp_par.SongSTDP_param.antiCausalSTDPmaxAbsValue *
                           lnp_par.SongSTDP_param.STDPmultiplier;
  
  for(i=1; i < DSD__maxSTDPspan_ms;i++) {    
    //LTPtable[i] = LTPtable[i-1] * 0.95;
    causalSTDPtable[i] = (1 - 1/lnp_par.SongSTDP_param.tauSTDP_ms)
      * causalSTDPtable[i-1];
  }
  for(i=1; i < DSD__maxSTDPspan_ms;i++) {
    //LTDtable[i] = LTDtable[i-1] * 0.95;
    antiCausalSTDPtable[i] = (1 - 1/lnp_par.SongSTDP_param.tauSTDP_ms)
      * antiCausalSTDPtable[i-1];
  }
};

float neuronClass::causalSTDP_ms(
  const int32_t localNeuralSpikeTime_ms,
  const int32_t previousSynActivationTime_ms) 
{
  //float fLTP_ms;
  //float fTimeDiff_ms;
  int timeDiff_ms;

  DPSNNverboseStart(true,1,0);
  if(localNeuralSpikeTime_ms <= previousSynActivationTime_ms) {
    printf(
  "ERROR causalSTDP time, locNeuSpike_ms %d should be later than synActiv_ms %d\n",
	   localNeuralSpikeTime_ms, previousSynActivationTime_ms);
    fflush(stdout);exit(0);
  };
  DPSNNverboseEnd();

  //  fTimeDiff_ms = float(
  //  localNeuralSpikeTime_ms - 1 - previousSynActivationTime_ms);
  //  fLTP_ms = 0.1 * exp(-fTimeDiff_ms/20.0);
  //if(previousSynActivationTime_ms != 0) {
    //excludes time before 0
  //   return(fLTP_ms);
  //}else{
  //  return(0.0);
  //}

  timeDiff_ms = localNeuralSpikeTime_ms - 1 -
                previousSynActivationTime_ms;

  if(previousSynActivationTime_ms != 0 ) {
    // excludes time before 0
    if(timeDiff_ms < DSD__maxSTDPspan_ms) {
      return(causalSTDPtable[timeDiff_ms]);
    }else{
      return(0.0);
    };
  }else{
    return(0.0);
  }
};

int LTPdebugCount;
int LTDdebugCount;

void neuronClass::causalSTDP_ms_ofBackwardSynapses(const uint32_t localSpikeTime_ms) {
  uint32_t synOffset;
  int32_t lastSynActivationTime_ms;
  synapseClass * pointSynapse;
  float fCausalSTDP;

  //LTP only for excitatory neurons

  for(synOffset=0;synOffset<N_pre;synOffset++) {
    pointSynapse = &(pointBackwardSynList[backwardSynOffset[synOffset]]);
    lastSynActivationTime_ms = pointSynapse->lastActivationTime;
    
    fCausalSTDP = causalSTDP_ms(localSpikeTime_ms, lastSynActivationTime_ms);

    DPSNNverboseStart(true,localSpikeTime_ms,lnp_par.debugPrintEnable_ms);
      if(fCausalSTDP != 0.0) {
        if(localSpikeTime_ms != 0.0) { 
          if(lnp_par.loc_h==0) {
	    if((pointSynapse->pre_glob_n < lnp_par.locN) &&
	        pointSynapse->post_glob_n < lnp_par.locN) {
                  pStat->statSTDPevent.write(localSpikeTime_ms/1000, 
		       lnp_par.moduloSec, localSpikeTime_ms,
		       pointSynapse->pre_glob_n, pointSynapse->post_glob_n, 
		       fCausalSTDP);
	    };
	  };
	};
      };
    DPSNNverboseEnd();

    #ifdef makeActiveLTP
      // (weightType)(pointSynapse->timeDerivative + (int16_t)(fLTP_ms*DSD__factorFloat_2_weightType));
      pointSynapse->timeDerivative = pointSynapse->timeDerivative + fCausalSTDP;
    //  };
    #else
      #warning "LTP disactivated: see makefile"
    #endif
  };
};

float neuronClass::antiCausalSTDP_ms(const int32_t synActivationT_ms) 
{
  float fAntiCausalSTDP_ms;
  //float fTimeDiff_ms;
  int timeDiff_ms;
  
  DPSNNverboseStart(true,0,1);
  if(synActivationT_ms < lastNeuralEmittedSpikeTime_ms) {
    printf("ERROR in antiCausalSTDP_ms synActivationT_ms %d < lastNeuralSpikeTime_ms %d\n",
	   synActivationT_ms, lastNeuralEmittedSpikeTime_ms);
    fflush(stdout);exit(0);
  };
  DPSNNverboseEnd();

  //fTimeDiff_ms = float( synActivationT_ms  -  
  //			lastNeuralEmittedSpikeTime_ms);
  //fLTD_ms = - 0.12 * exp(-fTimeDiff_ms/20.0);

  timeDiff_ms = synActivationT_ms - lastNeuralEmittedSpikeTime_ms;
  if(timeDiff_ms < DSD__maxSTDPspan_ms) {
    fAntiCausalSTDP_ms = antiCausalSTDPtable[timeDiff_ms];
  }else{
    return(0.0);
  }

  return(fAntiCausalSTDP_ms); 
};

void neuronClass::synapseSetLastActiveTime_antiCausalSTDP(
			    synapseClass *pointerActiveSynapse,
			    const uint32_t thisTimeStep_ms) {
  float antiCausalSTDPthis_ms;

  #ifdef makeActiveLTD

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
#ifdef LIFCAneuron
  if(lnp_par.loc_h==0) {
    printf("synapseSetLastActiveTime_antiCausalSTDP-START on neu=%d pre_neu=%d at %d ms\n",
	   glob_n, pointerActiveSynapse->pre_glob_n, thisTimeStep_ms);
  };
#else
    if(lnp_par.loc_h==0) {
    printf("synapseSetLastActiveTime_antiCausalSTDP-START on neu=%d I=%f v=%2.2f u=%2.2f addSynCurr at %d ms\n",
      glob_n, InputCurrent, v, u, thisTimeStep_ms);
    };
#endif
  DPSNNverboseEnd();


  pointerActiveSynapse->lastActivationTime = thisTimeStep_ms;
  antiCausalSTDPthis_ms = antiCausalSTDP_ms(thisTimeStep_ms);

  DPSNNverboseStart(true,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
    if(antiCausalSTDPthis_ms != 0.0) { 
      if(lnp_par.loc_h==0) {
	if((pointerActiveSynapse->pre_glob_n < lnp_par.locN) &&
	   pointerActiveSynapse->post_glob_n < lnp_par.locN) {
             pStat->statSTDPevent.write(pointerActiveSynapse->lastActivationTime/1000, 
		       lnp_par.moduloSec,
		       pointerActiveSynapse->lastActivationTime,
		       pointerActiveSynapse->pre_glob_n,
		       pointerActiveSynapse->post_glob_n, 
		       antiCausalSTDPthis_ms);
	};
      };
    };
  DPSNNverboseEnd();

      //  (weightType)(pointerActiveSynapse->timeDerivative + (int16_t)(LTDthis_ms*DSD__factorFloat_2_weightType));
      pointerActiveSynapse->timeDerivative = 
	pointerActiveSynapse->timeDerivative + antiCausalSTDPthis_ms;
  //};

    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
      printf("synapseSetLastActiveTime_antiCausalSTDP-EXISTS CURRENT on n=%d: totI=%f antiCausalSTDPms=%f at %d ms\n",
	     glob_n, InputCurrent, antiCausalSTDPthis_ms, thisTimeStep_ms);
    DPSNNverboseEnd();

    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
      printf("synapseSetLastActiveTime_antiCausalSTDP-C: called on neu=%d at %d ms\n",
	     glob_n, thisTimeStep_ms);
    DPSNNverboseEnd();  

    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
#ifdef LIFCAneuron
      printf("synapseSetLastActiveTime_antiCausalSTDP-END on neu=%d I=%.12f, v=%.12f, c=%.12f at %d ms\n",
	     glob_n, InputCurrent, v, c, thisTimeStep_ms);
#else
      printf("synapseSetLastActiveTime_antiCausalSTDP-END on neu=%d I=%.12f, v=%.12f, u=%.12f at %d ms\n",
	     glob_n, InputCurrent, v, u, thisTimeStep_ms);
#endif
    DPSNNverboseEnd();  

  #else
    #warning "LTD disactivated: see makefile"
  #endif

};
#endif //end of #if defined(makeActiveLTD) || defined (makeActiveLTP)


#ifdef LIFCAneuron
void neuronClass::clearInputCurrent() {
  inputCurrentPlus = 0.0;
  inputCurrentMinus = 0.0;
  inputCurrentPlusCount  = 0;
  inputCurrentMinusCount = 0;
  inputCurrentsCount = 0;
};
void neuronClass::addInputSpike(const float current,const double emissionTime) {
  double tempInt;
  DPSNNverboseStart(true,1,0); 
  if(inputCurrentsCount >= DSD__maxSimultaneousSpikesOnSameTarget){
    printf("WARNING: on h=%d on neuron %d the number of sinaptic input currents exceeded the maximum number allowed, fixed to the value miniSimSpikes=%d. Enlarge this value in the Makefile to allow more sinaptic currents.\n",lnp_par.loc_h,glob_n,DSD__maxSimultaneousSpikesOnSameTarget);
  }
  DPSNNverboseEnd();
  inputCurrents[inputCurrentsCount].inputCurrent = current;
  inputCurrents[inputCurrentsCount].originalEmissionTime = modf(emissionTime,&tempInt);
  //inputCurrents[inputCurrentsCount].originalEmissionTime = emissionTime;
  inputCurrentsCount++;
}
void neuronClass::addInputCurrent(const float inputValue) {
  float oldIPlus;
  float oldIMinus;

  if(inputValue>=0){
    inputCurrentPlus += inputValue;
    oldIPlus = inputCurrentPlus;
    inputCurrentPlusCount++;
    DPSNNverboseStart(false,1,lnp_par.debugPrintEnable_ms);
    printf(
      "on neuron %d adding %f to old IPlus=%f, new IPlus=%f for a total of %d positive current values in H=%d\n",
	   glob_n, inputValue, oldIPlus, inputCurrentPlus,inputCurrentPlusCount,lnp_par.loc_h);
    DPSNNverboseEnd();

  }
  else{
    inputCurrentMinus += inputValue;
    oldIMinus = inputCurrentMinus;
    inputCurrentMinusCount++;
    DPSNNverboseStart(false,1,lnp_par.debugPrintEnable_ms);
    printf(
      "on neuron %d adding %f to old IMinus=%f, new IMinus=%f for a total of %d positive current values in H=%d\n",
	   glob_n, inputValue, oldIMinus, inputCurrentMinus,inputCurrentMinusCount,lnp_par.loc_h);
    DPSNNverboseEnd();

  }
};

void neuronClass::sortCurrentsInPlace(inputCurrentArrayStruct *array, uint32_t count)
 {
   double tmp1,tmp2;
   uint32_t i,j;
      
   for (j = 0; j < count-1; j++ )
     { 
         for (i = 0; i < count-1-j; i++)
         {
           if (array[i].time>array[i+1].time) 
           { 
             tmp1 = array[i].time; 
             tmp2 = array[i].value; 
             array[i].time = array[i+1].time; 
             array[i].value = array[i+1].value; 
             array[i+1].time = tmp1;
             array[i+1].value = tmp2;
           } 
         }
     }
   
 }

void neuronClass::mergeCurrentsArrays(inputCurrentArrayStruct *aOut, inputCurrentArrayStruct *a1, inputCurrentArrayStruct *a2, uint32_t c1, uint32_t c2)
{
  uint32_t i, j, k;
  uint32_t cOut;

  j = k = 0;
  cOut = c1 + c2;
 
  for (i = 0; i < cOut;) {
    if (j < c1 && k < c2) {
      if (a1[j].time < a2[k].time) {
        aOut[i] = a1[j];
        j++;
      }
      else {
        aOut[i] = a2[k];
        k++;
      }
      i++;
    }
    else if (j == c1) {
      for (; i < cOut;) {
        aOut[i] = a2[k];
        k++;
        i++;
      }
    }
    else {
      for (; i < cOut;) {
        aOut[i] = a1[j];
        j++;
        i++;
      }
    }
  }
  
  return;
}

void neuronClass::mergeFunction(inputCurrentArrayStruct *A, uint32_t p, uint32_t q, uint32_t r) {
  uint32_t i, j, k;
  inputCurrentArrayStruct B[r+1];

  i = p;
  j = q+1;
  k = 0;
  while (i<=q && j<=r) {
    if (A[i].time<=A[j].time) {
      B[k] = A[i];
      i++;
    } else {
      B[k] = A[j];
      j++;
    }
    k++;
  }
  while (i<=q) {
    B[k] = A[i];
    i++;
    k++;
  }
  while (j<=r) {
    B[k] = A[j];
    j++;
    k++;
  }
  for (k=p; k<=r; k++)
    A[k] = B[k-p];
  return;
}

void neuronClass::mergeSort(inputCurrentArrayStruct *A, uint32_t p, uint32_t r) 
{
  uint32_t q;

  if (p<r) {
    q = (p+r)/2;
    mergeSort(A, p, q);
    mergeSort(A, q+1, r);
    mergeFunction(A, p, q, r);
  }
  return;
}

/* Function to merge the two haves arr[l..m] and arr[m+1..r] of array arr[] */
static inline void iter_merge(inputCurrentArrayStruct *arr, int l, int m, int r)
{
  int i, j, k;
  int n1 = m - l + 1;
  int n2 =  r - m;

  /* create temp arrays */
  inputCurrentArrayStruct L[n1], R[n2+1];

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++)
    L[i] = arr[l + i];
  /* memcpy(L, arr+l, n1*sizeof(inputCurrentArrayStruct)); */
  for (j = 0; j < n2; j++)
    R[j] = arr[m + 1+ j];
  /* memcpy(R, arr+m+1, n2*sizeof(inputCurrentArrayStruct)); */

  /* Merge the temp arrays back into arr[l..r]*/
  i = 0;
  j = 0;
  k = l;

  while (i < n1 && j < n2) {
    if (L[i].time <= R[j].time) {
      arr[k] = L[i];
      i++;
    } else {
      arr[k] = R[j];
      j++;
    }
    k++;
  }

  /* Copy the remaining elements of L[], if there are any */
  while (i < n1) {
    arr[k] = L[i];
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there are any */
  while (j < n2) {
    arr[k] = R[j];
    j++;
    k++;
  }
}

// Utility function to find minimum of two integers
static inline int min(int x, int y) { return (x<y) ? x : y; }

static void iter_mergeSort(inputCurrentArrayStruct *arr, int n)
{
  int curr_size;  // For current size of subarrays to be merged
                       // curr_size varies from 1 to n/2
  int left_start; // For picking starting index of left subarray
                       // to be merged
  
  // Merge subarrays in bottom up manner.  First merge subarrays of
  // size 1 to create sorted subarrays of size 2, then merge subarrays
  // of size 2 to create sorted subarrays of size 4, and so on.
  for (curr_size=1; curr_size<=n-1; curr_size = 2*curr_size) {

    // Pick starting point of different subarrays of current size
    for (left_start=0; left_start<n-1; left_start += 2*curr_size) {

      // Find ending point of left subarray. mid+1 is starting 
      // point of right
      int mid = left_start + curr_size - 1;

      int right_end = min(left_start + 2*curr_size - 1, n-1);
      /* if (right_end == mid) { */
      /* 	printf("cioccato n2 nullo: mid = %d, left_start = %d, curr_size = %d, n = %d\n", mid, left_start, curr_size, n); */
      /* 	printf("dump:\n"); */
      /* 	for (int i = 0; i<n; i++) printf("[%d]: [%f,%f] ", i, arr[i].time, arr[i].value); */
      /* 	printf("proc %u\n", h); */
      /* 	exit(0); */
      /* } */

      // Merge Subarrays arr[left_start...mid] & arr[mid+1...right_end]
      iter_merge(arr, left_start, mid, right_end);
    }
  }
}

void neuronClass::regularDynamicTimeStep(
  const uint32_t thisTimeStep_ms, const float *pInputDistrib, uint32_t * pThalInputCounter) {
    
  // double assignedInternalSpikeTime;
  //double internalInputCurrent;
  //double externalInputCurrent;
  //double thalamIn;
  double thisTimeStepF; 
  double nextTimeStepF;
  /* bool existInternalCurrent; */
  /* bool existExternalCurrent; */
  //double singleStep;
  //double sampleTime;
  //double thisSampleTime;
  inputCurrentArrayStruct inputCurrentArray[inputCurrentsCount+maxExternCurrentCount];
  //inputCurrentArrayStruct inputCurrentArray[inputCurrentsCount];
  //inputCurrentArrayStruct inputStimulusArray[maxExternCurrentCount];
  //inputCurrentArrayStruct inputToNeuronArray[inputCurrentsCount+maxExternCurrentCount];
  //double singleInputCurrent;
  uint32_t i,totalCurrentItems,inputStimulusCount;
  /* double randVal; */

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if (glob_n==168) 
    printf("START of current funnel on neu=%d at %u ms \n",glob_n,thisTimeStep_ms);
  DPSNNverboseEnd();

#undef reproducibleCode
  // Enable the following line to ensure code reproducibility over more than 1 process
  //#define reproducibleCode
#ifdef reproducibleCode
  uint32_t timeStepModule;
  uint32_t thisTab;
  uint32_t localSeed;

  if ((glob_n % lnp_par.tabNeuron) == 0){
    timeStepModule = 1 << 21;
    thisTab = (uint32_t)(glob_n / lnp_par.tabNeuron);
    localSeed = (0xAFDE75 + (thisTab + thisTimeStep_ms * DSD__maxTabNumber) % DSD__INT_MAX + 
	    (thisTimeStep_ms >> timeStepModule)) % DSD__INT_MAX;
    pLocalNetRandDev->SetRandomSeed(localSeed);
    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
      printf("TabSeed: on H=%d at ms=%u neu=%d timeStepModule=%d thisTab=%d seed=%d \n",
	     lnp_par.loc_h,thisTimeStep_ms,glob_n,timeStepModule,thisTab,localSeed);
    DPSNNverboseEnd();
  }
#endif

 // local step of LIFCA neuron dynamic
  thisTimeStepF = (double)thisTimeStep_ms;
  nextTimeStepF = (double)(thisTimeStepF + 1.0);
  /* existInternalCurrent = false; */
  /* existExternalCurrent = false; */

  totalCurrentItems=0;
  
  for(i=0;i<inputCurrentsCount;i++){
    inputCurrentArray[i].time = thisTimeStepF+inputCurrents[i].originalEmissionTime;
    inputCurrentArray[i].value = inputCurrents[i].inputCurrent;
  }

    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
    if (glob_n==168)
      if(inputCurrentsCount){
	for(i=0;i<inputCurrentsCount;i++)
	  printf("Real Input Current nosort: at %u ms, on neu %d: inputCurrent.value=%f inputCurrent.time=%f \n",
		 thisTimeStep_ms,glob_n,inputCurrentArray[i].value,inputCurrentArray[i].time);}
    DPSNNverboseEnd(); 

    if(inputCurrentsCount)
      //sortCurrentsInPlace(inputCurrentArray,inputCurrentsCount);
      /* iter_mergeSort(inputCurrentArray,inputCurrentsCount); */
      mergeSort(inputCurrentArray,0,inputCurrentsCount-1);

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
  if (glob_n==168)
    if(inputCurrentsCount){
      for(i=0;i<inputCurrentsCount;i++)
	printf("Real Input Current sorted: at %u ms, on neu %d: inputCurrent.value=%f inputCurrent.time=%f \n",
	       thisTimeStep_ms,glob_n,inputCurrentArray[i].value,inputCurrentArray[i].time);}
  DPSNNverboseEnd(); 
 
  inputStimulusCount=0;
if (invNuExt>0){
  if (assignedExternalSpikeTime < nextTimeStepF) {
    /* existExternalCurrent=true; */
    do 
      {
	inputCurrentArray[inputCurrentsCount+inputStimulusCount].time = assignedExternalSpikeTime;	
	inputCurrentArray[inputCurrentsCount+inputStimulusCount].value = 
	                  pTableLUT_synBath[uint32_t(pLocalNetRandDev->Random()*ANALOG_DEPTH)];
	DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
	if (glob_n==168)
	  printf("Bath Input Current: at %u ms, on neu %d: stimulusCurrent.value=%f stimulusCurrent.time=%f \n",
		 thisTimeStep_ms,glob_n,inputCurrentArray[inputCurrentsCount+inputStimulusCount].value,
		 inputCurrentArray[inputCurrentsCount+inputStimulusCount].time);
	DPSNNverboseEnd();  
	inputStimulusCount++;
	assignedExternalSpikeTime += getPoissonLapse();
      } while ((assignedExternalSpikeTime < nextTimeStepF) && 
	       (inputStimulusCount < maxExternCurrentCount)); //End of while

    if(inputStimulusCount >= maxExternCurrentCount)
      if(assignedExternalSpikeTime < nextTimeStepF){
	do
	  {
	    assignedExternalSpikeTime += getPoissonLapse();
	  } while (assignedExternalSpikeTime < nextTimeStepF);
	DPSNNverboseStart(true,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
	printf("WARNING: on h=%d at ms %u on neuron %d the number of generated bath currents exceeded the maximum number allowed, fixed to the value maxExternCurrentCount=%d. Enlarge this value to allow more external currents.\n",lnp_par.loc_h,thisTimeStep_ms,glob_n,maxExternCurrentCount);
	DPSNNverboseEnd();
      }
    
    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
    if (glob_n==168)
      printf("assignedExternalSpikeTime=%f nextTimeStepF=%f inputStimulusCount=%d maxExternCurrentCount=%d\n",assignedExternalSpikeTime,nextTimeStepF,inputStimulusCount,maxExternCurrentCount);
    DPSNNverboseEnd();  
  }
} //End of if (invNuExt>0)

    totalCurrentItems = inputCurrentsCount + inputStimulusCount;
      
    if(inputCurrentsCount && inputStimulusCount)
      //mergeCurrentsArrays(inputToNeuronArray, inputCurrentArray, inputStimulusArray, inputCurrentsCount, inputStimulusCount);
      mergeFunction(inputCurrentArray, 0, inputCurrentsCount-1, totalCurrentItems-1);
    DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
    if (glob_n==168)
      for(i=0;i<totalCurrentItems;i++)
	printf("Total input currents: at %u ms, on neu %d: inputCurrent.value=%f inputCurrent.time=%f \n",
	       thisTimeStep_ms,glob_n,inputCurrentArray[i].value,inputCurrentArray[i].time);
    DPSNNverboseEnd();  
    
    for(i=0;i<totalCurrentItems;i++)
      LIFCADynamic(thisTimeStep_ms,inputCurrentArray[i].time,inputCurrentArray[i].value);

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if (glob_n==168) 
    printf("END of current funnel and neu dynamic on neu=%d at %u ms \n",glob_n,thisTimeStep_ms);
  DPSNNverboseEnd();

};

void neuronClass::LIFCADynamic(const uint32_t thisTimeStep_ms, const double currentTime, const double thisInputCurrent) {
  double t;
  double deltaT;
  double TFLES; // Time From Last Emitted Spike
  double rc;
  /* double c0, rm, erm, erc; */
  double TrTe;
  double deltaTHalf, rmHalf;
  
 // local step of LIFCA neuron dynamic
  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
  if (glob_n==2424)
    printf("DYN START neu=%d, v=%f, c=%f, I=%f at %u ms with currentTime=%f \n",
	   glob_n, v, c, thisInputCurrent, thisTimeStep_ms, currentTime);
  DPSNNverboseEnd();
  
  DPSNNverboseStart(true,1,0);   
  if(((currentTime - (double)thisTimeStep_ms )>= 1.0 ) || ((currentTime - (double)thisTimeStep_ms)< 0 )) {
    printf("ERROR in neu dynamic: input currentTime in wrong ms: currentTime=%f, thisTimeStep_ms=%u on neu %d\n",
	   currentTime,thisTimeStep_ms,glob_n);fflush(stdout);exit(0);
  }
  DPSNNverboseEnd();

  t = currentTime;
  
  deltaT = t - Tr;    //Tr=Arriving time of the last pre-synaptic spike
  TFLES = t - Te;     //Te=Emission time of the last spike
  
  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms); 
  if (glob_n==2424)
    printf("LIFCA DYN: thisTimeStep_ms=%u currentTime=%f deltaT=%f TFLES=%f\n",
	 thisTimeStep_ms,t,deltaT,TFLES);
  DPSNNverboseEnd();

  if (TFLES > Tarp) {   //Updates the neuron state outside absolute refractory period (ARP)

    //when exiting refractory period, the dynamic is divided in two sub=times:
    //1- first period between previous incoming spike and end of refractory period: only c dynamic
    //2- second period between end of refractory perid and current incoming spike: normal dynamic
    TrTe = Tr - Te;
    if (TrTe < Tarp) {  //1- first period
      c *= exp(-(deltaT - TFLES + Tarp) / TauC);
      deltaT = TFLES - Tarp;
    }   //Now follow 2- second period with normal dynamic

    //Normal dynamic

#ifdef PerseoDyn
    c0 = c;

    rc = -deltaT / TauC;
    rm = -deltaT / Tau;

    erm = exp(rm);
    erc = exp(rc);

    v = v * erm - gC * (TauC * Tau) / (TauC - Tau) * c0 * (erc - erm);
    c *= erc;

#else

    // EPA - PSP : 27 July 2016
    // We assume that the "true" equation is:
    // v' = -v / Tau - gC * c / Capacity + I / Capacity
    // c' = -c / TauC
    // where Capacity is fixed to 1
    // In this implementation we integrated the ODE 
    // using a second order approximation of the analitic solution (exponential decay) for "c"
    // a second order approximation in deltaT/2 for "v"

    deltaTHalf = deltaT / 2;
    rmHalf = deltaTHalf / Tau;
    rc = deltaT / TauC;

    v = v - v * rmHalf - gC * c * deltaTHalf;
    v = v - v * rmHalf - gC * c * deltaTHalf;
    c = c * (1 - rc + rc * rc / 2);
    
    //The following would be better if we guess an exponentia decay of v dominated by the leakage 
    //v = v - v * rm + v * rm * rm / 2 - gC * c * deltaT;

#endif

    v += thisInputCurrent;

    if (v >= Theta) {               //Is it firing?
      v = H;	        // voltage after spike reset in LIFCA equation
      c += AlphaC;	// Ca concentration update after spike
      Te = t;           // Update of last emission time
      didItFire = true;

      DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
      //oscilloscope on some neurons
      if(glob_n==2480) {
	pStat->LIFCAivu0.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			       currentTime,glob_n,thisInputCurrent,60,c);};
      if(glob_n==135123) {
	pStat->LIFCAivu1.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			       currentTime,glob_n,thisInputCurrent,60,c);};
      if(glob_n==2) {
	pStat->LIFCAivu2.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			       currentTime,glob_n,thisInputCurrent,60,c);};
      if(glob_n==lnp_par.locN-1) {
	pStat->LIFCAivu999.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
				 currentTime,glob_n,thisInputCurrent,60,c);};
      DPSNNverboseEnd();   
    }
  
  } else { //Updates the neuron state during the absolute refractory period (ARP)  

#ifdef PerseoDyn

    c *= exp(-deltaT / TauC);

#else

    // EPA - PSP : 27 July 2016
    // We assume that the "true" equation is:
    // v' = -v / Tau - gC * c / Capacity + I / Capacity
    // c' = -c / TauC
    // where Capacity is fixed to 1
    // In this implementation we integrated the ODE 
    // using a second order approximation of the analitic solution (exponential decay) for "c"
    // a second order approximation in deltaT/2 for v

    rc = deltaT / TauC;
    c = c * (1 - rc + rc * rc / 2);

#endif

  }
  
  Tr = t;
  
  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if (glob_n==2424)
    printf("DYN END neu=%d, v=%f, c=%f, I=%f at %u ms\n",
	   glob_n, v, c, thisInputCurrent,thisTimeStep_ms);
  DPSNNverboseEnd();   

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
    pStat->inputCurrent.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms,  glob_n, thisInputCurrent);
  DPSNNverboseEnd();   

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
  //oscilloscope on some neurons
  if(glob_n==2480) {
    pStat->LIFCAivu0.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			   currentTime,glob_n,thisInputCurrent,v,c);};
  if(glob_n==135123) {
    pStat->LIFCAivu1.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			   currentTime,glob_n,thisInputCurrent,v,c);};
  if(glob_n==2) {
    pStat->LIFCAivu2.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			   currentTime,glob_n,thisInputCurrent,v,c);};
  if(glob_n==lnp_par.locN-1) {
    pStat->LIFCAivu999.write(thisTimeStep_ms/1000,lnp_par.moduloSec,thisTimeStep_ms,
			     currentTime,glob_n,thisInputCurrent,v,c);};
  DPSNNverboseEnd();   
  
}
/*
bool neuronClass::didItFire() {
  if (v >= Theta) {
  // did it fire?
    return(true);
  } else {
    return(false);
  };
};

void neuronClass::afterSpikeNeuralDynamic(const int thisTimeStep_ms) {
  v = H;	// voltage after spike reset in LIFCA equation
  //c += AlphaC;	// Ca concentration update after spike

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
  printf("after spike dynamic should switch off the neuron %d\n",glob_n);
  DPSNNverboseEnd();  
 };
*/
double neuronClass::getPoissonLapse()
{
  double temp;
  temp = (double)1.0 - getUniformLapse();
  if (temp < 0)
    temp = 0;
  return (-logf(temp) * invNuExt);
}

double neuronClass::getUniformLapse()
{
  //#define MAX_VAL 1024
  //return (float) getRandom(MAX_VAL) / MAX_VAL;
  double temp;
  temp = pLocalNetRandDev->Random();
  if (temp >= 1.0){
    printf("Error in getUniformLaps: random number >=1\n");
    fflush(stdout);exit(0);
  }
    
  return temp;
}

double neuronClass::getLastEmissionTime(){
  return Te;
}

#else // else on #ifdef LIFCAneuron

void neuronClass::clearInputCurrent() {
  InputCurrent = 0.0;
};
void neuronClass::addInputCurrent(const float inputValue) {
  float oldI;
  oldI = InputCurrent;
  InputCurrent += inputValue;

  DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
  printf("on neuron %d adding %f to old I=%f, new I=%f in H=%d\n",
	 glob_n, inputValue, oldI, InputCurrent,lnp_par.loc_h);
  DPSNNverboseEnd();
};

void neuronClass::regularDynamicTimeStep(
     const uint32_t thisTimeStep_ms, const float *pInputDistrib, uint32_t * pThalInputCounter) {

  // Adding the thalamic input when default_random_thalamicInput_1 selected  
  uint32_t this_CF, loc_CF;
  float thalamIn;
  float thalamicInputCurrent;

  if(lnp_par.thalamicInput==default_random_thalamicInput_1){
    //default value
    thalamicInputCurrent=20.0;
    this_CF = loc_n / lnp_par.neuronsPerCM;
    if(lnp_par.locCFT >= 1)
      loc_CF = this_CF % (uint32_t) lnp_par.locCFT;
    else
      loc_CF = 0;
    thalamIn = thalamicInputCurrent * pInputDistrib[hashId + loc_CF*lnp_par.neuronsPerCM];
    addInputCurrent(thalamIn);
    if (thalamIn != 0.0) {
      (*pThalInputCounter) ++;
      pStat->thalamicInput.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms, glob_n);
    }
  }
  // End of direct thalamic input current addition  

  // local step of Izhikevich equation updating...
  ///... the membrane potential v and auxiliary variable u

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
  printf("DYN START neu=%d, v=%f, u=%f, I=%f at %u ms\n",
	 glob_n, v, u, InputCurrent, thisTimeStep_ms);
  DPSNNverboseEnd();  

  v+=0.5*((0.04*v+5)*v+140-u+InputCurrent); // for numerical stability
  v+=0.5*((0.04*v+5)*v+140-u+InputCurrent); // here the numerical time step is 0.5 ms
  //b parameter of Izhikevic equations set to 0.2 for all neurons

  if(v >= 30.0)
    v = 30.0;

  u+=a*(0.2*v-u);
      
  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
  printf("DYN END neu=%d, v=%f, u=%f, I=%f at %u ms\n",
	 glob_n, v, u, InputCurrent,thisTimeStep_ms);
  DPSNNverboseEnd();   

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
  if(InputCurrent!=0.0)
    pStat->inputCurrent.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms,  glob_n, InputCurrent);
  DPSNNverboseEnd();   

  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);
  //oscilloscope on some neurons
  //float capV; //just for easier reading of printout 
  //if(v>35.0) {capV=35.0;} else {capV=v;}; //v cap added only on printing
  if(glob_n==5123) {
    pStat->ivu0.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms,  glob_n, InputCurrent, v, u);};
  if(glob_n==7662) {
    pStat->ivu1.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms,  glob_n, InputCurrent, v, u);};
  if(glob_n==2) {
    pStat->ivu2.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms,  glob_n, InputCurrent, v, u);};
  if(glob_n==lnp_par.locN-1) {
    pStat->ivu999.write(thisTimeStep_ms/1000, lnp_par.moduloSec, thisTimeStep_ms,  glob_n, InputCurrent, v, u);};
   DPSNNverboseEnd();   
};

bool neuronClass::didItFire() {
  if (v >= 30) {
  // did it fire? 
    return(true);
  } else {
    return(false);
  };
};

void neuronClass::afterSpikeNeuralDynamic(const uint32_t thisTimeStep_ms) {
  v = -65.0;	// voltage reset  (mV) c param in Izhikevic equation
  u += d;	// recovery variable reset
  if(lnp_par.fastDebugDyn == fastDebugDyn_1) {
      //false value to speed up simulation
      u=0.2*v;
  }
  DPSNNverboseStart(false,thisTimeStep_ms,lnp_par.debugPrintEnable_ms);  
  printf("after spike dynamic should switch off the neuron %d\n",glob_n);
  DPSNNverboseEnd();  
 };

#endif // end of else on #ifdef LIFCAneuron
