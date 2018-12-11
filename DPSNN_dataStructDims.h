// DPSNN_dataStructDims.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include "DPSNN_environmentSelection.h"
#ifdef DALonlyEnvironmentSelected
  #define miniLocN (32768)
  #define miniM (400)
  #define miniSimSpikes (4)
  #define makeInitExcitWeight (6.14)
  #define makeInitInhibWeight (5.12)
  #define makeMaxGlobH (32)
  #define makeActiveLTD
  #define makeActiveLTP
#endif

#define DSD__FOPEN_MAX 128

#ifndef DPSNN_dataStructDimsIncluded
#define DPSNN_dataStructDimsIncluded

//max values used to declare arrays before actual initialization
//surely working: locN=2048, M=512, globH=256, D=22,
#define DSD__maxLocN (miniLocN) // easyconf 16000 fits alternate 2048

#define DSD__maxGlobH (makeMaxGlobH) // easy conf 512 fits alternate 256
//#define DSD__maxGlobH 1
#define DSD__maxGlobN (DSD__maxLocN*DSD__maxGlobH)

#define DSD__maxM (miniM) // easyconf 1000 fits

#ifdef LIFCAneuron
#define DSD__maxD (80) 
#define DSD__delayNum (80)

#define DSD__stencilX_Max makeStencil_Max
#define DSD__stencilY_Max makeStencil_Max
#define DSD__bathEfficacyTemplateX_Max makeStencil_Max
#define DSD__bathEfficacyTemplateY_Max makeStencil_Max
/*
#if (makeOverallTopology==0)
// 0 old style explicit stencils and connectivity
#define DSD__maxModDeltaX (stencilDim-1)
#define DSD__maxModDeltaY (stencilDim-1)
#define DSD__stencilX_Max ((2 * DSD__maxModDeltaX) + 1)
#define DSD__stencilY_Max ((2 * DSD__maxModDeltaY) + 1)
#elif (makeOverallTopology==1)
  //1 homogeneous connectivity

#define DSD__stencilX_Max (makeMaxHomogeneousSide)
#define DSD__stencilY_Max (makeMaxHomogeneousSide)
#else
  #error "undefined makeOverallTopology"
#endif
//#define DSD__bathEfficacyTemplateX_Max (stencilDim)
//#define DSD__bathEfficacyTemplateY_Max (stencilDim)
*/

#define DSD__maxGlobCFT (96*96+1)
#define DSD__maxCMFract (makeMaxCMFract)
#define DSD__subPopTot (27)
#else
#define DSD__maxD (20) //10 easy fit max delay in ms: originally 22
#endif
#define DSD__cSm (10.0) //max synaptic value
//#warning "maxD reduced from 20ms to 2ms for faster dynamic during debug"
//#define DSD__maxD 2 //max delay in ms

//3* quite RISKY this assumption
#define DSD__maxM_pre (2 * DSD__maxM)
#define DSD__maxForwardLocSyn (DSD__maxLocN * DSD__maxM) 
#define DSD__maxBackwardLocSyn (DSD__maxLocN * DSD__maxM)
#define DSD__maxSizeOfSyn 24

/*
 * Don't use the following define because for large configurations (miniLocN>=262144)
 * maxMemoryPool needs an uint64_t value, that can't be correctly set with the #define macro.
 * In the code the maxMemoryPool value is recalculated time by time as follows:
 * maxMemoryPool = (uint64_t)DSD__maxSizeOfSyn * (uint64_t)DSD__maxBackwardLocSyn;
 *
 * #define DSD__maxMemoryPool (DSD__maxSizeOfSyn * DSD__maxBackwardLocSyn)
 */

#define DSD__maxSimultaneousSpikesOnSameTarget (miniSimSpikes)//512 works

#define DSD__circularBufferD (DSD__maxD+2)
//#define DSD__circularBufferD (16) //originally 32: must be a power of 2!!!!!
//#define DSD__circBufferBitwiseMaskD (DSD__circularBufferD - 1) 

//in mean, the number of spiking neurons every ms should be a fraction
#define DSD__maxAxonalSpike (DSD__maxLocN) //EPA: changed 20140204 with PSP
#define DSD__maxSpikesPerMs (DSD__maxLocN) //EPA: changed 20140204 with PSP
#define DSD__msPerFrame (1000)
#define DSD__maxSpikesPerFrame (DSD__maxSpikesPerMs*DSD__msPerFrame)

/*
 * The following section is used to define the size of the synapse weight and timeDerivative (weightType).
 * Please, change accordingly also the conversion factor from float to weightType and viceversa.
 * Note that here are defined also min and max synapse weight allowed.
 */
//#define weightType float
#define weightType int16_t
#define DSD__acceptable_minSynWeigth_f (-6.0)
#define DSD__INT16_MIN -32767
/*
the following section contains a few hint for GPU parallelization
//values used for actual initializations of data structures 
#define DSD__log2_blockDim 2 //change HERE -- for parallelization on GPUs we prepare blocks and grids
#define DSD__log2_gridDim  1  //change HERE -- for parallelization on GPUs we prepare blocks and grids
#define DSD__blockDim (1<<DSD__log2_blockDim)
#define DSD__gridDim (1<<DSD__log2_gridDim)
#define DSD__log2_locN (DSD__log2_blockDim + DSD__log2_gridDim)
*/

#define DSD__randTable (1000)
#define DSD__INT_MAX 2147483647
#define DSD__maxTabNumber 512

#endif
