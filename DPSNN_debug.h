// DPSNN_debug.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#ifndef DPSNN_debugIncluded
#define DPSNN_debugIncluded
 
//if the following define is active
//debug regions 

#define DPSNNinsertDebugRegions 1

//examples of use
//DPSNNverboseStart(false,0,lnp_par.debugPrintEnable_ms);
//DPSNNverboseStart(true,0,lnp_par.debugPrintEnable_ms);
//DPSNNverboseStart(false,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
//DPSNNverboseStart(true,thisSimTimeStep_ms,lnp_par.debugPrintEnable_ms);
//DPSNNverboseEnd();

#ifdef DPSNNinsertDebugRegions
  #define DPSNNverboseStart(cond,value1,value2) \
    if(cond)\
      if(value1>=value2) \
	do {
#else
  #define DPSNNverboseStart(cond,value1,value2) \
    if(false) \
      do {
#endif

#ifdef DPSNNinsertDebugRegions
  #define DPSNNverboseEnd() \
    fflush(stdout);} while (0)
#else
  #define DPSNNverboseEnd() \
      } while (0)
#endif

#endif //multiple inclusion guard
