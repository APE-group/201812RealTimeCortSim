// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
#include "DPSNN_environmentSelection.h"

#ifdef DALonlyEnvironmentSelected
  //#warning "dalEmulationUnderMPI.c: as expected this file is empty and not used in DALonly environment" 
#else
    #ifdef MPIandDALenvironmentSelected
 
	#include <stdio.h>
 
	extern int MPIglobH;
	extern int MPIprocessRank;

	int DAL_getIndex(int dimension, DALProcess *p) 
	{
	  return (MPIprocessRank);
	};

	int DAL_send_event(void *message, DALProcess *p) {
	  printf("dummy emulation of DAL_send_event\n");
	  fflush(stdout);
	  return (1);
	};
    #endif //MPIandDALenvironmentSelected
#endif //selection between DALonly and MPIandDAL environments
