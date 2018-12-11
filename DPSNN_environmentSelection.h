// DPSNN_environmentSelection.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#ifndef DPSNN_environmentSelectionIncluded
#define DPSNN_environmentSelectionIncluded

//next define must be active when dummy emulating DAL under MPI
#define MPIandDALenvironmentSelected 
//#warning "defined MPIandDALenvironmentSelected" 

//next line when using ETHZ single linux DAL environment
//#define DALonlyEnvironmentSelected
//#warning "defined DALonlyEnvironmentSelected" 

#ifdef MPIandDALenvironmentSelected
  #include <mpi.h>
  #include "dal.h" //different path for dal.h if also MPI used
#else 
  #include <dal.h> //usual path, MPI is not used
#endif 

#endif
