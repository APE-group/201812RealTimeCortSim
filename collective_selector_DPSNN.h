// DPSNN_messagePassing.h
// DPSNN_*.* Distribution/Parallelization 
// performed by Pier Stanislao Paolucci (Roma, Italy project start date 2011),
// starting from the sequential code
// SPNET.cpp: Spiking neural network with axonal conduction delays and STDP
// Created by Eugene M. Izhikevich, May 17, 2004, San Diego, CA

// reference paper, named [IzhPol] in the following 
// Eugene M. Izhikevich "Polychronization: Computation with Spikes" 
// Neural Computation 18, 245-282 (2006)

#ifndef DPSNN_collectiveSelector
#define DPSNN_collectiveSelector

int unrolled_A2A(const void *, int, MPI_Datatype,
		 void *, int, MPI_Datatype,
		 MPI_Comm);

#endif //multiple inclusion guard
