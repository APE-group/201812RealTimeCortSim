  // the displacements array for the MPI_Alltoallv
  // is done in DPSNN_localNet_init.c, when making
  // room for the synapses arrays
  MPI_Alltoallv((void *) pForwardAxonalSpikes[0].list,
		(int *)  axSpike_forwardCount,
		(int *)  axSpike_forwardOffset, MPI_spike,
		(void *) pBackwardAxonalSpikes[0].list,
		(int *)  axSpike_backwardCount,
		(int *)  axSpike_backwardOffset, MPI_spike,
		MPI_COMM_WORLD);
  // there should be no need for further copies...

  // ------------------
