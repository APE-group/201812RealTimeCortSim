/* Unrolled A2A */
/* int ompi_coll_base_alltoall_intra_pairwise(const void *sbuf, int scount, */
/* 					     struct ompi_datatype_t *sdtype, */
/* 					     void* rbuf, int rcount, */
/* 					     struct ompi_datatype_t *rdtype, */
/* 					     struct ompi_communicator_t *comm, */
/* 					     mca_coll_base_module_t *module) */

/* sbuf = axSpike_forwardBuffer */
/* scount = pktLength */
/* sdtype = MPI_spike */
/* rbuf = axSpike_backwardBuffer */
/* rcount = pktLength */
/* rdtype = MPI_spike */
/* comm = MPI_COMM_WORLD */

/* Perform pairwise exchange - starting from 1 so the local copy is last */
/* it needs rank (lnp_par.loc_h) and size (lnp_par.globH) of communicator */
/* which are global */

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_messagePassing.h"

extern const struct DPSNN_parameters lnp_par;
extern MPI_Aint sizeof_MPI_spike;
extern MPI_Datatype MPI_spike;

int unrolled_A2A(const void *sbuf, int scount, MPI_Datatype sdtype,
		 void *rbuf, int rcount, MPI_Datatype rdtype,
		 MPI_Comm comm) {

  for (uint32_t step = 1; step < lnp_par.globH + 1; step++) {

    /* Determine sender and receiver for this step. */
    int sendto   =                 (lnp_par.loc_h + step) % lnp_par.globH;
    int recvfrom = (lnp_par.globH + lnp_par.loc_h - step) % lnp_par.globH;

    /* Determine sending and receiving locations */
    void * tmpsend = (char*)sbuf + sizeof_MPI_spike * (ptrdiff_t)scount * (ptrdiff_t)sendto;
    void * tmprecv = (char*)rbuf + sizeof_MPI_spike * (ptrdiff_t)rcount * (ptrdiff_t)recvfrom;

    /* send and receive */
    /* this also is unrolled */
    /* err = ompi_coll_base_sendrecv( tmpsend, scount, sdtype, sendto, */
    /* 				   MCA_COLL_BASE_TAG_ALLTOALL, */
    /* 				   tmprecv, rcount, rdtype, recvfrom, */
    /* 				   MCA_COLL_BASE_TAG_ALLTOALL, */
    /* 				   comm, MPI_STATUS_IGNORE, rank); */

    MPI_Request req;

    /* post new irecv */
    MPI_Irecv(tmprecv, lnp_par.pktLength, MPI_spike, recvfrom, MPI_ANY_TAG, MPI_COMM_WORLD, &req);

    /* send data to children */
    MPI_Send(tmpsend, lnp_par.pktLength, MPI_spike, sendto, 0, MPI_COMM_WORLD);

    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }
  /* End of Unrolled A2A */

  return MPI_SUCCESS;
}
