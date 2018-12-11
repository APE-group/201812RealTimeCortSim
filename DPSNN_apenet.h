#ifndef DPSNN_apenet
#define DPSNN_apenet
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <malloc.h>
#include <time.h>
#include "apelib.h"
#include "DPSNN_dataStructDims.h"

#define MAX_NPROCS 256
#define MAX_STREAM_BUFS 6
#define APENET_ALIGNMENT 16

typedef struct{
  int dev; //apenet device (default 0)
  apelink_t link; //apenet device file descriptor
  am_portid_t net_port; //apenet port
  ape_net_address_t net_addr; //net addr (apenet coordinates in the torus)
  struct torus_size {
    int x;
    int y;
    int z;
  } torus_size;  //torus sizes
  int idx; //x_coord + x_size * y_coord + z_size*y_size * z_coord
  int rank; //this process' rank
  int nprocs; //number of processes in this 
  char ** stream_bufs;
  int tx_pkts;
  int rx_pkts;

  

  int pagesize;
  am_portid_t* ports;
  ape_net_address_t* coords;


} APE_RDMA_desc_t;

typedef  struct  __attribute__ ((__packed__)){
  int type;
  int sender_rank;
  am_vaddr_t recv_addr;
} APE_RDMA_hdr_t;


int sbuf_recv_cb(APE_RDMA_desc_t * descriptor, am_vaddr_t sbuf_addr);
//================================================================
// Utilities
//================================================================

int idx_from_coords_and_port(ape_net_address_t coords ,am_portid_t port);
#ifdef APENET_DEBUG
static void
handler(int sig, siginfo_t *si, void *uc)
{
  /* Note: calling printf() from a signal handler is not
     strictly correct, since printf() is not async-signal-safe;
     see signal(7) */
  
  printf("Caught signal %d\n", sig);
  print_siginfo(si);
  signal(sig, SIG_IGN);
}

#endif
//================================================================
// Library functions
//================================================================
//init apelink device and prepare some data structs
int APE_RDMA_Init(APE_RDMA_desc_t * descriptor, int rank, int size);

//clean up and close
int APE_RDMA_Finalize(APE_RDMA_desc_t * descriptor);



#endif
