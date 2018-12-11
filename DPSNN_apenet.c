#ifdef APENET_RDMA
#include <assert.h>
#include <unistd.h>
#include <malloc.h>
#include <string.h>


#include "DPSNN_apenet.h"
#include "DPSNN_debug.h"


int sbuf_recv_cb(APE_RDMA_desc_t * descriptor, am_vaddr_t sbuf_addr)
{
  am_vaddr_t vaddr;
  int retcode =  am_register_buf(descriptor->link, &vaddr, (void*)(sbuf_addr), AM_MAX_MSG_SIZE ,
		  (am_context_t)sbuf_recv_cb, AM_STREAMING_BUF);

  if(retcode)
    {
      printf("[%d] ERROR while registering streaming buf addr=%p vaddr=%p size=%d\n", descriptor->rank, sbuf_addr , vaddr, AM_MAX_MSG_SIZE);
      return 1;
    }
  else
    {
      return 0;
    }
}

int APE_RDMA_Init(APE_RDMA_desc_t * descriptor, int rank, int size)
{
  int i;
  int ret = 0;
  int dev = 0; //default is dev 0
  am_portid_t net_myport;
  ape_net_address_t net_myaddr;
  unsigned int n_proc_x, n_proc_y, n_proc_z;
  int pagesize = sysconf(_SC_PAGESIZE);

  descriptor->tx_pkts = 0;
  descriptor->rx_pkts = 0;

  //printf("[%d] APE Init\n", rank);

  descriptor->link = ape_open(dev, &net_myport);
  
  if(!descriptor->link) {
    printf("can't open APEnet device\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }
  descriptor->dev = dev;
  
  //get torus sizes
  if(APELINK_SUCCESS != ape_get_topology(descriptor->link, &n_proc_x, &n_proc_y, &n_proc_z)) {
    printf("can't get APEnet topology\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }
  net_myaddr = ape_get_self_address(descriptor->link);
  
  descriptor->torus_size.x = n_proc_x;
  descriptor->torus_size.y = n_proc_y;
  descriptor->torus_size.z = n_proc_z;
  
  descriptor->idx = net_myport + (net_myaddr.s.x * AM_NPORTS) + (net_myaddr.s.y * AM_NPORTS * n_proc_x) + (net_myaddr.s.z * AM_NPORTS * n_proc_x * n_proc_y); 
  
  //I'm going to need some mpi data exchange here to fill a/few struct/s with ranks and coordinates
  
  
  //allocate arrays
  descriptor->coords = (ape_net_address_t*)calloc(size, sizeof(ape_net_address_t));
  descriptor->ports =  (am_portid_t*)calloc(size, sizeof(am_portid_t));
#ifdef APENET_DEBUG
 descriptor->stream_bufs = (char**) calloc(MAX_STREAM_BUFS, sizeof(char*)); 
#endif

  //copy this process's info
  (descriptor->coords)[rank].u = net_myaddr.u;
  (descriptor->ports)[rank] = net_myport;

  printf("APE Init: rank=%d coords=%d %d %d port %d\n",rank, net_myaddr.s.x, net_myaddr.s.y, net_myaddr.s.z, net_myport );  
 
  //exchange
  
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
		descriptor->coords, sizeof(ape_net_address_t), MPI_BYTE,
		MPI_COMM_WORLD);
  
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
		descriptor->ports, sizeof(am_portid_t), MPI_BYTE,
		MPI_COMM_WORLD);
  
  //MPI_Barrier(MPI_COMM_WORLD);
  
  //register some streaming for use in EAGER protocol
  
#ifdef APENET_DEBUG
  
  for(i = 0; i<MAX_STREAM_BUFS; i++)
    {
      descriptor->stream_bufs[i] = (char*)memalign(pagesize, AM_MAX_MSG_SIZE);
      am_vaddr_t vaddr = 0;
      
      
      //filled with zeros
      memset(descriptor->stream_bufs[i], 0, AM_MAX_MSG_SIZE);
      
      am_register_buf(descriptor->link, &vaddr, (void*)(descriptor->stream_bufs[i]), AM_MAX_MSG_SIZE ,
  		      (am_context_t)sbuf_recv_cb, AM_STREAMING_BUF);
      printf("[%d] Registered streaming buf addr=%p vaddr=%p size=%d\n", rank,  descriptor->stream_bufs[i], vaddr, AM_MAX_MSG_SIZE);

      descriptor->stream_bufs[i] = (char*)vaddr;
      
    }

#endif
  
  descriptor->rank = rank;
  descriptor->nprocs = size;
  descriptor->pagesize = pagesize;

  //printf("[%d] nprocs=%d, APE Init is done.\n", rank, size);


  return ret;
  
}

//clean up and close
int APE_RDMA_Finalize(APE_RDMA_desc_t * descriptor)
{
  int i;
  int ret = 0;
  int retcode;
  am_vaddr_t vaddr;

  //  DBG("[%d] APE Finalize\n", descriptor->rank);

  //do we need a barrier here?
  MPI_Barrier(MPI_COMM_WORLD);

  assert(descriptor->link);

#ifdef APENET_DEBUG
  for(i = 0; i<MAX_STREAM_BUFS; i++)
    {

      am_unregister_buf(descriptor->link,(am_vaddr_t)(descriptor->stream_bufs[i]) ,
  			(void*)(descriptor->stream_bufs[i]), AM_MAX_MSG_SIZE, AM_STREAMING_BUF);

      free((void*) (descriptor->stream_bufs[i]));
    }

  free((void*)(descriptor->stream_bufs));


  timer_delete(descriptor->timerid);
#endif
  /* Not many things to do, just close dev */
  retcode = ape_close(descriptor->link);
  if(retcode != APELINK_NOERROR)
    printf("[%d] error during ape_close(dev=%d)\n", descriptor->rank, descriptor->dev);



  descriptor->link = NULL;
  free(descriptor->ports);
  free(descriptor->coords);


  //    DBG("[%d] APE Finalize DONE\n", descriptor->rank);
  return ret;
}

#endif
