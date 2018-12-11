// DPSNN_main.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
 
#include "DPSNN_environmentSelection.h"
#ifdef DALonlyEnvironmentSelected
    //#warning "DPSNN_main.c: as expected this MPI main is empty under DALonly environment" 
#else
#ifdef MPIandDALenvironmentSelected
//#warning "compiling the MPI main"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/sysinfo.h>

#include "DPSNN_localNet.h"
#include "localNetProcess.h"
#include "DPSNN_chrono.h"

int MPIglobH;
int MPIprocessRank;
localNetProcess_localState lnp_ls;
#define lnp_ls_par p->local->parameters
DALProcess theDALProcessDeclaredByEachMPIMain,*p;

chronoClass chronoInit, chronoFire;

int main(int argc, char*argv[]){
  int stopFiring;

    //MPI_Status status;
    char hostName[256];
    struct sysinfo si;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIprocessRank);

    MPI_Comm_size(MPI_COMM_WORLD, (int *)&(MPIglobH));
    if(MPIprocessRank == 0) {
      printf("001- MPIprocessRank= %03d, COMM_WORLD of %d processes\n",
	       MPIprocessRank,
	       MPIglobH);
      fflush(stdout);
    };

    p=&theDALProcessDeclaredByEachMPIMain;
    p->local = &lnp_ls;

    //strcpy(hostName, getenv("HOSTNAME"));//name of computer
    gethostname(hostName,sizeof(hostName));

    if(MPIprocessRank==0){
      //freeMemory=system ("free");
      sysinfo (&si);
      printf ("Before simulation the available free memory on %s is %lu kB\n",
	      hostName?hostName:"NULL",si.freeram/1024 );
      }

    chronoInit.clearChrono();
    chronoInit.startChrono();

    localNetProcess_init(p);

    chronoInit.stopChrono();
    if(MPIprocessRank<=3 || MPIprocessRank>=(MPIglobH-2))
    { printf("CHRONO:MPI-MAIN h=%d on node %s - localNetProcess_init = %f sec \n",
	   MPIprocessRank, hostName?hostName:"NULL",
	   chronoInit.getAccumulatedChrono());
    };
    //stopFiring = localNetProcess_fire(p);

    chronoFire.clearChrono();

    do {
      chronoFire.startChrono();
         stopFiring = localNetProcess_fire(p);
      chronoFire.stopChrono();
    } while (stopFiring == 0);

    localNetProcess_finish(p);
 
    if(MPIprocessRank<2 || MPIprocessRank>=(MPIglobH-2)){  
      printf(
      "CHRONO:MPI-MAIN h=%d on node %s - mean localNetProcess_fire = %f sec on %d samples total fire time = %f\n",
      MPIprocessRank, hostName?hostName:"NULL",
           (double)(  chronoFire.getAccumulatedChrono()/
	             (double) chronoFire.getNumLaps() ),
       chronoFire.getNumLaps(),chronoFire.getAccumulatedChrono());
    };

    MPI_Finalize();
    if(MPIprocessRank == 0) {   
    printf("999- main() h=%03d after MPI_finalize()\n", 
	      MPIprocessRank);
    fflush(stdout);
    }
}

#endif //MPIandDALenvironmentSelected
#endif //selection between DALonly and MPIandDAL environments
