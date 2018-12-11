// DPSNN_sort.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#include "DPSNN_debug.h"
#include "DPSNN_dataStructDims.h"
#include "DPSNN_localNet.h"

void localNetClass::sortTargetHostInForwardSynList(
     synapseClass *synList, synapseClass *q1,
     uint32_t synCount)
{
  //if globH==1 no sort needed
  if(lnp_par.globH > 1) {//sort needed
    DPSNNverboseStart(false,1,0);  
      printf("h=%03d, BEFORE radixSort START on %d synapses\n",
	     lnp_par.loc_h, synCount);
      fflush(stdout);
    DPSNNverboseEnd();

    #define LSDradixSort
    #ifdef LSDradixSort
    radixSortBinForwardSynList(synList, q1, synCount);    
    #else
      #define bubbleSort
      #ifdef bubbleSort 
      //BEGIN sort of synapses according to the target host
      //this version is absolutely not optimized 
      forwardSynClass swapLocation;
      uint32_t i,j;
      for (i=0;i<synCount;i++) {
        for(j=0;j<(synCount-1-i);j++) {
	  if ((synList[j].post_glob_n   / lnp_par.locN) >
	      (synList[j+1].post_glob_n / lnp_par.locN)) {
	    swapLocation=synList[j];
	    synList[j]=synList[j+1];
	    synList[j+1]=swapLocation;
	  }
        }
      }
      #endif
      #undef bubbleSort
    #endif
  };

}; 
 
void localNetClass::radixSortBinForwardSynList(synapseClass *a, synapseClass *q1, uint32_t n)
{
  //the sort key is post_loc_h
  //NOTE: initially, it sorts from key max to key min
  //then it reverts the order
  uint32_t i;
  //synapseClass q0[DSD__maxForwardLocSyn];
  //forwardSynClass q1[DSD__maxForwardLocSyn];
  uint32_t c0, c1; 
  uint32_t m; 
  uint32_t exp = 1;

  DPSNNverboseStart(false,1,0);  
      printf("h=%03d, radixSort START\n",
       lnp_par.loc_h);
  DPSNNverboseEnd();

  m = a[0].post_glob_n / lnp_par.locN;

  for (i = 0; i < n; i++)
  {
    if ((a[i].post_glob_n/lnp_par.locN) > m)
      m = a[i].post_glob_n/lnp_par.locN;
  }

  DPSNNverboseStart(false,1,0);  
      printf("h=%03d radixSort max post_glob_n/lnp_par.locN detected=%d\n",
	     lnp_par.loc_h, m);
  DPSNNverboseEnd(); 

  while (m / exp > 0)
  {
    c0=0;
    c1=0;
    for (i = 0; i < n; i++) {
      if((((a[i].post_glob_n/lnp_par.locN) / exp) % 2) == 0) {
	a[c0++]=a[i];
      }else{
	q1[c1++]=a[i];
      }
    };

    if(c0+c1!=n) {
      printf("ERROR c0+c1 in radix synSort\n");
      fflush(stdout);
      exit(0);
    };

    for (i = c0; i < n; i++)
      a[i]=q1[i-c0];

    exp *= 2;
  };
};

/* basic algo inserted as doc 
#define MAX 8
#define SHOWPASS

struct r {int key; float value;};

void print(struct r *a, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    printf("{%d %03.0f}\t", a[i].key, a[i].value);
  }
}

void radixSortBinStruct(r *a, int n)
{
  int i;
  r q0[MAX], q1[MAX], b[MAX];
  int c0, c1; 
  int m = a[0].key; 
  int exp = 1;

  for (i = 0; i < n; i++)
  {
    if (a[i].key > m)
      m = a[i].key;
  }
 
  while (m / exp > 0)
  {
    c0=0;
    c1=0;
    for (i = 0; i < n; i++) {
      if((a[i].key / exp) % 2) {
	q0[c0++]=a[i];
      }else{
	q1[c1++]=a[i];
      }
    }
    if(c0+c1!=n) printf("ERROR c0+c1\n");

    for (i = 0; i < c0; i++)
      a[i]=q0[i];
    for (i = c0; i < n; i++)
      a[i]=q1[i-c0];

    exp *= 2; 
    #ifdef SHOWPASS
      printf("\nPASS   : ");
      print(a, n);
    #endif
  }
 
}

int main()
{
  r arr[MAX];
  int i, n;
 
  printf("Enter total elements (n < %d) : ", MAX);
  scanf("%d", &n);
 
  printf("Enter %d Elements : ", n);

  for (i = 0; i < n; i++)
    scanf("%d", &(arr[i].key));

  for (i = 0; i < n; i++)
    arr[i].value = (float) (arr[i].key * 10);
 
  printf("\nARRAY  : ");
  print(arr, n);
 
  radixSortBinStruct(arr, n);
 
  printf("\nSORTED : ");
  print(arr, n);
  printf("\n");
 
  return 0;
}
*/
