#include <stdio.h>
#include "DPSNN_messagePassing.h"

//Could not find the way to pass 'lnp_par' via EXTERN thus passing its
//'loc_h' and 'globH' via parameter copy.
//Type of ?Buf is chosen passing the size:
//4 is int (forwardPrep and backwardPrep)
//16 is axonalSpikeDataOnlyClass (forwardBuffer and backwardBuffer)
//this latter seems to be non-packed...
void save_and_tell(const char* filename,
		   void* fBuf,   void* bBuf,
		    int* fCount,  int* bCount,
		    int* fOffset, int* bOffset,
		   const uint32_t loc_h,
		   const uint32_t globH,
		   const size_t size) {
  FILE * fp;
  fp = fopen(filename, "w");

  fprintf(fp, "forward %p\n",  fBuf);
  fprintf(fp, "backward %p\n", bBuf);

  fprintf(fp, "forwardCount[0]=%u\n",  fCount[0]);
  fprintf(fp, "backwardCount[0]=%u\n", bCount[0]);

  fprintf(fp, "forwardOffset[%u]=%u\n",  loc_h, fOffset[loc_h]);
  fprintf(fp, "backwardOffset[%u]=%u\n", loc_h, bOffset[loc_h]);

  switch (size) {
  case 4: //pointer to uint32_t
    {
      int *pF=(int*)fBuf, *pB=(int*)bBuf;

      for(uint32_t h=0; h < globH; h++) {
	fprintf(fp, "[forward+%u] x %u:", fOffset[h], fCount[h]);
	for (int k=0; k < fCount[h]; k++)
	  fprintf(fp, " %u", pF[fOffset[h]+k]);
	fprintf(fp, "\n");
      }

      for(uint32_t h=0; h < globH; h++) {
	fprintf(fp, "[backward+%u] x %u:", bOffset[h], bCount[h]);
	for (int k=0; k < bCount[h]; k++)
	  fprintf(fp, " %u", pB[bOffset[h]+k]);
	fprintf(fp, "\n");
      }
    }
    break;

  case 16: //pointer to axonalSpikeDataOnlyClass
    {
      axonalSpikeDataOnlyClass *pF=(axonalSpikeDataOnlyClass*)fBuf;
      axonalSpikeDataOnlyClass *pB=(axonalSpikeDataOnlyClass*)bBuf;

      for(uint32_t h=0; h < globH; h++) {
	fprintf(fp, "[fwdBuffer+%u] x %u:", fOffset[h], fCount[h]);
	for (int k=0; k < fCount[h]; k++)
	  fprintf(fp, " {%f,%u}",
		  pF[fOffset[h]+k].originalEmissionTime,
		  pF[fOffset[h]+k].pre_glob_n);
	fprintf(fp, "\n");
      }

      for(uint32_t h=0; h < globH; h++) {
	fprintf(fp, "[bwdBuffer+%u] x %u:", bOffset[h], bCount[h]);
	for (int k=0; k < bCount[h]; k++)
	  fprintf(fp, " {%f,%u}",
		  pB[bOffset[h]+k].originalEmissionTime,
		  pB[bOffset[h]+k].pre_glob_n);
	fprintf(fp, "\n");
      }
    }
    break;

  default: //unknown
    exit(0);
  }

  fclose(fp);
}
