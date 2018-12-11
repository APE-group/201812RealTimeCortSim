#ifndef SAVE_AND_TELL
#define SAVE_AND_TELL

#include <stdio.h>
#include "DPSNN_messagePassing.h"

void save_and_tell(const char* filename,
		   void* fBuf,   void* bBuf,
		    int* fCount,  int* bCount,
		    int* fOffset, int* bOffset,
		   const uint32_t loc_h,
		   const uint32_t globH,
		   const size_t size);
#endif //inclusion guard end
