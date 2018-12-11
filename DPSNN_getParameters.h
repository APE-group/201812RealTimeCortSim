#ifndef DPSNN_getParametersIncluded
#define DPSNN_getParametersIncluded

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"

struct DPSNN_parameters getParameters(const uint32_t loc_h_initValue);
void readColumnFile(uint32_t *pSubPopNumber,uint32_t *pNeuronsPerCM);
void readPlasticityParametersFile(DPSNN_parameters *pLnp_par);

#endif
