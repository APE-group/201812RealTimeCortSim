// DPSNN_LIFCAconnectome.cpp
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include "DPSNN_LIFCAconnectome.h"
#include "DPSNN_localNet.h"
#include "DPSNN_neuron.h"
#include "DPSNN_debug.h"
#include <math.h>
#include <stdio.h>

void simpleCM_connectomeClass::initConnectivityParam(
  neuSubPopParamStruct *neuSubPopParam,
  struct DPSNN_parameters *p_lnp_par)
{
  
#define bufferSize (10*DSD__stencilX_Max*DSD__stencilY_Max)
  
  char fileName[30];
  FILE *fp;
  uint32_t x,y,j,p;
  char buffer[bufferSize];
  float initValue[6+DSD__stencilX_Max*DSD__stencilY_Max];
  char *token;
  uint32_t numOfTokenInLine;
  uint32_t sourceSubPop;
  uint32_t targetSubPop;
  uint32_t index;
  uint32_t maxNumOfParam;
  

  switch(p_lnp_par->overallConnectivity) {
    case explicitStencil:
      maxNumOfParam = 6 + (p_lnp_par->stencilX_Max+1)*(p_lnp_par->stencilY_Max+1);
    break;
    case homogeneous:
      maxNumOfParam=6+1;
    break;
    default: {
     printf("ERROR: unkown overallConnectivity in DPSNN_LIFCAconnectome\n");
     fflush(stdout);exit(0); }
  };

  if(p_lnp_par->loc_h==0) 
  printf("DEBUG 1138 LIFCA Connectome p_lnp_par->stencilX_Max = %d\n", p_lnp_par->stencilX_Max);

  for(sourceSubPop=0;sourceSubPop<p_lnp_par->subPopNumber;sourceSubPop++)
    for(targetSubPop=0;targetSubPop<p_lnp_par->subPopNumber;targetSubPop++)
      for(x=0;x<=p_lnp_par->maxModDeltaX;x++)
	for(y=0;y<=p_lnp_par->maxModDeltaY;y++)
	  specConnectivityMatrix[sourceSubPop][targetSubPop][x][y] = 0.0;
  
  //sprintf(fileName,"connectivity_h%d.txt",p_lnp_par->loc_h);
  sprintf(fileName,"connectivity.txt");
  fp = fopen(fileName,"r");
  if (fp == NULL) {
    printf("ERROR opening file connectivity.txt \n");
    fflush(stdout);exit(0);
  }
 
  p=0;
  while(fgets(buffer,bufferSize,fp)!=NULL) {
    initValue[0] = p_lnp_par->subPopNumber;
    initValue[1] = p_lnp_par->subPopNumber;
    for(j=2;j<maxNumOfParam;j++)
      initValue[j] = 0.0;
    numOfTokenInLine = 0;
    token = strtok(buffer," \t\n");
    y=0;
    while (token != NULL) {
      if (token[0] == '#') break;
      numOfTokenInLine++;
      initValue[y] = atof(token);
      y++;
      token = strtok (NULL," \t\n");
    }

    if(numOfTokenInLine>0) {
      if((numOfTokenInLine-6)!=((p_lnp_par->maxModDeltaX+1)*(p_lnp_par->maxModDeltaY+1))){
	printf("ERROR: inconsistency in the number of elements of connectivity.txt file. Please check if the connectivity matrix size matches the stencilDim set in the Makefile.\n");fflush(stdout);exit(0);
    }else{
	p++;
	sourceSubPop = (uint32_t) initValue[0];
	targetSubPop = (uint32_t) initValue[1];
	if((sourceSubPop >= p_lnp_par->subPopNumber) || (targetSubPop >= p_lnp_par->subPopNumber)) {
	  printf("ERROR in connectivity initialization: wrong subPop identifier in connectivity.txt file \n");
	  fflush(stdout);exit(0);
	} else {
	  index=2;
	  neuSubPopParam[sourceSubPop].J[targetSubPop] = initValue[index++];
	  neuSubPopParam[sourceSubPop].DJ[targetSubPop] = initValue[index++];
	  neuSubPopParam[sourceSubPop].DMin[targetSubPop] = initValue[index++];
	  neuSubPopParam[sourceSubPop].DMax[targetSubPop] = initValue[index++];

	  for(y = 0; y <= p_lnp_par->maxModDeltaY; y++)
	    for(x = 0; x <= p_lnp_par->maxModDeltaX; x++)
	      specConnectivityMatrix[sourceSubPop][targetSubPop][x][y] = initValue[index++];
	}  
      }
    }
  }
  fclose(fp);

  DPSNNverboseStart(false,1,0);
  if(p!=(p_lnp_par->subPopNumber*p_lnp_par->subPopNumber)){
    printf("WARNING in connectivity initialization: wrong number of connections described in connectivity.txt file x=%d subPopNumber=%d\n",x,p_lnp_par->subPopNumber);
  }
  DPSNNverboseEnd();

  DPSNNverboseStart(true,1,0);
  if(p_lnp_par->loc_h==0) {
    for(sourceSubPop = 0; sourceSubPop < p_lnp_par->subPopNumber; sourceSubPop++)
      for(targetSubPop = 0; targetSubPop < p_lnp_par->subPopNumber; targetSubPop++) {
	printf("\n specConnectivityMatrix for sourceSubPop=%d  targetSubPop=%d: \n",
	       sourceSubPop,targetSubPop);
	for(y = 0; y <= p_lnp_par->maxModDeltaY; y++) {
	  for(x = 0; x <= p_lnp_par->maxModDeltaX; x++)  
	    printf("%f   ",specConnectivityMatrix[sourceSubPop][targetSubPop][x][y]);
	  printf("\n");
	}

      }
  }
  DPSNNverboseEnd();

#undef bufferSize
}

// Initialization function for the Probability Connection Matrix
// probConnectMatrix[sourceSubPop][X][Y]   <-- this is the stencil

void simpleCM_connectomeClass::initProbConnectMatrix(struct DPSNN_parameters *p_lnp_par) {
  uint32_t x,y;
  uint32_t sourceSubPop;
  uint32_t targetSubPop;

  // The Probability Connection Matrix is initialized to 0
  for(sourceSubPop=0;sourceSubPop<p_lnp_par->subPopNumber;sourceSubPop++)
    for(targetSubPop=0;targetSubPop<p_lnp_par->subPopNumber;targetSubPop++)
      for(x=0;x<DSD__stencilX_Max;x++)
	for(y=0;y<DSD__stencilY_Max;y++)
	  probConnectMatrix[sourceSubPop][targetSubPop][x][y] = 0.0;

  for(sourceSubPop=0;sourceSubPop<p_lnp_par->subPopNumber;sourceSubPop++)
    for(targetSubPop=0;targetSubPop<p_lnp_par->subPopNumber;targetSubPop++)     

    switch(p_lnp_par->overallConnectivity) {
    case explicitStencil: 
       for(x=0;x<=p_lnp_par->maxModDeltaX;x++)
  	 for(y=0;y<=p_lnp_par->maxModDeltaY;y++) {
	    probConnectMatrix[sourceSubPop][targetSubPop][p_lnp_par->maxModDeltaX + x][p_lnp_par->maxModDeltaY + y] = 
	      specConnectivityMatrix[sourceSubPop][targetSubPop][x][y];
	    probConnectMatrix[sourceSubPop][targetSubPop][p_lnp_par->maxModDeltaX - x][p_lnp_par->maxModDeltaY + y] = 
	      specConnectivityMatrix[sourceSubPop][targetSubPop][x][y];
	    probConnectMatrix[sourceSubPop][targetSubPop][p_lnp_par->maxModDeltaX - x][p_lnp_par->maxModDeltaY - y] = 
	      specConnectivityMatrix[sourceSubPop][targetSubPop][x][y];
	    probConnectMatrix[sourceSubPop][targetSubPop][p_lnp_par->maxModDeltaX + x][p_lnp_par->maxModDeltaY - y] = 
	      specConnectivityMatrix[sourceSubPop][targetSubPop][x][y];
	    
	 };
    break;
  case homogeneous:
        for(x=0;x<p_lnp_par->stencilX_Max;x++)
	  for(y=0;y<p_lnp_par->stencilY_Max;y++)
	    probConnectMatrix[sourceSubPop][targetSubPop][x][y] = 
	      specConnectivityMatrix[sourceSubPop][targetSubPop][0][0];
  break;
  default:
    {printf("ERROR: unknown overallTopology in DPSNN_LIFCAconnectome.cpp\n");
    fflush(stdout);
    exit(0);
    }
  }
  DPSNNverboseStart(true,1,0);
  {
    if(p_lnp_par->loc_h==0) {
      for(sourceSubPop = 0; sourceSubPop < p_lnp_par->subPopNumber; sourceSubPop++)
	for(targetSubPop = 0; targetSubPop < p_lnp_par->subPopNumber; targetSubPop++)
	  {
	    printf("\n probConnectivityMatrix for sourceSubPop=%d  targetSubPop=%d: \n",
		   sourceSubPop,targetSubPop);
	    for(y = 0; y < p_lnp_par->stencilY_Max; y++) {
	      for(x = 0; x < p_lnp_par->stencilX_Max; x++) 
		printf("%f   ",probConnectMatrix[sourceSubPop][targetSubPop][x][y]);
	      printf("\n");
	    }
	    printf("\n");
	  }

      fflush(stdout);
    }
  }
  DPSNNverboseEnd();
}

// Initialization of the forward connection matrix for each sourceSubPop:
// connectMatrix[sourceSubPop][targetSubPop][X][Y]
// Each matrix element contains the exact number of synapses
// that must be generated for the couple sourceSubPop-targetSubPop 
// for each specific neighbourhood in fixedNum synapses generation

void simpleCM_connectomeClass::initConnectMatrix(struct DPSNN_parameters *p_lnp_par) {
  uint32_t x,y;
  uint32_t sourceSubPop,targetSubPop;

  for(sourceSubPop = 0; sourceSubPop < p_lnp_par->subPopNumber; sourceSubPop++)
    for(targetSubPop = 0; targetSubPop < p_lnp_par->subPopNumber; targetSubPop++){
      for(x = 0; x < p_lnp_par->stencilX_Max; x++)
	for(y = 0; y < p_lnp_par->stencilY_Max; y++)
	  connectMatrix[sourceSubPop][targetSubPop][x][y] = 
	    (uint32_t)ceil((float)probConnectMatrix[sourceSubPop][targetSubPop][x][y] * 
			   (float)neuSubPopInCM_count[targetSubPop]);
    }

  DPSNNverboseStart(false,1,0);
    printf("\n======================================================== \n\n");
    for(sourceSubPop = 0; sourceSubPop < p_lnp_par->subPopNumber; sourceSubPop++){
      printf(" Connectivity Matrix (number of synapses among subPops) for sourceSubPop %d: \n",
	     sourceSubPop);
      printf("-------------------------------------------------------- \n");
      for(targetSubPop = 0; targetSubPop < p_lnp_par->subPopNumber; targetSubPop++){
	printf("\n sourceSubPop:   %d   ---->    targetSubPop:   %d\n\n",
	       sourceSubPop,targetSubPop);
	for(y = 0; y < p_lnp_par->stencilY_Max; y++) {
	  for(x = 0; x < p_lnp_par->stencilX_Max; x++) 
	    printf("%8d   ",connectMatrix[sourceSubPop][targetSubPop][x][y]);
	  printf("\n");
	}
      }
      printf("\n\n");
    }
    fflush(stdout); 
  DPSNNverboseEnd();
}

void simpleCM_connectomeClass::initBathEfficacyTemplate(struct DPSNN_parameters *p_lnp_par) {
#define bufferSize 100000

  char fileName[30];
  FILE *fp;
  uint32_t x,y;
  char buffer[bufferSize];
  char *token;
  uint32_t numOfTokenInLine;

  for(x=0;x<p_lnp_par->bathEfficacyTemplateX_Max;x++)
    for(y=0;y<p_lnp_par->bathEfficacyTemplateY_Max;y++)
      bathEfficacyTemplate[x][y] = 0;
  
  //sprintf(fileName,"bathEfficacy_h%d.txt",p_lnp_par->loc_h);
  sprintf(fileName,"bathEfficacy.txt");
  fp = fopen(fileName,"r");
  if (fp == NULL) {
    printf("ERROR opening file bathEfficacy.txt \n");
    fflush(stdout);exit(0);
  }

  y=0;
  while(fgets(buffer,bufferSize,fp)!=NULL) {
    numOfTokenInLine = 0;
    token = strtok(buffer," \t\n");
    x=0;
    while (token != NULL){
      if (token[0] == '#') break;
      numOfTokenInLine++;
      bathEfficacyTemplate[x][y] = atoi(token);
      x++;
      token = strtok (NULL," \t\n");
    }
    if(numOfTokenInLine>0){
      if((numOfTokenInLine!=p_lnp_par->bathEfficacyTemplateX_Max) || 
	 (numOfTokenInLine!=p_lnp_par->bathEfficacyTemplateY_Max)){
	printf("ERROR: number of elements in bathEfficacy.txt file doesn't match the stencil size stencilDim set in the Makefile.\n");fflush(stdout);exit(0);
      }else{
      y++;
      }
    }
  }

  fclose(fp);

  DPSNNverboseStart(false,1,0);
  {
      printf("\n Bath Efficacy Matrix: \n");
      for(y = 0; y < p_lnp_par->bathEfficacyTemplateY_Max; y++) {
	for(x = 0; x < p_lnp_par->bathEfficacyTemplateX_Max; x++) 
	  printf("%4d   ",bathEfficacyTemplate[x][y]);
	printf("\n");
      }
      printf("\n");
    fflush(stdout);
  }
  DPSNNverboseEnd();

#undef bufferSize
}

void simpleCM_connectomeClass::initSubPopParamStruct(neuSubPopParamStruct *neuSubPopParam,  struct DPSNN_parameters *p_lnp_par) {
#define maxNumOfParam 12
#define bufferSize 200

  char fileName[30];
  FILE *fp;
  uint32_t x,y,j;
  char buffer[bufferSize];
  char *token;
  uint32_t numOfTokenInLine;
  uint32_t subPop;
  float initValue[maxNumOfParam];

  //Do we need some initialization to default values??
  
  //sprintf(fileName,"column_h%d.txt",p_lnp_par->loc_h);
  sprintf(fileName,"column.txt");
  fp = fopen(fileName,"r");
  if (fp == NULL) {
    printf("ERROR opening file column.txt \n");
    fflush(stdout);exit(0);
  }

  x=0;
  while(fgets(buffer,bufferSize,fp)!=NULL) {
    initValue[0] = p_lnp_par->subPopNumber;
    for(j=1;j<maxNumOfParam;j++)
      initValue[j] = 0.0;
    numOfTokenInLine = 0;
    token = strtok(buffer," \t\n");
    y=0;
    while (token != NULL){
      if (token[0] == '#') break;
      numOfTokenInLine++;
      initValue[y] = atof(token);
      y++;
      token = strtok (NULL," \t\n");
    }
    if(numOfTokenInLine>0){
      x++;
      subPop = (uint32_t) initValue[0];
      if(subPop >= neuSubPopTotal){
	printf("ERROR in subPop initialization: wrong subPop identifier in column.txt file \n");
	fflush(stdout);exit(0);
      } else {
	neuSubPopParam[subPop].count = (uint32_t) initValue[1];
	if(subPop==0)
	  neuSubPopParam[subPop].offset = 0;
	else
	  neuSubPopParam[subPop].offset = neuSubPopParam[subPop-1].offset + neuSubPopParam[subPop-1].count;
	neuSubPopParam[subPop].JExt = initValue[2];
	neuSubPopParam[subPop].DJExt = initValue[3];
	neuSubPopParam[subPop].NuExt = initValue[4];
	neuSubPopParam[subPop].Tau = initValue[5];
	neuSubPopParam[subPop].Theta = initValue[6];
	neuSubPopParam[subPop].H = initValue[7];
	neuSubPopParam[subPop].Tarp = initValue[8];
	neuSubPopParam[subPop].AlphaC = initValue[9];
	neuSubPopParam[subPop].TauC = initValue[10];
	neuSubPopParam[subPop].gC = initValue[11];
      }  
    }
  }
  DPSNNverboseStart(false,1,0);
  if(x!=neuSubPopTotal){
    printf("WARNING in subPop initialization: number of subPops in column.txt file doesn't match the predefined number %d \n",neuSubPopTotal);
  }
  DPSNNverboseEnd();

  fclose(fp);

  DPSNNverboseStart(true,1,0);
  {if(p_lnp_par->loc_h==0) {
      printf("\n SubPopulation parameters: \n");
      for(subPop = 0; subPop < p_lnp_par->subPopNumber; subPop++) 
	printf("%d   %4d   %4d   %g   %g   %7g   %g   %g   %g   %g   %g   %4g   %4g \n",
	       subPop,
	       neuSubPopParam[subPop].count,
	       neuSubPopParam[subPop].offset,
	       neuSubPopParam[subPop].JExt,
	       neuSubPopParam[subPop].DJExt,
	       neuSubPopParam[subPop].NuExt,
	       neuSubPopParam[subPop].Tau,
	       neuSubPopParam[subPop].Theta,
	       neuSubPopParam[subPop].H,
	       neuSubPopParam[subPop].Tarp,
	       neuSubPopParam[subPop].AlphaC,
	       neuSubPopParam[subPop].TauC,
	       neuSubPopParam[subPop].gC);
      fflush(stdout);
    }
  }
  DPSNNverboseEnd();

#undef maxNumOfParam
#undef bufferSize
}

void simpleCM_connectomeClass::describeConnectome(
  struct DPSNN_parameters *p_lnp_par, 
  struct neuSubPopParamStruct *neuSubPopParam) {
  uint32_t iNeuSubPop;
  uint32_t subPop;

  //connectRandDev->SetRandomSeed(uint32_t(time(NULL)));
  //connectRandDev->SetRandomSeed(1);

  //KEY POINT change connectome constants only inside the initialization files 
  
  initSubPopParamStruct(neuSubPopParam, p_lnp_par);
  initBathEfficacyTemplate(p_lnp_par);
 
  //setting offset and counts for neural subpopulations
  for(subPop = 0; subPop < p_lnp_par->subPopNumber; subPop++) {
    neuSubPopInCM_offset[subPop] = neuSubPopParam[subPop].offset;
    neuSubPopInCM_count[subPop] = neuSubPopParam[subPop].count;
  } 
  
  //reporting neural subPop offset and count 
  DPSNNverboseStart(false,1,0);
  for(iNeuSubPop = 0; iNeuSubPop < p_lnp_par->subPopNumber; iNeuSubPop++) {
    printf("\n NEW NEURON SUBPOP: iNeuSubPop=%d \n",
	   iNeuSubPop); fflush(stdout); 
    printf("neuSubPopInCM_count[iNeuSubPop] = %d\n ",
	   neuSubPopInCM_count[iNeuSubPop]); 
    printf("neuSubPopInCM_offset[iNeuSubPop] = %d\n \n",
	   neuSubPopInCM_offset[iNeuSubPop]);
    fflush(stdout);
  };
  DPSNNverboseEnd();

  //defining neural count AT GLOBAL SCALE
  //valid only if there are 3 subPop (Exc-Exc-Inh) for each level !!!
  //TO BE GERALIZED
  {
    uint32_t totalNeuLbExc;
    uint32_t totalNeuLaExc;
    uint32_t totalNeuLaInh;
    totalNeuLbExc = 0;
    totalNeuLaExc = 0;
    totalNeuLaInh = 0;
    for(subPop = 0; subPop < p_lnp_par->subPopNumber; subPop+=3) {
      totalNeuLbExc += neuSubPopInCM_count[subPop];
    }
    for(subPop = 1; subPop < p_lnp_par->subPopNumber; subPop+=3) {
      totalNeuLaExc += neuSubPopInCM_count[subPop];
    }
    for(subPop = 2; subPop < p_lnp_par->subPopNumber; subPop+=3) {
      totalNeuLaInh += neuSubPopInCM_count[subPop];
    }
    p_lnp_par->globNe = (totalNeuLbExc + totalNeuLaExc) * p_lnp_par->globCFT;
    p_lnp_par->globNi = totalNeuLaInh * p_lnp_par->globCFT;
  }
  if((p_lnp_par->globNe + p_lnp_par->globNi) != p_lnp_par->globN) {
    printf(
	   "ERROR simpleCM: globNe =%d +globNi=%d != globN=%d \n",
	   p_lnp_par->globNe, p_lnp_par->globNi, p_lnp_par->globN);
    fflush(stdout);exit(0);
  }; 
   
  //  initSpecConnectivityMatrix();
  initConnectivityParam(neuSubPopParam, p_lnp_par);
  initProbConnectMatrix(p_lnp_par);
  initConnectMatrix(p_lnp_par);
  
  //calculate delayMin and delayMax for the whole problem
  {
    uint32_t sourceSubPop, targetSubPop;
    double delayMin = 1e36;
    double delayMax = 0;
    for(sourceSubPop = 0; sourceSubPop < p_lnp_par->subPopNumber; sourceSubPop++)
      for(targetSubPop = 0; targetSubPop < p_lnp_par->subPopNumber; targetSubPop++){
	if (neuSubPopParam[sourceSubPop].DMax[targetSubPop] > delayMax) 
	  delayMax = neuSubPopParam[sourceSubPop].DMax[targetSubPop];
	if (neuSubPopParam[sourceSubPop].DMin[targetSubPop] < delayMin) 
	  delayMin = neuSubPopParam[sourceSubPop].DMin[targetSubPop];
      }
    p_lnp_par->delayMin = delayMin;
    p_lnp_par->delayMax = delayMax;
    DPSNNverboseStart(false,1,0);
    if(p_lnp_par->loc_h == 0)
      printf("The current problem has a delay range from delayMin=%f to delayMax=%f \n",delayMin,delayMax);
    DPSNNverboseEnd();
  }

};

simpleCM_neuCoordinatesStruct 
simpleCM_connectomeClass::convert_loc_n_h_to_neuCMCoordinates(uint32_t loc_n, uint32_t loc_h, struct DPSNN_parameters *p_lnp_par)
{ //compute coordinates for this neuron i in
      //other coordinate systems: global system and cortical modules system
  simpleCM_neuCoordinatesStruct neuCMCoordinates;
      neuCMCoordinates.loc_n = loc_n;
      neuCMCoordinates.loc_h = loc_h;
      neuCMCoordinates.glob_n = 
	neuCMCoordinates.loc_n + p_lnp_par->locN * loc_h;
      neuCMCoordinates.cfy_n = 
        (neuCMCoordinates.glob_n / p_lnp_par->neuronsPerCM) / 
	p_lnp_par->globCFX;
      neuCMCoordinates.cfx_n = 
        (neuCMCoordinates.glob_n / p_lnp_par->neuronsPerCM) % 
	p_lnp_par->globCFX;
      neuCMCoordinates.inCM_n = 
	neuCMCoordinates.glob_n % p_lnp_par->neuronsPerCM;
      neuCMCoordinates.subPop = 
	getNeuralSubPop(neuCMCoordinates.inCM_n, p_lnp_par);
      neuCMCoordinates.neuralKind = 
        getNeuralKind(neuCMCoordinates.inCM_n, p_lnp_par);
      return(neuCMCoordinates);
};

uint32_t simpleCM_connectomeClass::generateTargetNeuList(
  const simpleCM_neuCoordinatesStruct sourceNeu, 
  simpleCM_targetNeuListStruct *targetNeuList,
  struct DPSNN_parameters *p_lnp_par)
{
  uint32_t x,y;
  uint32_t targetSubPop;
  uint32_t totSynNum;
  simpleCM_neuCoordinatesStruct targetNeu;
  int32_t naive_icfx_target,naive_icfy_target;
  uint32_t lastNeuInSubPop;
  uint32_t avilableTargetNeurons;
  
#ifdef fixedNumSyn
   uint32_t s,k,SynNum;
   int32_t npre;
   double r, dx, offset;
   double SynExtraction[2000]; // Support array in which are collected the poissonian
                                        // extractions which set the post-synaptic neurons with
                                        // a synaptic contact.
   uint32_t exists;
   uint32_t glob_r; 
#endif

  totSynNum = 0;
  // Enable the following line to ensure code reproducibility over more than 1 process
  //connectRandDev->SetRandomSeed(sourceNeu.glob_n);
  //connectRandDev->SetRandomSeed(uint32_t(time(NULL)));  

    //use srand initialization when in the following is used getRandom()
    //srand((glob_n+1)*(thisTimeStep_ms+1));  
    //use neuRandDev->SetRandomSeed initialization when in the following is used neuRandDev->Random()
    // Enable the following line to ensure code reproducibility over more than 1 process
    //connectRandDev->SetRandomSeed(sourceNeu.glob_n);

  //i = 0;


  for(targetSubPop=0; targetSubPop<p_lnp_par->subPopNumber; targetSubPop++)
    for(x = 0; x < p_lnp_par->stencilX_Max; x++)
      for(y = 0; y < p_lnp_par->stencilY_Max; y++){
	//deltaCM.cfx = x - 3;
	//deltaCM.cfy = (y - 3) * (-1);

	targetNeu.cfx_n = p_lnp_par->globCFX;  //Initialized to an absourd value
	targetNeu.cfy_n = p_lnp_par->globCFY;  //Initialized to an absourd value

	//naive_icfx_target = deltaCM.cfx + sourceNeu.cfx_n;
	//naive_icfy_target = deltaCM.cfy + sourceNeu.cfy_n;
	if(p_lnp_par->overallConnectivity==explicitStencil) {
	//0 old style explicit stencils and connectivity
	  naive_icfx_target = sourceNeu.cfx_n + x - p_lnp_par->maxModDeltaX;
	  naive_icfy_target = sourceNeu.cfy_n + y - p_lnp_par->maxModDeltaY;
	} else {
	  if (p_lnp_par->overallConnectivity==homogeneous) {
	    //1 homogeneous connectivity
            naive_icfx_target=x;
            naive_icfy_target=y;
	  } else {
	    printf("ERROR in LIFCAconnectome generateTargetNeuList, unknown value of overallConnectivity\n");
	    fflush(stdout);
	    exit(0);
	  }
	}
    
	//naive, because there should be exceptions at the boundaries
	if((naive_icfx_target >= 0) &&
	   (naive_icfx_target < (int32_t)p_lnp_par->globCFX))
	  targetNeu.cfx_n = (uint32_t)naive_icfx_target;

        #ifdef periodicBoundaryConditions
	else 
        #warning " periodic boundary conditions on CFX"
	  targetNeu.cfx_n = (uint32_t)(naive_icfx_target % p_lnp_par->globCFX);
        #endif

	//naive, because there should be exceptions at the boundaries
	if((naive_icfy_target >= 0) &&
	   (naive_icfy_target < (int32_t)p_lnp_par->globCFY))
	  targetNeu.cfy_n = (uint32_t) naive_icfy_target;
        #ifdef periodicBoundaryConditions
	else 
          #warning " periodic boundary conditions on CFY"
	  targetNeu.cfy_n = (uint32_t)(naive_icfy_target % p_lnp_par->globCFY);
        #endif

	if(((targetNeu.cfx_n) < p_lnp_par->globCFX) &&
	   ((targetNeu.cfy_n) < p_lnp_par->globCFY)){
	  if((targetNeu.cfx_n<0) || (targetNeu.cfx_n>=p_lnp_par->globCFT) ||
	     (targetNeu.cfy_n<0) || (targetNeu.cfy_n>=p_lnp_par->globCFT)){
	    printf("ERROR in generateTargetNeuList: targetNeu coord (cfx_n,cfy_n) not valid");
	    fflush(stdout);exit(1);
	  }

#ifdef fixedNumSyn
	  
	  SynNum = connectMatrix[sourceNeu.subPop][targetSubPop][x][y];
	  if (SynNum > 0) {      

      /* The following is the code for fixed syn generation */
      /* that uses directly a random number generation      */
      /* for a NOT faster version.                              */
	    
	      /* START OF NOT FAST CODE FOR FIXED-SYN*/
	    /*
	      for (s=0; s<SynNum; s++) {
		do { 
		  exists = 0;// to avoid multiple synapses

		  r = neuSubPopInCM_offset[targetSubPop] + 
		      (uint32_t(pLocalNetRandDev->Random()*neuSubPopInCM_count[targetSubPop]))
		      %neuSubPopInCM_count[targetSubPop];
		  glob_r = conv_neuIdInCM_to_glob_n(r,targetNeu.cfx_n,targetNeu.cfy_n);
		  for (k=0;k<totSynNum;k++) {
		    if (targetNeuList[k].glob_n==glob_r) {
		      //no duplicated synapses with same target allowed
		      exists = 1; //already existing found
		      DPSNNverboseStart(false,0,1);
		      printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
			     p_lnp_par->loc_h, sourceNeu.glob_n, r);
		      DPSNNverboseEnd();
		    }; // synapse already exists
		  }; 
		} while (exists == 1); // if exist==1 try using another random r

		targetNeu.inCM_n = r;
		targetNeu.glob_n = glob_r;

		//targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
		targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
		targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
		targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
		
		targetNeuList[totSynNum].glob_n = targetNeu.glob_n;
		targetNeuList[totSynNum].loc_h = targetNeu.loc_h;
		targetNeuList[totSynNum].subPop = targetNeu.subPop;
		targetNeuList[totSynNum].neuralKind = targetNeu.neuralKind;
	      
		totSynNum++;
	      }
	    */
	      /* END OF NOT FAST CODE FOR FIXED-SYN*/

      /* The following is the code for fixed syn generation as in Perseo */
      /* Please, comment the previouse one and use the following         */
      /* to reproduce the same behaviour as in Perseo.                   */

	      /* START OF PERSEO-LIKE CODE FOR FIXED-SYN*/
	      
	      for (s=0, r=0.0; s<SynNum; s++) {
		  r += pLocalNetRandDev->ExpDev();
		  SynExtraction[s] = r;
		}
	      dx = SynExtraction[SynNum-1] / (neuSubPopInCM_count[targetSubPop] - SynNum);
	      offset = -(pLocalNetRandDev->Random()) * (SynExtraction[SynNum-1] + SynNum * dx);
	      for (s=0, k=0; s<SynNum; s++) {
		npre = (int)floor(s+(offset+SynExtraction[s])/dx);
		if (npre >= 0) {
		  targetNeu.inCM_n = npre + neuSubPopInCM_offset[targetSubPop];
		  targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n, 
				     targetNeu.cfx_n, targetNeu.cfy_n, p_lnp_par);
		  //targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
		  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
		  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n, p_lnp_par);
		  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n, p_lnp_par);
		  
		  targetNeuList[totSynNum].glob_n = targetNeu.glob_n;
		  targetNeuList[totSynNum].loc_h = targetNeu.loc_h;
		  targetNeuList[totSynNum].subPop = targetNeu.subPop;
		  targetNeuList[totSynNum].neuralKind = targetNeu.neuralKind;

		  k++;
		  totSynNum++;
		}
	      }
	      for (s=0; k<SynNum; k++, s++) {
		npre = neuSubPopInCM_count[targetSubPop] + (int)floor(s+(offset+SynExtraction[s])/dx);
		targetNeu.inCM_n = npre + neuSubPopInCM_offset[targetSubPop];
		targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n, 
							    targetNeu.cfx_n, targetNeu.cfy_n);
		//targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
		targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
		targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n);
		targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n);
		
		targetNeuList[totSynNum].glob_n = targetNeu.glob_n;
		targetNeuList[totSynNum].loc_h = targetNeu.loc_h;
		targetNeuList[totSynNum].subPop = targetNeu.subPop;
		targetNeuList[totSynNum].neuralKind = targetNeu.neuralKind;
	      
		totSynNum++;
	      }
	    
	      /* END OF PERSEO-LIKE CODE FOR FIXED-SYN*/
	  }
	      

#else
	  if(probConnectMatrix[sourceNeu.subPop][targetSubPop][x][y]>0){
	    if(probConnectMatrix[sourceNeu.subPop][targetSubPop][x][y]>0.00001){
	      targetNeu.inCM_n = neuSubPopInCM_offset[targetSubPop] - 1;
	      lastNeuInSubPop = neuSubPopInCM_offset[targetSubPop] + neuSubPopInCM_count[targetSubPop];
	      avilableTargetNeurons = lastNeuInSubPop - targetNeu.inCM_n;
	      while((targetNeu.inCM_n += getEmptyRandomSynapses(sourceNeu.subPop,targetSubPop,x,y,
					 avilableTargetNeurons)) < lastNeuInSubPop){
		targetNeu.glob_n = conv_neuIdInCM_to_glob_n(targetNeu.inCM_n, targetNeu.cfx_n, targetNeu.cfy_n, p_lnp_par);
		//targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
		targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
		targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n, p_lnp_par);
		targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n, p_lnp_par);

		targetNeuList[totSynNum].glob_n = targetNeu.glob_n;
		targetNeuList[totSynNum].loc_h = targetNeu.loc_h;
		targetNeuList[totSynNum].subPop = targetNeu.subPop;
		targetNeuList[totSynNum].neuralKind = targetNeu.neuralKind;

		totSynNum++;
	      }
	    } else {
	      // At the moment this branch of if-statement is not used because desn't produce 
	      // correct results. The calculation of SynNum in the following must be adjusted 
	      uint32_t s,k,SynNum;
	      double r;
	      uint32_t exists;
	      uint32_t glob_r;
	      uint32_t percent;
	      uint32_t randVar;
	      float connectProb;
	      //float percentF;
	      //float randVarF;
	     
	      connectProb = probConnectMatrix[sourceNeu.subPop][targetSubPop][x][y];
	      SynNum = connectMatrix[sourceNeu.subPop][targetSubPop][x][y];
	      
	      // With the following, too many syn generated
	      percent = (uint32_t)(SynNum * 0.05);
	      randVar = (uint32_t)(pLocalNetRandDev->Random() * (float)(2 * percent)) - percent;
	      SynNum += randVar;

	      // With the following, too few syn generated
	      //percentF = connectProb * 0.05;
	      //randVarF = (pLocalNetRandDev->Random() * (2 * percentF)) - percentF;
	      //SynNum = (uint32_t)((connectProb + randVar) * (float)neuSubPopInCM_count[targetSubPop]);

	      // something better then the previous two solution must be found!!!

	      if(SynNum>0){
		uint32_t synList[SynNum];
		r = neuSubPopInCM_offset[targetSubPop] + 
		  (uint32_t(pLocalNetRandDev->Random()*neuSubPopInCM_count[targetSubPop]))
		  %neuSubPopInCM_count[targetSubPop];
		glob_r = conv_neuIdInCM_to_glob_n(r,targetNeu.cfx_n,targetNeu.cfy_n, p_lnp_par);
		synList[0] = glob_r;
 		targetNeu.inCM_n = r;
		targetNeu.glob_n = glob_r;

		//targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
		targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
		targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n, p_lnp_par);
		targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n, p_lnp_par);
		
		targetNeuList[totSynNum].glob_n = targetNeu.glob_n;
		targetNeuList[totSynNum].loc_h = targetNeu.loc_h;
		targetNeuList[totSynNum].subPop = targetNeu.subPop;
		targetNeuList[totSynNum].neuralKind = targetNeu.neuralKind;
	      
		totSynNum++;

		for (s=1; s<SynNum; s++) {
		  do { 
		    exists = 0;// to avoid multiple synapses
		    r = neuSubPopInCM_offset[targetSubPop] + 
		      (uint32_t(pLocalNetRandDev->Random()*neuSubPopInCM_count[targetSubPop]))
		      %neuSubPopInCM_count[targetSubPop];
		    glob_r = conv_neuIdInCM_to_glob_n(r,targetNeu.cfx_n,targetNeu.cfy_n, p_lnp_par);
		    for (k=0;k<s;k++) {
		      if (synList[k]==glob_r) {
			//no duplicated synapses with same target allowed
			exists = 1; //already existing found
			DPSNNverboseStart(false,1,0);
			printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
			       p_lnp_par->loc_h, sourceNeu.glob_n, glob_r);
			DPSNNverboseEnd();
		      } // synapse already exists
		    }
		  } while (exists == 1); // if exist==1 try using another random r
	      
		  synList[s] = glob_r;
		  targetNeu.inCM_n = r;
		  targetNeu.glob_n = glob_r;

		  //targetNeu.loc_n = targetNeu.glob_n % p_lnp_par->locN;
		  targetNeu.loc_h = targetNeu.glob_n / p_lnp_par->locN;
		  targetNeu.subPop = getNeuralSubPop(targetNeu.inCM_n, p_lnp_par);
		  targetNeu.neuralKind = getNeuralKind(targetNeu.inCM_n, p_lnp_par);
		
		  targetNeuList[totSynNum].glob_n = targetNeu.glob_n;
		  targetNeuList[totSynNum].loc_h = targetNeu.loc_h;
		  targetNeuList[totSynNum].subPop = targetNeu.subPop;
		  targetNeuList[totSynNum].neuralKind = targetNeu.neuralKind;
	      
		  totSynNum++;
		}
	      }
	    }
	  }
#endif
	} //end on if on being inside CF
      } //end for on target, on x, on y
  return(totSynNum);
}

uint32_t simpleCM_connectomeClass::generateTargetNeu(
    uint32_t mId,
    simpleCM_neuCoordinatesStruct sourceNeu,
    uint32_t subPopId,
    struct DPSNN_parameters *p_lnp_par) 
{
  uint32_t k, r, exists;

  //generation of target neural id in target CM
  do { 
    //search for a valid post-synaptic neuron index r
    if(countRandSynGen++ > (p_lnp_par->M*8)) {
      printf("ERROR: in gen.For.Conn. more than M*8 random attempts on neu=%d syn=%d \n",sourceNeu.glob_n,mId); 
      fflush(stdout);exit(0);
    };
    exists = 0;// to avoid multiple synapses

    //r = neuSubPopInCM_offset[subPopId] + (sourceNeu.inCM_n + 
    //getRandom_r(&seedForSynapses,neuSubPopInCM_count[subPopId]))%neuSubPopInCM_count[subPopId];
    r = neuSubPopInCM_offset[subPopId] + (sourceNeu.inCM_n + 
	  uint32_t(pLocalNetRandDev->Random()*neuSubPopInCM_count[subPopId]))%neuSubPopInCM_count[subPopId];

    if (sourceNeu.glob_n==r) { // no self-synapses allowed
      exists=1;
      DPSNNverboseStart(false,0,1);
      printf("gen.Forw.Conn. R-h%d REJECTED source glob_n %d == target %d\n",
	     p_lnp_par->loc_h, sourceNeu.glob_n, r);
      DPSNNverboseEnd();
    }; 
    for (k=0;k<mId;k++) {
      if (synListOfThisNeu[k]==r) {
	//no duplicated synapses with same target allowed
	exists = 1; //already existing found
	DPSNNverboseStart(false,0,1);
	printf("gen.Forw.Conn. S-h%d already existing %d->%d REJECTED \n",
	       p_lnp_par->loc_h, sourceNeu.glob_n, r);
	DPSNNverboseEnd();
      }; // synapse already exists
    };  
  } while (exists == 1); // if exist==1 try using another random r

  synListOfThisNeu[mId] = r;

  return(r);
};

float simpleCM_connectomeClass::generateSourceNeuCExt(uint32_t cfx_val, uint32_t cfy_val, struct DPSNN_parameters *p_lnp_par) 
{
  uint32_t cfx = cfx_val;
  uint32_t cfy = cfy_val;

  if(cfx >= p_lnp_par->globCFX/2)
    cfx = -1 * (cfx - p_lnp_par->globCFX + 1);
  if(cfy >= p_lnp_par->globCFY/2)
    cfy = -1 * (cfy - p_lnp_par->globCFY + 1);

  if(cfx >= p_lnp_par->bathEfficacyTemplateX_Max)
    cfx = p_lnp_par->bathEfficacyTemplateX_Max - 1;
  if(cfy >= p_lnp_par->bathEfficacyTemplateY_Max)
    cfy = p_lnp_par->bathEfficacyTemplateY_Max - 1;

  return bathEfficacyTemplate[cfx][cfy];
}

neuSubPopEnum simpleCM_connectomeClass::getNeuralSubPop(uint32_t iNeuInCM, struct DPSNN_parameters *p_lnp_par) 
{
  uint32_t subPop;
  //  if(neuSubPopInCM_count[RS] != 0) // Check not needed now. To be reintroduced when more subPop available	  
  for(subPop = 0; subPop < p_lnp_par->subPopNumber; subPop++) {
    if((iNeuInCM >= neuSubPopInCM_offset[subPop])  && 
       (iNeuInCM < (neuSubPopInCM_count[subPop] + 
		    neuSubPopInCM_offset[subPop]) ) )
      return(getSubPopEnum(subPop));
  }
  /*
  if((iNeuInCM >= neuSubPopInCM_offset[LbExc])  && 
     (iNeuInCM < (neuSubPopInCM_count[LbExc] + 
		  neuSubPopInCM_offset[LbExc]) ) )
    return(LbExc);
  else if((iNeuInCM >= neuSubPopInCM_offset[LaExc])  && 
     (iNeuInCM < (neuSubPopInCM_count[LaExc] + 
		  neuSubPopInCM_offset[LaExc]) ) )
    return(LaExc);
  else
    return(LaInh);
  */

  printf("ERROR: in gettin neuron subpopulation - subPopNumber=%d iNeuInCM=%d\n",p_lnp_par->subPopNumber,iNeuInCM);
  fflush(stdout);exit(0);
  return(neuSubPopTotal);
}

neuralKindEnum simpleCM_connectomeClass::getNeuralKind(uint32_t iNeuInCM, struct DPSNN_parameters *p_lnp_par) 
{
  neuSubPopEnum subPopInCM;
  neuralKindEnum neuralKind;
  subPopInCM = getNeuralSubPop(iNeuInCM,p_lnp_par);
  switch(subPopInCM) {
  //case LbExc:
  case L11:
  case L21:
  case L31:
  case L41:
  case L51:
  case L61:
  case L71:
  case L81:
  case L91:
    neuralKind = excitatoryLbExc;
    break;
  //case LaExc:
  case L12:
  case L22:
  case L32:
  case L42:
  case L52:
  case L62:
  case L72:
  case L82:
  case L92:
    neuralKind = excitatoryLaExc;
    break;
  //case LaInh:
  case L13:
  case L23:
  case L33:
  case L43:
  case L53:
  case L63:
  case L73:
  case L83:
  case L93:
    neuralKind = inhibitoryLaInh;
    break;
  default:
    printf("ERROR: unrecogn. subPop->neural kind\n");
    fflush(stdout);exit(0);
    break;
  };
return(neuralKind);
};

uint32_t simpleCM_connectomeClass::conv_neuIdInCM_to_glob_n(
  uint32_t iNeuIdInCM, uint32_t icfx, uint32_t icfy, struct DPSNN_parameters *p_lnp_par) {
  uint32_t neuGlobId;
  if((iNeuIdInCM >= p_lnp_par->neuronsPerCM) ||
     (icfx >= p_lnp_par->globCFX)            ||
     (icfy >= p_lnp_par->globCFY))
    {printf(
     "ERROR in conv_neuIdInCM_to_glob_n, iNeuIdInCM=%d,icfx=%d,icfy=%d\n",  
     iNeuIdInCM,icfx,icfy);fflush(stdout);exit(0);};

  neuGlobId = iNeuIdInCM + 
      icfx * p_lnp_par->neuronsPerCM +
      icfy * p_lnp_par->neuronsPerCM *p_lnp_par->globCFX;

  if(neuGlobId>=p_lnp_par->globN) {printf(
   "ERROR conv_neuIdInCM_to_glob_n, would ret. neuGlobId=%d out of range\n",  
     neuGlobId);fflush(stdout);exit(0);};

  return(neuGlobId);
};

void simpleCM_connectomeClass::report() {
};

char* simpleCM_connectomeClass::getSubPopName(enum neuSubPopEnum neuSubPop) 
{
   switch (neuSubPop) 
     {
       //case LbExc: return "LbExc";break;
       //case LaExc: return "LaExc";break;
       //case LaInh: return "LaInh";break;
     case L11: return (char*)"L11";break;
     case L12: return (char*)"L12";break;
     case L13: return (char*)"L13";break;
     case L21: return (char*)"L21";break;
     case L22: return (char*)"L22";break;
     case L23: return (char*)"L23";break;
     case L31: return (char*)"L31";break;
     case L32: return (char*)"L32";break;
     case L33: return (char*)"L33";break;
     case L41: return (char*)"L41";break;
     case L42: return (char*)"L42";break;
     case L43: return (char*)"L43";break;
     case L51: return (char*)"L51";break;
     case L52: return (char*)"L52";break;
     case L53: return (char*)"L53";break;
     case L61: return (char*)"L61";break;
     case L62: return (char*)"L62";break;
     case L63: return (char*)"L63";break;
     case L71: return (char*)"L71";break;
     case L72: return (char*)"L72";break;
     case L73: return (char*)"L73";break;
     case L81: return (char*)"L81";break;
     case L82: return (char*)"L82";break;
     case L83: return (char*)"L83";break;
     case L91: return (char*)"L91";break;
     case L92: return (char*)"L92";break;
     case L93: return (char*)"L93";break;
     default: return (char*)"ERROR: wrong neuSubPop name";break;
   }
}

neuSubPopEnum simpleCM_connectomeClass::getSubPopEnum(uint32_t neuSubPop) 
{
   switch (neuSubPop) 
     {
      //case 0: return LbExc;break;
      //case 1: return LaExc;break;
      //case 2: return LaInh;break;
    case 0: return L11;break;
    case 1: return L12;break;
    case 2: return L13;break;
    case 3: return L21;break;
    case 4: return L22;break;
    case 5: return L23;break;
    case 6: return L31;break;
    case 7: return L32;break;
    case 8: return L33;break;
    case 9: return L41;break;
    case 10: return L42;break;
    case 11: return L43;break;
    case 12: return L51;break;
    case 13: return L52;break;
    case 14: return L53;break;
    case 15: return L61;break;
    case 16: return L62;break;
    case 17: return L63;break;
    case 18: return L71;break;
    case 19: return L72;break;
    case 20: return L73;break;
    case 21: return L81;break;
    case 22: return L82;break;
    case 23: return L83;break;
    case 24: return L91;break;
    case 25: return L92;break;
    case 26: return L93;break;
    default: printf("ERROR in getSubPopEnum\n");fflush(stdout);exit(0);break;
   }
}

/*------------------------*
 *  getEmptySynapses_RAN  *
 *------------------------*/

/**
 *  Return the number plus one of non-existent (empty) synapses
 *  between two consecutive post-synaptic neurons, following 
 *  a Bernoulli distribution with probability p.
 */

//#define INT_MAX       2147483647    /* maximum (signed) int value */

uint32_t simpleCM_connectomeClass::getEmptyRandomSynapses (uint32_t sourceSubPop, uint32_t targetSubPop, 
							   uint32_t X, uint32_t Y,uint32_t maxN)
{
  static uint32_t n;
  static double  r, P, C;
  static double connectivity;

  connectivity = (double)probConnectMatrix[sourceSubPop][targetSubPop][X][Y];

  DPSNNverboseStart(true,1,0);
  if (sourceSubPop < 0 && X < 0 && Y < 0){
    printf("ERROR 1 in getEmptyRandomSynapses\n");
    fflush(stdout);exit(0);
    //return 0;
  }
  DPSNNverboseEnd();

  DPSNNverboseStart(true,1,0);
  if (probConnectMatrix[sourceSubPop][targetSubPop][X][Y] <= 0.0){
    printf("ERROR 2 in getEmptyRandomSynapses\n");
    fflush(stdout);exit(0);
    //n = DSD__INT_MAX;
  }
  DPSNNverboseEnd();

  r = pLocalNetRandDev->Random();
  DPSNNverboseStart(true,1,0);
  if(r == 1.0){
    printf("ERROR 4 in getEmptyRandomSynapses\n");
    fflush(stdout);exit(0);
  }
  DPSNNverboseEnd();
  
  P = connectivity;
  C = P;
  n = 1;
  while ((C < r)&&(n<=maxN)) {
    n++;
    P *= 1 - connectivity;
    C += P;
  }
  DPSNNverboseStart(true,1,0);
  if (n<=0) {
    printf("ERROR 3 in getEmptyRandomSynapses\n");
    fflush(stdout);exit(0);
    //n = DSD__INT_MAX;
  }
  DPSNNverboseEnd();
  
  DPSNNverboseStart(false,1,0);
  printf("getEmptyRandomSynapses: generated new delta from next random syn is %d\n",n);
  DPSNNverboseEnd();

  return n;
}

void simpleCM_connectomeClass::initLocalNetRandDevPointer(randdevClass *pRandDevPoint_initValue) {
  pLocalNetRandDev = pRandDevPoint_initValue; 
};

float simpleCM_connectomeClass::getConnectProb
(uint32_t sourceSubPop,uint32_t targetSubPop,uint32_t stencilX,uint32_t stencilY){
return probConnectMatrix[sourceSubPop][targetSubPop][stencilX][stencilY];
};
