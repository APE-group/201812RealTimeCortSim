// DPSNN_parameters.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#ifndef DPSNN_parametersIncluded
#define DPSNN_parametersIncluded
#include <stdint.h>
#include "DPSNN_dataStructDims.h"

#define ANALOG_DEPTH (256)
 
enum thalamicInputEnum { 
      no_thalamicInput_0, 
      default_random_thalamicInput_1, 
      linear_glob_n_vs_time_thalamicInput_2,
      linear_neuronsPerCM_vs_time_thalamicInput_3,
      random_perCM_perMs_thalamicInput_4,
      pseudo_random_neuronsPerCM_vs_time_thalamicInput_5,
      randTable_thalamicInput_6};

enum synGenEnum { 
      nonRandom_synGen_0, 
      default_random_synGen_1,
      simpleCorticalModule_synGen_2,
      randTable_simpleCorticalModule_synGen_3,
      LIFCACorticalModule_synGen_4,
      undefinedSynGen
};

enum fastDebugDynEnum { 
      default_nonActiveFastDebugDyn_0, 
      fastDebugDyn_1, };

enum howManyOutputFilesEnum {
      default_shortOutput_0,
      longOutput_1,
};

enum SongSTDP_paramFileEnum {
  antiCausalCoeffEnum=1, STDPmultiplierEnum, tauSTDP_msEnum,
  causalSTDPmaxAbsValueEnum, antiCausalSTDPmaxAbsValueEnum,
  tauDerivativeDecay_sEnum, derivativeDriftEnum,numOfSTDPparamsEnum
};

struct SongSTDP_paramStruct {
  float STDPmultiplier;
  float antiCausalCoeff;
  float tauSTDP_ms;
  float causalSTDPmaxAbsValue;
  float antiCausalSTDPmaxAbsValue;
  float tauDerivativeDecay_s;
  float derivativeDrift;
  
};
  
enum plasticityAlgoEnum {
  SongSTDPenum=1, numOfPlasticityModelsEnum
};

enum overallConnectivityEnum {
  explicitStencil=0,homogeneous=1,undefinedConnectivity
};

enum collectivesAlgoEnum {
  standardMPI=0, unrolledMPI=1, undefinedCollectivesAlgo
};

struct DPSNN_parameters {
  uint32_t globH; //total number of processes
  uint32_t loc_h; //local process if
  uint32_t locN;  //number of neurons managed by each process
  uint32_t M;
  uint32_t D;
  uint32_t globN;//total number of neurons
  uint32_t globNe;
  uint32_t globNi;
  uint32_t totalSimTime_ms;
  uint32_t thalamicInput;
  uint32_t thalamicInputFreq;
  synGenEnum synGen;
  fastDebugDynEnum fastDebugDyn;
  howManyOutputFilesEnum howManyOutputFiles;
  uint32_t debugPrintEnable_ms;
  uint32_t moduloSec;
  uint32_t ratesSampling;
  uint32_t chrono;
  bool sequentialCode;
  char *codeRev;
  void *thisNet;
  uint32_t globalSeed;
  //needed to describe CorticalModules
  uint32_t globCFT; //CF := Cortical Field composed of CM (Cortical Modules)
  uint32_t globCFX; //number of CM (Cortical Modules) along X
  uint32_t globCFY; //number of CM (Cortical Modules) along Y
  uint32_t neuronsPerCM; //neurons per Cortical Module (CM) 
  uint32_t RSPerCM;
  uint32_t first_locCFX, last_locCFX; // first and last CM x on this loc_h
  uint32_t first_locCFY, last_locCFY; // first and last CM y on this loc_h
  float locCFT; //number of CM on this loc_h
  uint32_t first_glob_n, last_glob_n;
  uint32_t subPopNumber;
  uint32_t tabNeuron;
  //statistical parameters
  uint32_t freeMemAtStart;
  //uint32_t totalSyn;
  char hostName[128];
  uint32_t delayMin;
  uint32_t delayMax;
  plasticityAlgoEnum plasticityAlgo;
  struct SongSTDP_paramStruct SongSTDP_param;
  uint32_t startPlasticity_ms;
  uint32_t stopPlasticity_ms;
  bool spikes_writeDisable;
  bool statSTDPevent_writeDisable;
  bool statSynPeriodicProbe_writeDisable;
  bool statRatesPerPop_writeDisable;
  uint32_t startStatFiles_s;
  uint32_t stopStatFiles_s;
  uint32_t startRatesPerPop_s;
  uint32_t stopRatesPerPop_s;
  uint32_t startSynPeriodicProbe_s;
  uint32_t stopSynPeriodicProbe_s;
  uint32_t startSTDPevent_s;
  uint32_t stopSTDPevent_s;
  float minSynWeight_f;
  float factorFloat_2_weightType;
  float factorWeightType_2_Float;
  overallConnectivityEnum overallConnectivity;
  uint32_t maxHomogeneousSide;
  uint32_t stencilDim;
  uint32_t maxModDeltaX, maxModDeltaY;
  uint32_t stencilX_Max, stencilY_Max;
  uint32_t bathEfficacyTemplateX_Max;
  uint32_t bathEfficacyTemplateY_Max;

  collectivesAlgoEnum collectivesAlgo;   // choice of ordinary A2A and A2AV or loops of Send/Recv emulating them
                                         // if collectivesAlgo == standardMPI -> ordinary A2A
                                         // if collectivesAlgo == unrolledMPI -> loop of Send/Recv
  uint32_t pktLength;  // choice of how many spikes to fit into a fixed-packet length
                       // if pktLength == 0, use ordinary algo with A2AV(lengths) + A2AV(payload)
                       // if pktLength > 0 , use A2A(packet) + A2AV(remainder)
};


#ifdef LIFCAneuron
enum neuralKindEnum {excitatoryLbExc, excitatoryLaExc, inhibitoryLaInh,
		     excitatoryRS, inhibitoryFS,
		     unknown};
struct inputCurrentArrayStruct{
  double time;
  double value;
};

class inputCurrentClass {

 public:
  float inputCurrent;
  double originalEmissionTime;
};

//enum neuSubPopEnum {Lb23Exc,La23Exc,La23Inh,neuSubPopTotal};
enum neuSubPopEnum {L11,L12,L13,
		    L21,L22,L23,
		    L31,L32,L33,
		    L41,L42,L43,
		    L51,L52,L53,
		    L61,L62,L63,
		    L71,L72,L73,
		    L81,L82,L83,
		    L91,L92,L93,
		    neuSubPopTotal};

struct neuSubPopParamStruct {
  uint32_t count;               //Total number of neurons in this SubPop
  uint32_t offset;              //Offset of the first neuoron of this SubPop in the Cortical Module
  double DMin[DSD__subPopTot];  //Minimum delay of this SubPop
  double DMax[DSD__subPopTot];  //Maximum delay of this SubPop
  double J[DSD__subPopTot];     //Synaptic efficacy of this SubPop with all the SubPops (including itself)
  double DJ[DSD__subPopTot];    //Standard deviation of J
  double JExt;                  //Synaptic efficacy of the synaptic bath of this SubPop
  double DJExt;                 //Standard deviation of JExt
  double NuExt;                 //Standard deviation of extern bath connection
  double Tau;                   //Decay time constant for the membrane potential
  double Theta;                 //The emission threshold
  double H;                     //The reset potential
  double Tarp;                  //The absolute refractory period
  double AlphaC;                //Increase of Ca concentration after the emission of a spike
  double TauC;                  //Characteristic time of Ca channel inactivation
  double gC;                    //K-current due to a unitary Ca concentration
};

struct simpleCM_neuCoordinatesStruct {
  uint32_t inCM_n;
  uint32_t cfx_n;
  uint32_t cfy_n;
  uint32_t glob_n;
  uint32_t loc_n;
  uint32_t loc_h;
  neuSubPopEnum subPop;
  neuralKindEnum neuralKind;
};

struct simpleCM_targetNeuListStruct {
  uint32_t glob_n;
  uint32_t loc_h;
  neuSubPopEnum subPop;
  neuralKindEnum neuralKind;
};

#else
enum neuralKindEnum {excitatoryRS, inhibitoryFS, 
		     excitatoryLbExc, excitatoryLaExc, inhibitoryLaInh,
		     unknown};
#endif

#endif //multiple inclusion guard
