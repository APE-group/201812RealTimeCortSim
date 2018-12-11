// localNetProcess.h
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011-...),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef localNet_H
#define localNet_H

#include "DPSNN_environmentSelection.h"
#include "DPSNN_parameters.h"
#include "DPSNN_random.h"
#include "DPSNN_neuron.h"
#include "DPSNN_messagePassing.h"
#include "DPSNN_spike.h"
#include "DPSNN_localNet.h"
#include "DPSNN_stopWatch.h"
#include "DPSNN_stat.h"

#define PORT_I "i"
#define PORT_O "o"

#ifdef DALonlyEnvironmentSelected
//#define GLOBH 8
#else
  #ifdef MPIandDALenvironmentSelected
  #define GLOBH MPIglobH
  #else
    #error "missing selection between DALonlyEnvironmentSelected and MPIandDALenvironmentSelected" 
  #endif  
#endif

#define EVENT_stop_fsm  0
#define MAX_FIRE_COUNT  2
#define MSG_SYNAPSE_INIT_COMPLETED 75640000
#define MSG_SYNC_COMPLETED 75640001
#define MSG_FIRE 75640002

class localNetClass;

typedef struct _local_states
{
  struct DPSNN_parameters parameters;
  localNetClass localNeuralNetObject; //only 1 cluster of neurons per process
  messagePassingClass messagePassing;
  stopWatchClass stopWatch;
  statClass stat;
  //several of next data could be removed once running under DAL
  int id;
  char name[32];
  int globH;
  int fireCount;
  int syncMessage;
} localNetProcess_localState;

#define lnp_ls_par p->local->parameters
#define lnp_ls_neuralNet p->local->localNeuralNetObject
#define lnp_ls_messagePassing p->local->messagePassing
#define lnp_ls_stopWatch p->local->stopWatch
#define lnp_ls_stat p->local->stat

#define EVENT_stop_fsm  0

void localNetProcess_init(DALProcess *);
int localNetProcess_fire(DALProcess *);
void localNetProcess_finish(DALProcess *);

#endif
