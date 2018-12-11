/*! \file dal.h
    \brief Define the DAL process handler scheme.

    - Local variables are defined in structure LocalState. Local
      variables may vary from different processes.
    - The ProcessInit function pointer points to a function which
      initializes a process.
    - The ProcessFire function pointer points to a function which
      performs the actual computation. The communication between
      processes is inside the ProcessFire function.
    - The WPTR is a placeholder for callback. One can just
      leave it blank.

    \authors Lothar Thiele, Lars Schor, Devendra Rai
*/

/************************************************************************
 * do not add code to this header
 ************************************************************************/

#ifndef DAL_H
#define DAL_H

//structure for local memory of process
typedef struct _local_states *LocalState;

//process handler
struct _process;

//additional behavioral functions could be declared here
typedef void (*ProcessInit)(struct _process*);
typedef void (*ProcessFinish)(struct _process*);
typedef int (*ProcessFire)(struct _process*);
typedef void *WPTR;

typedef struct _process {
    LocalState     local;
    ProcessInit    init;
    ProcessFire    fire;
    ProcessFinish  finish;
    WPTR           wptr; //placeholder for wrapper instance
} DALProcess;

//process interfaces
extern int DAL_write(void *port, void *buf, int len, DALProcess *p);
extern int DAL_read(void *port, void *buf, int len, DALProcess *p);
extern int DAL_send_event(void *message, DALProcess *p);
extern int *createPort(int *port, int base, int number_of_indices, int index_range_pairs, ...);
extern int DAL_getIndex(int dimension, DALProcess *p);
extern int DAL_printf(const char *in, ...);

// checkpointing interface
extern void DAL_create_checkpoint(void *buf, int len, char* name, DALProcess *p);
extern int DAL_read_checkpoint(void *buf, int len, char* name, DALProcess *p);
extern void DAL_destroy_checkpoint(char* name, DALProcess *p);

#define CREATEPORTVAR(name) \
    int name

#define CREATEPORT(port, base, number_of_indices, index_range_pairs...) \
    createPort(&port, base, number_of_indices, index_range_pairs)

#define GETINDEX(dimension) \
		DAL_getIndex(dimension, p)

#endif
