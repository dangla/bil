#ifndef BUFFERS_H
#define BUFFERS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct Buffers_s     ; typedef struct Buffers_s     Buffers_t ;


extern Buffers_t*  (Buffers_Create)(size_t) ;
extern void        (Buffers_Delete)(void*) ;


#define Buffers_GetNbOfBuffers(BFS)   ((BFS)->nbofbuffers)
#define Buffers_GetBuffer(BFS)        ((BFS)->buffer)

#include "SharedMS.h"

#define Buffers_MaxNbOfBuffers   SharedMS_MaxNbOfThreads


#define Buffers_GetBufferOfCurrentThread(BFS) \
        ((SharedMS_CurrentThreadId < Buffers_GetNbOfBuffers(BFS)) ? \
        (Buffers_GetBuffer(BFS) + SharedMS_CurrentThreadId) : \
        NULL)


#include "Buffer.h"

struct Buffers_s {
  int nbofbuffers ;
  Buffer_t* buffer ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
