#ifndef DISTRIBUTEDMS_H
#define DISTRIBUTEDMS_H

/* vacuous declarations and typedef names */

/* class-like structure */
//struct DistributedMS_s     ; typedef struct DistributedMS_s     DistributedMS_t ;


//extern int DistributedMS_CurrentThreadId(void) ;



//#define DistributedMS_GetRankOfCallingProcess(DM)        ((DM)->rank)
//#define DistributedMS_GetNbOfProcesses(DM)               ((DM)->size)



#define  DistributedMS_None      0
#define  DistributedMS_MPI       1



#include "BilConfig.h"


#if defined HAVE_MPI
  #include <mpi.h>
  #define DistributedMS_API  MPI
#else
  #define DistributedMS_API  None
#endif


#include "Utils.h"


/* Test the DistributedMS API */
#define DistributedMS_APIis(TYP) \
        (Utils_CAT(DistributedMS_,DistributedMS_API) == Utils_CAT(DistributedMS_,TYP))

#define DistributedMS_APIisNot(TYP) \
        (!DistributedMS_APIis(TYP))


/* The number of processes */
#if DistributedMS_APIis(None)
  #define DistributedMS_NbOfProcessors 1
#elif DistributedMS_APIis(MPI)
  /* We use a C extension provided by GNU C:
   * A compound statement enclosed in parentheses may appear 
   * as an expression in GNU C.
   * (https://gcc.gnu.org/onlinedocs/gcc/Statement-Exprs.html#Statement-Exprs) */
  #define DistributedMS_NbOfProcessors \
  ({ \
    int DistributedMS_size ; \
    MPI_Comm_size(MPI_COMM_WORLD,&DistributedMS_size); \
    DistributedMS_size ; \
  })
#else
  #error "Distributed memory system not available"
#endif


/* The rank of the calling process */
#if DistributedMS_APIis(None)
  #define DistributedMS_RankOfCallingProcess 0
#elif DistributedMS_APIis(MPI)
  /* We use a C extension provided by GNU C:
   * A compound statement enclosed in parentheses may appear 
   * as an expression in GNU C.
   * (https://gcc.gnu.org/onlinedocs/gcc/Statement-Exprs.html#Statement-Exprs) */
  #define DistributedMS_RankOfCallingProcess \
  ({ \
    int DistributedMS_rank ; \
    MPI_Comm_rank(MPI_COMM_WORLD,&DistributedMS_rank); \
    DistributedMS_rank ; \
  })
#else
  #error "Distributed memory system not available"
#endif


/* Barrier */
#if DistributedMS_APIis(None)
  #define DistributedMS_Barrier
#elif DistributedMS_APIis(MPI)
  #define DistributedMS_Barrier \
          MPI_Barrier(MPI_COMM_WORLD)
#else
  #error "Distributed memory system not available"
#endif


#if 0
struct DistributedMS_s {
  int rank ;
  int size ;
} ;
#endif

#endif
