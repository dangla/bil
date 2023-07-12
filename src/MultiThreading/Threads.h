#ifndef THREADS_H
#define THREADS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Threads_s     ; typedef struct Threads_s     Threads_t ;


//extern int Threads_CurrentThreadId(void) ;



#define Threads_GetNbOfCores(TH)         ((TH)->numberofcores)
#define Threads_GetMaxNbOfThreads(TH)    ((TH)->maxnumberofthreads)



#define  Threads_None      0
#define  Threads_Pthread   1
#define  Threads_OpenMP    2



#include "BilConfig.h"


#if defined HAVE_PTHREAD
  #include <pthread.h>
  #define Threads_API  Pthread
#elif defined HAVE_OPENMP
  #include <omp.h>
  #define Threads_API  OpenMP
#else
  #define Threads_API  None
#endif


#include "Utils.h"

/* Test the Multithreading API */
#define Threads_APIis(TYP) \
        (Utils_CAT(Threads_,Threads_API) == Utils_CAT(Threads_,TYP))

#define Threads_APIisNot(TYP) \
        (!Threads_APIis(TYP))


/* The nb of logical cores (n_sockets*n_cores*n_threads) 
 * n_sockets = number of CPU sockets 
 * n_cores   = number of physical cores per socket
 * n_threads = number of physical threads per core
 */
#include <unistd.h>
#include <sys/sysinfo.h>


#define Threads_NbOfLogicalCores \
        get_nprocs()
        //sysconf(_SC_NPROCESSORS_ONLN)



/* The max nb of threads that can be used */
#if Threads_APIis(None)
  #define Threads_MaxNbOfThreads   1
  #define Threads_DefaultNbOfRequestedThreads   1
#else
  #define Threads_MaxNbOfThreads   (4*Threads_NbOfLogicalCores)
  #define Threads_DefaultNbOfRequestedThreads   Threads_NbOfLogicalCores
#endif


/* The current thread Id */
#if Threads_APIis(Pthread) // Not yet implemented
  #define Threads_CurrentThreadId \
          #error "Pthread not available"
          //((int) pthread_self())
#elif Threads_APIis(OpenMP)
  #define Threads_CurrentThreadId \
          ((int) omp_get_thread_num())
#else
  #define Threads_CurrentThreadId    0
#endif


/* Set the number of threads */
#if Threads_APIis(Pthread)
  #define Threads_SetTheNbOfThreads(N) \
          #error "Pthread not available"
          //((int) pthread_self())
#elif Threads_APIis(OpenMP)
  #define Threads_SetTheNbOfThreads(N) \
          omp_set_num_threads(N)
#else
  #define Threads_SetTheNbOfThreads(N) 1
#endif


/* Get the number of threads in the current team */
#if Threads_APIis(Pthread)
  #define Threads_GetTheNbOfThreads() \
          #error "Pthread not available"
          //((int) pthread_self())
#elif Threads_APIis(OpenMP)
  #define Threads_GetTheNbOfThreads() \
          omp_get_num_threads()
#else
  #define Threads_GetTheNbOfThreads() 1
#endif



struct Threads_s {          /* Threads */
  int numberofcores ;
  int numberofphysicalthreads ;
  int maxnumberofthreads ;
} ;

#endif
