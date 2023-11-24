#ifndef SHAREDMS_H
#define SHAREDMS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct SharedMS_s     ; typedef struct SharedMS_s     SharedMS_t ;


//extern int SharedMS_CurrentThreadId(void) ;



#define SharedMS_GetNbOfCores(SM)         ((SM)->numberofcores)
#define SharedMS_GetMaxNbOfThreads(SM)    ((SM)->maxnumberofthreads)



#define  SharedMS_None      0
#define  SharedMS_Pthread   1
#define  SharedMS_OpenMP    2



#include "BilConfig.h"


#if defined HAVE_PTHREAD
  #include <pthread.h>
  #define SharedMS_API  Pthread
#elif defined HAVE_OPENMP
  #include <omp.h>
  #define SharedMS_API  OpenMP
#else
  #define SharedMS_API  None
#endif


#include "Utils.h"

/* Test the Multithreading API */
#define SharedMS_APIis(TYP) \
        (Utils_CAT(SharedMS_,SharedMS_API) == Utils_CAT(SharedMS_,TYP))

#define SharedMS_APIisNot(TYP) \
        (!SharedMS_APIis(TYP))


/* The nb of logical cores (n_sockets*n_cores*n_threads) 
 * n_sockets = number of CPU sockets 
 * n_cores   = number of physical cores per socket
 * n_threads = number of physical threads per core
 */
#include <unistd.h>
#include <sys/sysinfo.h>


#define SharedMS_NbOfLogicalCores \
        get_nprocs()
        //sysconf(_SC_NPROCESSORS_ONLN)



/* The max nb of threads that can be used */
#if SharedMS_APIis(None)
  #define SharedMS_MaxNbOfThreads   1
  #define SharedMS_DefaultNbOfRequestedThreads   1
#else
  #define SharedMS_MaxNbOfThreads   (4*SharedMS_NbOfLogicalCores)
  #define SharedMS_DefaultNbOfRequestedThreads   SharedMS_NbOfLogicalCores
#endif


/* The current thread Id */
#if SharedMS_APIis(Pthread) // Not yet implemented
  #define SharedMS_CurrentThreadId \
          #error "Pthread not available"
          //((int) pthread_self())
#elif SharedMS_APIis(OpenMP)
  #define SharedMS_CurrentThreadId \
          ((int) omp_get_thread_num())
#else
  #define SharedMS_CurrentThreadId    0
#endif


/* Set the number of threads */
#if SharedMS_APIis(Pthread)
  #define SharedMS_SetTheNbOfThreads(N) \
          #error "Pthread not available"
          //((int) pthread_self())
#elif SharedMS_APIis(OpenMP)
  #define SharedMS_SetTheNbOfThreads(N) \
          omp_set_num_threads(N)
#else
  #define SharedMS_SetTheNbOfThreads(N) 1
#endif


/* Get the number of threads in the current team */
#if SharedMS_APIis(Pthread)
  #define SharedMS_GetTheNbOfThreads() \
          #error "Pthread not available"
          //((int) pthread_self())
#elif SharedMS_APIis(OpenMP)
  #define SharedMS_GetTheNbOfThreads() \
          omp_get_num_threads()
#else
  #define SharedMS_GetTheNbOfThreads() 1
#endif



struct SharedMS_s {
  int numberofcores ;
  int numberofphysicalthreads ;
  int maxnumberofthreads ;
} ;

#endif
