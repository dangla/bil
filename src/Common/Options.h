#ifndef OPTIONS_H
#define OPTIONS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Options_s      ; typedef struct Options_s      Options_t ;


#include "Context.h"

extern Options_t*  (Options_Create)(Context_t*) ;
extern void        (Options_Delete)(void*) ;


#define Options_MaxLengthOfKeyWord               (30)

#define Options_GetPrintedInfos(OPT)           ((OPT)->debug)
#define Options_GetResolutionMethod(OPT)       ((OPT)->method)
#define Options_GetPrintLevel(OPT)             ((OPT)->level)
#define Options_GetModule(OPT)                 ((OPT)->module)
#define Options_GetGraphMethod(OPT)            ((OPT)->graph)
#define Options_GetElementOrderingMethod(OPT)  ((OPT)->eordering)
#define Options_GetNodalOrderingMethod(OPT)    ((OPT)->nordering)
#define Options_GetPostProcessingMethod(OPT)   ((OPT)->postprocess)
#define Options_GetContext(OPT)                ((OPT)->context)
#define Options_GetNbOfThreads(OPT)            ((OPT)->nthreads)



#define Options_GetDebug(OPT) \
        Options_GetPrintedInfos(OPT)
        

#include <string.h>
#include "String_.h"
        
#define Options_IsToPrintOutAtEachIteration(OPT) \
        (!strcmp(Options_GetPrintLevel(OPT),"2"))


/* Resolution method */
#define Options_ResolutionMethodIs(OPT,...) \
        String_Is(Options_GetResolutionMethod(OPT),__VA_ARGS__)
        
//#define Options_SetResolutionMethodTo(OPT,M) \
        (strncpy(Options_GetResolutionMethod(OPT),M,Options_MaxLengthOfKeyWord))


/* Module */
#define Options_ModuleIs(OPT,M) \
        String_Is(Options_GetModule(OPT),M)


/* Crout method */
#define Options_ResolutionMethodIsCrout(OPT) \
        Options_ResolutionMethodIs(OPT,"crout")


/* SuperLU method */
#define Options_ResolutionMethodIsSuperLU(OPT) \
        (Options_ResolutionMethodIs(OPT,"slu") || \
        (Options_ResolutionMethodIs(OPT,"superlu") && \
        !Options_ResolutionMethodIs(OPT,"superlumt") && \
        !Options_ResolutionMethodIs(OPT,"superludist")))


/* SuperLUMT method */
#define Options_ResolutionMethodIsSuperLUMT(OPT) \
        Options_ResolutionMethodIs(OPT,"superlumt")


/* SuperLUDist method */
#define Options_ResolutionMethodIsSuperLUDist(OPT) \
        Options_ResolutionMethodIs(OPT,"superludist")


/* MA38 method */
#define Options_ResolutionMethodIsMA38(OPT) \
        Options_ResolutionMethodIs(OPT,"ma38")


/* PetscKSP method */
#define Options_ResolutionMethodIsPetscKSP(OPT) \
        Options_ResolutionMethodIs(OPT,"petscksp")
    

#define Options_GetFillFactor(OPT) \
        (((!strcmp(((char**) Context_GetSolver(Options_GetContext(OPT)))[2],"-ff")) && atof(((char**) Context_GetSolver(Options_GetContext(OPT)))[3])) ? \
        atof(((char**) Context_GetSolver(Options_GetContext(OPT)))[3]) : Options_DefaultFillFactor(OPT))
        

#define Options_DefaultFillFactor(OPT) \
        (Options_ResolutionMethodIsMA38(OPT) ? 2 : 0)



/* Number of solutions for the monolithic approach */
#define Options_GetNbOfSolutions(OPT) \
        ((Options_ModuleIs(OPT,"Monolithic") && Context_IsUseModule(Options_GetContext(OPT))) ? \
        atoi(((char**) Context_GetUseModule(Options_GetContext(OPT)))[2]) : 2)


/* Number of sequences for the sequential non iterative approach */
#define Options_GetNbOfSequences(OPT) \
        ((Options_ModuleIs(OPT,"SNIA") && Context_IsUseModule(Options_GetContext(OPT))) ? \
        atoi(((char**) Context_GetUseModule(Options_GetContext(OPT)))[2]) : 1)


/* Number of threads for multi-threaded calculations */
#define Options_NbOfThreads(OPT) \
        atoi(Options_GetNbOfThreads(OPT))


/* Number of threads for multi-threaded solver */
#define Options_NbOfThreadsInSolver(OPT) \
        Options_NbInResolutionMethod(OPT)


/* Number of threads for multi-threaded solver */
#define Options_NbOfThreadsInSharedMemorySolver(OPT) \
        Options_NbInResolutionMethod(OPT)


/* Number of processors for distributed-memory solver */
#define Options_NbOfProcessorsInDistributedMemorySolver(OPT) \
        Options_NbInResolutionMethod(OPT)


/* Implementation */
#define Options_NbInResolutionMethod(OPT) \
        ((String_pchar = String_FindAnyChar(Options_GetResolutionMethod(OPT),"0123456789")) ? atoi(String_pchar) : 1)



struct Options_s {            /* options */
  char*   debug ;             /* what to be printed */
  char*   method ;            /* resolution method */
  /* Print level
   *   0:  minimal, 
   *   1:  normal at each step , 
   *   2:  normal at each iteration, 
   *   3-: not used) */
  char*   level ;             /* print level */
  char*   module ;            /* module */
  char*   graph ;             /* Graph method */
  char*   eordering ;         /* Element ordering method */
  char*   nordering ;         /* Nodal ordering method */
  char*   postprocess ;       /* Post-processing method */
  char*   nthreads ;          /* Nb of requested threads */
  Context_t* context ;
} ;


#endif
