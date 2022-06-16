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



#define Options_GetDebug(OPT) \
        Options_GetPrintedInfos(OPT)
        

#include <string.h>
#include "String_.h"
        
#define Options_IsToPrintOutAtEachIteration(OPT) \
        (!strcmp(Options_GetPrintLevel(OPT),"2"))


/* Crout method */
#define Options_ResolutionMethodIsCrout(OPT) \
        Options_ResolutionMethodIs(OPT,"crout")
        
#define Options_SetResolutionMethodToCrout(OPT) \
        Options_SetResolutionMethodTo(OPT,"crout")


/* Superlu method */
#define Options_ResolutionMethodIsSuperLU(OPT) \
        Options_ResolutionMethodIs(OPT,"slu")
        
#define Options_SetResolutionMethodToSuperLU(OPT) \
        Options_SetResolutionMethodTo(OPT,"slu")


/* MA38 method */
#define Options_ResolutionMethodIsMA38(OPT) \
        Options_ResolutionMethodIs(OPT,"ma38")
        
#define Options_SetResolutionMethodToMA38(OPT) \
        Options_SetResolutionMethodTo(OPT,"ma38")
    

#define Options_GetFillFactor(OPT) \
        (((!strcmp(((char**) Context_GetSolver(Options_GetContext(OPT)))[2],"-ff")) && atof(((char**) Context_GetSolver(Options_GetContext(OPT)))[3])) ? \
        atof(((char**) Context_GetSolver(Options_GetContext(OPT)))[3]) : Options_DefaultFillFactor(OPT))
        

#define Options_DefaultFillFactor(OPT) \
        (Options_ResolutionMethodIsMA38(OPT) ? 2 : 0)


/* Implementations */

#define Options_ResolutionMethodIs(OPT,M) \
        String_Is(Options_GetResolutionMethod(OPT),M)
        
#define Options_SetResolutionMethodTo(OPT,M) \
        (strncpy(Options_GetResolutionMethod(OPT),M,Options_MaxLengthOfKeyWord))


/* Number of sequences for the sequential non iterative approach */
#define Options_GetNbOfSequences(OPT) \
        (String_Is(Options_GetModule(OPT),"SNIA") ? atoi(((char**) Context_GetUseModule(Options_GetContext(OPT)))[2]) : 1)


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
  Context_t* context ;
} ;


#endif
