#ifndef OPTIONS_H
#define OPTIONS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Options_s      ; typedef struct Options_s      Options_t ;


#include "Context.h"

extern Options_t*  (Options_Create)(Context_t*) ;
extern void        (Options_Delete)(void*) ;


#define Options_MaxLengthOfKeyWord               (30)

#define Options_GetPrintData(OPT)              ((OPT)->debug)
#define Options_GetResolutionMethod(OPT)       ((OPT)->method)
#define Options_GetPrintLevel(OPT)             ((OPT)->level)
#define Options_GetModule(OPT)                 ((OPT)->module)
#define Options_GetGraphMethod(OPT)            ((OPT)->graph)
#define Options_GetElementOrderingMethod(OPT)  ((OPT)->eordering)
#define Options_GetNodalOrderingMethod(OPT)    ((OPT)->nordering)
#define Options_GetPostProcessingMethod(OPT)   ((OPT)->postprocess)
#define Options_GetContext(OPT)                ((OPT)->context)



#define Options_GetDebug(OPT) \
        Options_GetPrintData(OPT)
        

#include <string.h>
        
#define Options_IsToPrintOutAtEachIteration(OPT) \
        (!strcmp(Options_GetPrintLevel(OPT),"2"))


#define Options_GetFillFactor(OPT) \
        (((!strcmp(((char**) Context_GetSolver(Options_GetContext(OPT)))[2],"-ff")) && atof(((char**) Context_GetSolver(Options_GetContext(OPT)))[3])) ? \
        atof(((char**) Context_GetSolver(Options_GetContext(OPT)))[3]) : Options_DefaultFillFactor)
        

#define Options_DefaultFillFactor (2)
        


struct Options_s {            /* options */
  char*   debug ;             /* data to be printed */
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
