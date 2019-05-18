#ifndef CONTEXT_H
#define CONTEXT_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Context_s  ; typedef struct Context_s  Context_t ;



extern Context_t*    (Context_Create)(int,char**) ;
extern void          (Context_Delete)(void*) ;


#define Context_GetHelpOnline(CTX)         ((CTX)->help)
#define Context_GetPrintInfo(CTX)          ((CTX)->info)
#define Context_GetPrintModel(CTX)         ((CTX)->printmodel)
#define Context_GetPrintModule(CTX)        ((CTX)->printmodule)
#define Context_GetPrintUsage(CTX)         ((CTX)->printusage)
#define Context_GetReadOnly(CTX)           ((CTX)->readonly)
#define Context_GetInversePermutation(CTX) ((CTX)->invperm)
#define Context_GetPostProcessing(CTX)     ((CTX)->postprocess)
#define Context_GetSolver(CTX)             ((CTX)->solver)
#define Context_GetDebug(CTX)              ((CTX)->debug)
#define Context_GetPrintLevel(CTX)         ((CTX)->printlevel)
#define Context_GetUseModule(CTX)          ((CTX)->usemodule)
#define Context_GetGraph(CTX)              ((CTX)->graph)
#define Context_GetInputFileName(CTX)      ((CTX)->inputfilename)
#define Context_GetMiscellaneous(CTX)      ((CTX)->misc)
#define Context_GetElementOrdering(CTX)    ((CTX)->eorder)
#define Context_GetNodalOrdering(CTX)      ((CTX)->norder)
#define Context_GetCommandLine(CTX)        ((CTX)->commandline)
#define Context_GetOptions(CTX)            ((CTX)->options)
#define Context_GetTest(CTX)               ((CTX)->test)



#include "CommandLine.h"
#include "Options.h"

struct Context_s {        /* Context */
  CommandLine_t* commandline ;
  Options_t*     options ;
  void*   help ;
  void*   info ;
  void*   printmodel ;
  void*   printmodule ;
  void*   printusage ;
  void*   readonly ;
  void*   invperm ;
  void*   postprocess ;
  void*   solver ;
  void*   usemodule ;
  void*   debug ;
  void*   printlevel ;
  void*   graph ;
  void*   inputfilename ;
  void*   misc ;
  void*   eorder ;
  void*   norder ;
  void*   test ;
} ;


#endif
