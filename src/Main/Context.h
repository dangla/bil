#ifndef CONTEXT_H
#define CONTEXT_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Context_s  ; typedef struct Context_s  Context_t ;


extern Context_t*    (Context_GetInstance)(void) ;


#define Context_GetHelpOnline(CTX)         ((CTX)->help)
#define Context_GetPrintInfo(CTX)          ((CTX)->info)
#define Context_GetPrintModel(CTX)         ((CTX)->printmodel)
#define Context_GetPrintModule(CTX)        ((CTX)->printmodule)
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


#define Context_IsHelpOnline \
        Context_GetHelpOnline(Context_GetInstance())
        
#define Context_IsPrintInfo \
        Context_GetPrintInfo(Context_GetInstance())

#define Context_IsPrintModel \
        Context_GetPrintModel(Context_GetInstance())

#define Context_IsPrintModule \
        Context_GetPrintModule(Context_GetInstance())

#define Context_IsReadOnly \
        Context_GetReadOnly(Context_GetInstance())

#define Context_IsGraph \
        Context_GetGraph(Context_GetInstance())

#define Context_IsInversePermutation \
        Context_GetInversePermutation(Context_GetInstance())

#define Context_IsPostProcessing \
        Context_GetPostProcessing(Context_GetInstance())

#define Context_IsMiscellaneous \
        Context_GetMiscellaneous(Context_GetInstance())

#define Context_IsElementOrdering \
        Context_GetElementOrdering(Context_GetInstance())

#define Context_IsNodalOrdering \
        Context_GetNodalOrdering(Context_GetInstance())


#define Context_HasInputFileName \
        (Context_GetInputFileName(Context_GetInstance()))

#define Context_HasSolverName \
        Context_GetSolver(Context_GetInstance())


#define Context_InputFileName \
        ((char**) Context_GetInputFileName(Context_GetInstance()))[0]
        
#define Context_SolverName \
        ((char**) Context_GetSolver(Context_GetInstance()))[0]




struct Context_s {        /* options */
  void*   help ;
  void*   info ;
  void*   printmodel ;
  void*   printmodule ;
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
} ;


#endif
