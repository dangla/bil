#ifndef MODULE_H
#define MODULE_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Vacuous declarations and typedef names */

/* class-like structure "Modules_t" */
struct Module_s        ; typedef struct Module_s        Module_t ;



/* Declaration of Macros, Methods and Structures */

extern Module_t*  (Module_New)(void) ;
extern void       (Module_Delete)(void*) ;
extern Module_t*  (Module_Initialize)(Module_t*,const char*) ;


#if 0
#include "Mesh.h"
#include "Solution.h"
#include "BConds.h"
#include "Solver.h"

//extern int    Module_StoreCurrentSolution(Mesh_t*,double,char*) ;
//extern int    Module_LoadCurrentSolution(Mesh_t*,double*,char*) ;
//extern void   Module_InitializeCurrentPointers(Mesh_t*,Solution_t*) ;
//extern void   Module_SetCurrentUnknownsWithBoundaryConditions(Mesh_t*,BConds_t*,double) ;
//extern void   Module_UpdateCurrentUnknowns(Mesh_t*,Solver_t*) ;
#endif


#define Module_MaxLengthOfKeyWord        (30)
#define Module_MaxLengthOfFileName       (60)
#define Module_MaxLengthOfShortTitle     (80)
#define Module_MaxLengthOfAuthorNames    (80)


#define Module_GetCodeNameOfModule(MOD)   ((MOD)->codename)
#define Module_GetShortTitle(MOD)         ((MOD)->shorttitle)
#define Module_GetNameOfAuthors(MOD)      ((MOD)->authors)
#define Module_GetSetModuleProp(MOD)      ((MOD)->setmoduleprop)
#define Module_GetComputeProblem(MOD)     ((MOD)->computeproblem)
#define Module_GetSolveProblem(MOD)       ((MOD)->solveproblem)
#define Module_GetIncrement(MOD)          ((MOD)->increment)
#define Module_GetStepForward(MOD)        ((MOD)->stepforward)
#define Module_GetInitializeProblem(MOD)  ((MOD)->initializeproblem)
#define Module_GetSequentialIndex(MOD)    ((MOD)->sequentialindex)
#define Module_GetNbOfSequences(MOD)      ((MOD)->nbofsequences)


/* Copy operations */
#define Module_CopyCodeNameOfModule(MOD,codename) \
        (strcpy(Module_GetCodeNameOfModule(MOD),codename))

#define Module_CopyShortTitle(MOD,title) \
        (strcpy(Module_GetShortTitle(MOD),title))

#define Module_CopyNameOfAuthors(MOD,authors) \
        (strcpy(Module_GetNameOfAuthors(MOD),authors))


/* Short hands */
#define Module_SetModuleProp(MOD) \
        (Module_GetSetModuleProp(MOD)) ? \
        Module_GetSetModuleProp(MOD)(MOD) : -1

#define Module_ComputeProblem(MOD,...) \
        (Module_GetComputeProblem(MOD)) ? \
        Module_GetComputeProblem(MOD)(__VA_ARGS__) : -1

#define Module_SolveProblem(MOD,...) \
        (Module_GetSolveProblem(MOD)) ? \
        Module_GetSolveProblem(MOD)(__VA_ARGS__) : -1

#define Module_Increment(MOD,...) \
        (Module_GetIncrement(MOD)) ? \
        Module_GetIncrement(MOD)(__VA_ARGS__) : -1

#define Module_StepForward(MOD,...) \
        (Module_GetStepForward(MOD)) ? \
        Module_GetStepForward(MOD)(__VA_ARGS__) : -1

#define Module_InitializeProblem(MOD,...) \
        (Module_GetInitializeProblem(MOD)) ? \
        Module_GetInitializeProblem(MOD)(__VA_ARGS__) : -1


/*  Typedef names of Methods */
#include "DataSet.h"
#include "OutputFiles.h"
#include "Solutions.h"
#include "Solver.h"

typedef int   Module_SetModuleProp_t(Module_t*) ;
typedef int   Module_ComputeProblem_t(DataSet_t*) ;
typedef int   Module_SolveProblem_t(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*) ;
typedef int   Module_StepForward_t(DataSet_t*,Solutions_t*,Solver_t*,double,double) ;
typedef int   Module_Increment_t(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*,double,double) ;
typedef int   Module_InitializeProblem_t(DataSet_t*,Solutions_t*) ;


struct Module_s {              /* module */
  Module_SetModuleProp_t*      setmoduleprop ;
  Module_ComputeProblem_t*     computeproblem ;
  Module_SolveProblem_t*       solveproblem ;
  Module_Increment_t*          increment ;
  Module_StepForward_t*        stepforward ;
  Module_InitializeProblem_t*  initializeproblem ;
  unsigned int nbofsequences ;
  int sequentialindex ;
  
  char*  codename ;          /* Code name of the module */
  char*  authors ;           /* Authors of this module */
  char*  shorttitle ;        /* Short title of the module */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
