#ifndef MODULES_H
#define MODULES_H


/* Vacuous declarations and typedef names */

/* class-like structure "Modules_t" */
struct Modules_s       ; typedef struct Modules_s       Modules_t ;
/*   1. Modules_t attributes */
struct Module_s        ; typedef struct Module_s        Module_t ;


/* Declaration of Macros, Methods and Structures */

/* 1. Modules_t */

extern Modules_t* Modules_Create(void) ;
extern void       Modules_Delete(Modules_t**) ;
extern void       Modules_Print(char*) ;
extern Module_t*  Modules_FindModule(Modules_t*,const char*) ;

#define Modules_GetNbOfModules(MODS)    ((MODS)->n_modules)
#define Modules_GetModule(MODS)         ((MODS)->module)


struct Modules_s {              /* modules */
  unsigned int n_modules ;      /* nb of modules */
  Module_t* module ;            /* module */
} ;



/* 2. Module_t */

extern Module_t*  Module_Create(int) ;

#include "Solutions.h"
#include "Mesh.h"
#include "IConds.h"
#include "BConds.h"
#include "Solver.h"

extern int    Module_StoreCurrentSolution(Mesh_t*,double,char*) ;
extern int    Module_LoadCurrentSolution(Mesh_t*,double*,char*) ;
extern void   Module_InitializeCurrentPointers(Mesh_t*,Solution_t*) ;
extern void   Module_SetCurrentUnknownsWithBoundaryConditions(Mesh_t*,BConds_t*,double) ;
extern void   Module_UpdateCurrentUnknowns(Mesh_t*,Solver_t*) ;


#define Module_MaxLengthOfKeyWord        (30)
#define Module_MaxLengthOfFileName       (60)
#define Module_MaxLengthOfShortTitle     (80)
#define Module_MaxLengthOfAuthorNames    (80)


#define Module_GetCodeNameOfModule(MOD) ((MOD)->codename)
#define Module_GetShortTitle(MOD)       ((MOD)->shorttitle)
#define Module_GetNameOfAuthors(MOD)    ((MOD)->authors)
#define Module_GetSetModuleProp(MOD)    ((MOD)->setmoduleprop)
#define Module_GetComputeProblem(MOD)   ((MOD)->computeproblem)
#define Module_GetSolveProblem(MOD)     ((MOD)->solveproblem)


/* Copy operations */
#define Module_CopyCodeNameOfModule(MOD,codename) \
        (strcpy(Module_GetCodeNameOfModule(MOD),codename))

#define Module_CopyShortTitle(MOD,title) \
        (strcpy(Module_GetShortTitle(MOD),title))

#define Module_CopyNameOfAuthors(MOD,authors) \
        (strcpy(Module_GetNameOfAuthors(MOD),authors))


/* Short hands */
#define Module_SetModuleProp(MOD) \
        Module_GetSetModuleProp(MOD)(MOD)

#define Module_ComputeProblem(MOD,...) \
        Module_GetComputeProblem(MOD)(__VA_ARGS__)

#define Module_SolveProblem(MOD,...) \
        Module_GetSolveProblem(MOD)(__VA_ARGS__)


/*  Typedef names of Methods */

#include "DataSet.h"
#include "OutputFiles.h"

typedef int   Module_SetModuleProp_t(Module_t*) ;
typedef int   Module_ComputeProblem_t(DataSet_t*) ;
typedef int   Module_SolveProblem_t(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*) ;

struct Module_s {            /* module */
  Module_SetModuleProp_t*  setmoduleprop ;
  Module_ComputeProblem_t* computeproblem ;
  Module_SolveProblem_t*   solveproblem ;
  
  char*  codename ;          /* Code name of the module */
  char*  authors ;           /* Authors of this module */
  char*  shorttitle ;        /* Short title of the module */
} ;

#endif
