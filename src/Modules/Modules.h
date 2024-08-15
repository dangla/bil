#ifndef MODULES_H
#define MODULES_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Vacuous declarations and typedef names */

/* class-like structure "Modules_t" */
struct Modules_s       ; typedef struct Modules_s       Modules_t ;


/* Declaration of Macros, Methods and Structures */

#include "Module.h"

extern void       (Modules_Delete)(void*) ;
extern void       (Modules_Print)(char*) ;
extern Module_t*  (Modules_FindModule)(Modules_t*,const char*) ;


#include "ListOfModules.h"

#define Modules_NbOfModules               (ListOfModules_Nb)
#define Modules_ListOfNames               ListOfModules_Names
#define Modules_ListOfSetModuleProp       ListOfModules_Methods(_SetModuleProp)


#define Modules_GetNbOfModules(MODS)    ((MODS)->n_modules)
#define Modules_GetModule(MODS)         ((MODS)->module)


struct Modules_s {              /* modules */
  unsigned int n_modules ;      /* nb of modules */
  Module_t* module ;            /* module */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
