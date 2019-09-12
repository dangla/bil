#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Functions_s    ; typedef struct Functions_s    Functions_t ;


#include "DataFile.h"
//#include "Materials.h"


extern Functions_t* (Functions_New)     (const int) ;
extern Functions_t* (Functions_Create)  (DataFile_t*) ;


#define Functions_GetNbOfFunctions(FCTS)      ((FCTS)->n_fn)
#define Functions_GetFunction(FCTS)           ((FCTS)->fn)


#include "Function.h"


struct Functions_s {          /* time functions */
  unsigned int n_fn ;         /* nb of functions */
  Function_t* fn ;            /* function */
} ;

#endif
