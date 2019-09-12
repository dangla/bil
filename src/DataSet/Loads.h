#ifndef LOADS_H
#define LOADS_H

/* vacuous declarations and typedef names */

/* class-like structure and attributes */
struct Loads_s        ; typedef struct Loads_s        Loads_t ;



#include "DataFile.h"
#include "Functions.h"
#include "Fields.h"

extern Loads_t* Loads_Create(DataFile_t*,Fields_t*,Functions_t*) ;

#define Loads_GetNbOfLoads(LOADS)        ((LOADS)->n_cg)
#define Loads_GetLoad(LOADS)             ((LOADS)->cg)



#include "Load.h"

struct Loads_s {              /* loadings */
  unsigned int n_cg ;         /* nb */
  Load_t* cg ;                /* loading */
} ;


#endif
