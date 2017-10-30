#ifndef LOADS_H
#define LOADS_H

/* vacuous declarations and typedef names */

/* class-like structure and attributes */
struct Loads_s        ; typedef struct Loads_s        Loads_t ;
/*   1. Loads_t attributes */
struct Load_s         ; typedef struct Load_s         Load_t ;



/* 1. Loads_t 
 * ----------*/
#include "DataFile.h"
#include "Functions.h"
#include "Fields.h"

extern Loads_t* Loads_Create(DataFile_t*,Fields_t*,Functions_t*) ;

#define Loads_GetNbOfLoads(LOADS)        ((LOADS)->n_cg)
#define Loads_GetLoad(LOADS)             ((LOADS)->cg)



/* 2. Load_t 
 * ---------*/
#include "Elements.h"
#include "IntFcts.h"


#define Load_MaxLengthOfKeyWord               (30)

#define Load_GetRegionIndex(LOAD)        ((LOAD)->reg)
#define Load_GetType(LOAD)               ((LOAD)->t)
#define Load_GetNameOfEquation(LOAD)     ((LOAD)->eqn)
#define Load_GetFunction(LOAD)           ((LOAD)->fn)
#define Load_GetField(LOAD)              ((LOAD)->ch)

#define Load_TypeIs(LOAD,TYPE)           (!strcmp(Load_GetType(LOAD),TYPE))


struct Loads_s {              /* loadings */
  unsigned int n_cg ;         /* nb */
  Load_t *cg ;                /* loading */
} ;

struct Load_s {               /* chargement */
  int    reg ;                /* numero de la region */
  char   *t ;                 /* type de chargement */
  char   *eqn ;               /* nom de l'equation */
  Function_t *fn ;            /* fonction du temps */
  Field_t *ch ;               /* champ */
} ;


/* Old notations which I try to eliminate little by little */
#define char_t    Load_t

#endif
