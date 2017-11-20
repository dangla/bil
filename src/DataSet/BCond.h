#ifndef BCOND_H
#define BCOND_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct BCond_s        ; typedef struct BCond_s        BCond_t ;


#define BCond_MaxLengthOfKeyWord        (30)

#define BCond_GetRegionIndex(BC)         ((BC)->reg)
#define BCond_GetNameOfUnknown(BC)       ((BC)->inc)
#define BCond_GetNameOfEquation(BC)      ((BC)->eqn)
#define BCond_GetFunction(BC)            ((BC)->fn)
#define BCond_GetField(BC)               ((BC)->ch)


#include "Function.h"
#include "Field.h"


struct BCond_s {              /* condition a la limite */
  int    reg ;                /* numero de la region */
  char*  inc ;                /* nom de l'inconnue a imposer */
  char*  eqn ;                /* nom de l'equation a eliminer */
  Function_t* fn ;            /* fonction du temps */
  Field_t* ch ;               /* champ */
} ;


/* Old notations which I try to eliminate little by little */
#define cond_t    BCond_t

#endif
