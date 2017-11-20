#ifndef INTFCTS_H
#define INTFCTS_H


/* vacuous declarations and typedef names */

/* class-like structure "DataSet_t and attributes */
struct IntFcts_s      ; typedef struct IntFcts_s      IntFcts_t ;


/* Declaration of Macros, Methods and Structures */


/* 1. IntFcts_t 
 * ------------*/
extern IntFcts_t*  (IntFcts_Create)(void) ;
extern int         (IntFcts_FindIntFct)(IntFcts_t*,int,int,const char*) ;
extern int         (IntFcts_AddIntFct)(IntFcts_t*,int,int,const char*) ;


#define IntFcts_MaxNbOfIntFcts             (4)

#define IntFcts_GetNbOfIntFcts(IFCTS)    ((IFCTS)->n_fi)
#define IntFcts_GetIntFct(IFCTS)         ((IFCTS)->fi)


#include "IntFct.h"

struct IntFcts_s {            /* interpolations */
  unsigned int n_fi ;         /* nb of interpolation function */
  IntFct_t *fi ;              /* interpolation function */
} ;


#endif
