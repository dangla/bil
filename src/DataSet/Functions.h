#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Functions_s    ; typedef struct Functions_s    Functions_t ;
/*   1. Functions_t attributes */
struct Function_s     ; typedef struct Function_s     Function_t ;

#include "DataFile.h"
#include "Materials.h"


/* 1. Functions */

extern Functions_t* (Functions_Create)(DataFile_t*,Materials_t*) ;
extern double       (Function_ComputeValue)(Function_t*,double) ;

//#define Function_ComputeValue(A,B) \
        ((A) ? (Function_ComputeValue)(A,B) : 1.)


#define Functions_GetNbOfFunctions(FCTS)      ((FCTS)->n_fn)
#define Functions_GetFunction(FCTS)           ((FCTS)->fn)



/* 2. Function */

extern Function_t*  (Function_New)(const int) ;

#define Function_MaxLengthOfFileName       (200)
#define Function_MaxLengthOfTextLine       (500)

#define Function_GetNbOfPoints(FCT)      ((FCT)->n)
#define Function_GetXValue(FCT)          ((FCT)->t)
#define Function_GetFValue(FCT)          ((FCT)->f)


struct Functions_s {          /* time functions */
  unsigned int n_fn ;         /* nb of functions */
  Function_t* fn ;            /* function */
} ;

struct Function_s {           /* fonction du temps */
  int    n ;                  /* nombre de points */
  double* t ;                 /* temps */
  double* f ;                 /* valeurs f(t) */
} ;


/* Old notations which I try to eliminate little by little */
#define fonc_t           Function_t
#define fonction(a,b)    Function_ComputeValue(&(b),(a))

#endif
