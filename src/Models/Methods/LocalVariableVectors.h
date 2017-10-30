#ifndef LOCALVARIABLEVECTORS_H
#define LOCALVARIABLEVECTORS_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct LocalVariableVectors_s     ; 
typedef struct LocalVariableVectors_s     LocalVariableVectors_t ;
/*     1. LocalVariableVectors_t attributes */
struct LocalVariableVector_s      ; 
typedef struct LocalVariableVector_s      LocalVariableVector_t ;


/* Declaration of Macros, Methods and Structures */


/* 1. LocalVariableVectors_t */
#include "Elements.h"
#include "IntFcts.h"


extern LocalVariableVectors_t*     LocalVariableVectors_Create(int) ;
extern void LocalVariableVectors_Delete(LocalVariableVectors_t**) ;


#define LocalVariableVectors_MaxNbOfVariableVectors              (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))
//#define LocalVariableVectors_MaxNbOfVariables                    (100)


#define LocalVariableVectors_GetNbOfVariableVectors(lvv)           ((lvv)->nbofvariablevectors)
#define LocalVariableVectors_GetNbOfVariables(lvv)                 ((lvv)->nbofvariables)
#define LocalVariableVectors_GetLocalVariableVector(lvv)           ((lvv)->localvariablevector)

#define LocalVariableVectors_GetVariable(lvv,i)              (LocalVariableVector_GetVariable(LocalVariableVectors_GetLocalVariableVector(lvv) + (i)))
#define LocalVariableVectors_GetVariableDerivative(lvv,i)    (LocalVariableVector_GetVariableDerivative(LocalVariableVectors_GetLocalVariableVector(lvv) + (i)))



struct LocalVariableVectors_s {
  unsigned int  nbofvariablevectors ;
  unsigned int  nbofvariables ;
  LocalVariableVector_t* localvariablevector ;
} ;



/* 2. LocalVariableVector_t */
#define LocalVariableVector_GetVariable(lvv)                ((lvv)->variable)
#define LocalVariableVector_GetVariableDerivative(lvv)      ((lvv)->varderiv)


struct LocalVariableVector_s {                /* Local Secondary Variables */
  double*   variable ;                        /* Variable */
  double*   varderiv ;                        /* Variable derivative */
} ;

#endif
