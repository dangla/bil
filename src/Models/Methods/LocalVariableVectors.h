#ifndef LOCALVARIABLEVECTORS_H
#define LOCALVARIABLEVECTORS_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct LocalVariableVectors_s     ; 
typedef struct LocalVariableVectors_s     LocalVariableVectors_t ;


/* Declaration of Macros, Methods and Structures */


extern LocalVariableVectors_t*     LocalVariableVectors_Create(int) ;
extern void                        LocalVariableVectors_Delete(LocalVariableVectors_t**) ;



#include "Element.h"
#include "IntFct.h"

#define LocalVariableVectors_MaxNbOfVariableVectors              (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))
//#define LocalVariableVectors_MaxNbOfVariables                    (100)



#define LocalVariableVectors_GetNbOfVariableVectors(lvv)           ((lvv)->nbofvariablevectors)
#define LocalVariableVectors_GetNbOfVariables(lvv)                 ((lvv)->nbofvariables)
#define LocalVariableVectors_GetLocalVariableVector(lvv)           ((lvv)->localvariablevector)



#define LocalVariableVectors_GetVariable(lvv,i) \
        LocalVariableVector_GetVariable(LocalVariableVectors_GetLocalVariableVector(lvv) + (i))
        
#define LocalVariableVectors_GetVariableDerivative(lvv,i) \
        LocalVariableVector_GetVariableDerivative(LocalVariableVectors_GetLocalVariableVector(lvv) + (i))



#include "LocalVariableVector.h"

struct LocalVariableVectors_s {
  unsigned int  nbofvariablevectors ;
  unsigned int  nbofvariables ;
  LocalVariableVector_t* localvariablevector ;
} ;

#endif
