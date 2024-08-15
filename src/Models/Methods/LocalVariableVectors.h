#ifndef LOCALVARIABLEVECTORS_H
#define LOCALVARIABLEVECTORS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct LocalVariableVectors_s     ; 
typedef struct LocalVariableVectors_s     LocalVariableVectors_t ;


/* Declaration of Macros, Methods and Structures */


extern LocalVariableVectors_t*     (LocalVariableVectors_Create)(int) ;
extern void                        (LocalVariableVectors_Delete)(void*) ;



#include "Element.h"
#include "IntFct.h"

#define LocalVariableVectors_MaxNbOfVariableVectors              (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))
//#define LocalVariableVectors_MaxNbOfVariables                    (100)



#define LocalVariableVectors_GetNbOfVariableVectors(LVV)           ((LVV)->nbofvariablevectors)
#define LocalVariableVectors_GetNbOfVariables(LVV)                 ((LVV)->nbofvariables)
#define LocalVariableVectors_GetLocalVariableVector(LVV)           ((LVV)->localvariablevector)



#define LocalVariableVectors_GetVariable(LVV,i) \
        LocalVariableVector_GetVariable(LocalVariableVectors_GetLocalVariableVector(LVV) + (i))
        
#define LocalVariableVectors_GetCurrentVariable(LVV,i) \
        LocalVariableVectors_GetVariable(LVV,i)
        
#define LocalVariableVectors_GetPreviousVariable(LVV,i) \
        LocalVariableVector_GetPreviousVariable(LocalVariableVectors_GetLocalVariableVector(LVV) + (i))
        
#define LocalVariableVectors_GetVariableDerivative(LVV,i) \
        LocalVariableVector_GetVariableDerivative(LocalVariableVectors_GetLocalVariableVector(LVV) + (i))



#include "LocalVariableVector.h"

struct LocalVariableVectors_s {
  unsigned int  nbofvariablevectors ;
  unsigned int  nbofvariables ;
  LocalVariableVector_t* localvariablevector ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
