#ifndef LOCALVARIABLEVECTOR_H
#define LOCALVARIABLEVECTOR_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct LocalVariableVector_s      ; 
typedef struct LocalVariableVector_s      LocalVariableVector_t ;


/* Declaration of Macros, Methods and Structures */

extern LocalVariableVector_t*     (LocalVariableVector_Create)(int) ;
extern void                       (LocalVariableVector_Delete)(void*) ;


#define LocalVariableVector_GetNbOfVariables(LVV)           ((LVV)->nbofvariables)
#define LocalVariableVector_GetVariable(LVV)                ((LVV)->variable)
#define LocalVariableVector_GetPreviousVariable(LVV)        ((LVV)->previousvariable)
#define LocalVariableVector_GetVariableDerivative(LVV)      ((LVV)->varderiv)



#define LocalVariableVector_GetCurrentVariable(LVV) \
        LocalVariableVector_GetVariable(LVV)


struct LocalVariableVector_s {
  unsigned int  nbofvariables ;
  double*   variable ;               /* Variables */
  double*   previousvariable ;       /* Variables at the previous time */
  double*   varderiv ;               /* Variable derivatives */
} ;

#endif
