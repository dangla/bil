#ifndef LOCALVARIABLEVECTOR_H
#define LOCALVARIABLEVECTOR_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct LocalVariableVector_s      ; 
typedef struct LocalVariableVector_s      LocalVariableVector_t ;


/* Declaration of Macros, Methods and Structures */


#define LocalVariableVector_GetVariable(lvv)                ((lvv)->variable)
#define LocalVariableVector_GetVariableDerivative(lvv)      ((lvv)->varderiv)



struct LocalVariableVector_s {                /* Local Secondary Variables */
  double*   variable ;                        /* Variable */
  double*   varderiv ;                        /* Variable derivative */
} ;

#endif
