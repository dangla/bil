#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "LocalVariableVector.h"
#include "Mry.h"


static void LocalVariableVector_Free(void*) ;



LocalVariableVector_t* LocalVariableVector_Create(int n)
{
  LocalVariableVector_t* lvv = (LocalVariableVector_t*) Mry_New(LocalVariableVector_t) ;
  
  
  LocalVariableVector_GetNbOfVariables(lvv) = n ;
  
  
  /* Space allocation for the variables */
  {
    double* val = (double*) Mry_New(double[n]) ;
    
    LocalVariableVector_GetVariable(lvv) = val ;
  }
  
  
  /* Space allocation for the variables at the previous time */
  {
    double* val = (double*) Mry_New(double[n]) ;
    
    LocalVariableVector_GetPreviousVariable(lvv) = val ;
  }
  
  
  /* Space allocation for the variable derivatives */
  {
    double* val = (double*) Mry_New(double[n]) ;
    
    LocalVariableVector_GetVariableDerivative(lvv) = val ;
  }
  
  return(lvv) ;
}



void LocalVariableVector_Delete(void* self,const int nvec)
{
  LocalVariableVector_t** plvv = (LocalVariableVector_t**) self ;
  LocalVariableVector_t*   lvv =  *plvv ;
  
  {
    int i ;
    
    for(i = 0 ; i < nvec ; i++) {
      LocalVariableVector_Free(lvv + i) ;
    }
  }
  
  free(lvv) ;
  *plvv = NULL ;
}



void LocalVariableVector_Free(void* self)
{
  LocalVariableVector_t* lvv = (LocalVariableVector_t*) self ;
  
  free(LocalVariableVector_GetVariable(lvv)) ;
  free(LocalVariableVector_GetPreviousVariable(lvv)) ;
  free(LocalVariableVector_GetVariableDerivative(lvv)) ;
}
