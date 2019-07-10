#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "LocalVariableVectors.h"
#include "Message.h"
#include "Tools/Math.h"
#include "Mry.h"



static LocalVariableVector_t*     LocalVariableVector_Create(int,int) ;
static void                       LocalVariableVector_Delete(void*) ;


LocalVariableVectors_t* LocalVariableVectors_Create(int NbOfVariables)
{
  LocalVariableVectors_t* lvvs = (LocalVariableVectors_t*) Mry_New(LocalVariableVectors_t) ;
  
  //LocalVariableVectors_GetNbOfVariableVectors(lvvs) = 0 ;
  
  LocalVariableVectors_GetNbOfVariables(lvvs) = NbOfVariables ;
  
  {
    int nvec = LocalVariableVectors_MaxNbOfVariableVectors ;
    int nvar = NbOfVariables ;
    LocalVariableVector_t* lvv = LocalVariableVector_Create(nvec,nvar) ;
    
    LocalVariableVectors_GetLocalVariableVector(lvvs) = lvv ;
  }
  
  return(lvvs) ;
}



void LocalVariableVectors_Delete(void* self)
{
  LocalVariableVectors_t** plvvs = (LocalVariableVectors_t**) self ;
  LocalVariableVectors_t*   lvvs =  *plvvs ;
  
  {
    LocalVariableVector_t* lvv = LocalVariableVectors_GetLocalVariableVector(lvvs) ;
    
    LocalVariableVector_Delete(&lvv) ;
  }
  
  free(lvvs) ;
  *plvvs = NULL ;
}






LocalVariableVector_t* LocalVariableVector_Create(int NbOfVectors,int NbOfVariables)
{
  LocalVariableVector_t* lvv = (LocalVariableVector_t*) Mry_New(LocalVariableVector_t,NbOfVectors) ;
  
  
  /* Space allocation for the variables */
  {
    int i ;
    double* val = (double*) Mry_New(double,NbOfVectors*NbOfVariables) ;
    
    for(i = 0 ; i < NbOfVectors ; i++) {
      LocalVariableVector_GetVariable(lvv + i) = val + i*NbOfVariables ;
    }
  }
  
  
  /* Space allocation for the variable derivatives */
  {
    int i ;
    double* val = (double*) Mry_New(double,NbOfVectors*NbOfVariables) ;
    
    for(i = 0 ; i < NbOfVectors ; i++) {
      LocalVariableVector_GetVariableDerivative(lvv + i) = val + i*NbOfVariables ;
    }
  }
  
  return(lvv) ;
}



void LocalVariableVector_Delete(void* self)
{
  LocalVariableVector_t** plvv = (LocalVariableVector_t**) self ;
  LocalVariableVector_t*   lvv =  *plvv;
  
  free(LocalVariableVector_GetVariable(lvv)) ;
  free(LocalVariableVector_GetVariableDerivative(lvv)) ;
  
  free(lvv) ;
  *plvv = NULL ;
}
