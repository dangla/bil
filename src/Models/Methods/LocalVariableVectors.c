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



LocalVariableVectors_t* LocalVariableVectors_Create(int nvar)
{
  LocalVariableVectors_t* lvvs = (LocalVariableVectors_t*) Mry_New(LocalVariableVectors_t) ;
  
  
  {
    int nvec = LocalVariableVectors_MaxNbOfVariableVectors ;
    LocalVariableVector_t* lvv = (LocalVariableVector_t*) Mry_New(LocalVariableVector_t[nvec]) ;
    int i ;
    
    for(i = 0 ; i < nvec ; i++) {
      LocalVariableVector_t* lvvi = LocalVariableVector_Create(nvar) ;
      
      lvv[i] = lvvi[0] ;
    }
    
    LocalVariableVectors_GetLocalVariableVector(lvvs) = lvv ;
    LocalVariableVectors_GetNbOfVariableVectors(lvvs) = nvec ;
    LocalVariableVectors_GetNbOfVariables(lvvs)       = nvar ;
  }
  
  
  #if 0
  {
    int nvec = LocalVariableVectors_MaxNbOfVariableVectors ;
    LocalVariableVector_t* lvv = LocalVariableVector_Create(nvec,nvar) ;
    
    LocalVariableVectors_GetLocalVariableVector(lvvs) = lvv ;
  }
  #endif
  
  return(lvvs) ;
}



void LocalVariableVectors_Delete(void* self)
{
  LocalVariableVectors_t** plvvs = (LocalVariableVectors_t**) self ;
  LocalVariableVectors_t*   lvvs =  *plvvs ;
  
  {
    int nvec = LocalVariableVectors_GetNbOfVariableVectors(lvvs) ;
    LocalVariableVector_t* lvv = LocalVariableVectors_GetLocalVariableVector(lvvs) ;
    
    LocalVariableVector_Delete(&lvv,nvec) ;
  }
  
  free(lvvs) ;
  *plvvs = NULL ;
}





#if 0
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
#endif
