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



LocalVariableVectors_t* (LocalVariableVectors_Create)(int nvar)
{
  LocalVariableVectors_t* lvvs = (LocalVariableVectors_t*) Mry_New(LocalVariableVectors_t) ;

  {
    int nvec = LocalVariableVectors_MaxNbOfVariableVectors ;
    LocalVariableVector_t* lvv = (LocalVariableVector_t*) Mry_New(LocalVariableVector_t[nvec]) ;
    int i ;
    
    for(i = 0 ; i < nvec ; i++) {
      LocalVariableVector_t* lvvi = LocalVariableVector_Create(nvar) ;
      
      lvv[i] = lvvi[0] ;
      free(lvvi) ;
    }
    
    LocalVariableVectors_GetLocalVariableVector(lvvs) = lvv ;
    LocalVariableVectors_GetNbOfVariableVectors(lvvs) = nvec ;
    LocalVariableVectors_GetNbOfVariables(lvvs)       = nvar ;
  }
  
  return(lvvs) ;
}



void (LocalVariableVectors_Delete)(void* self)
{
  LocalVariableVectors_t* lvvs = (LocalVariableVectors_t*) self ;
  
  {
    int nvec = LocalVariableVectors_GetNbOfVariableVectors(lvvs) ;
    LocalVariableVector_t* lvv = LocalVariableVectors_GetLocalVariableVector(lvvs) ;
    int i ;
    
    for(i = 0 ; i < nvec ; i++) {
      LocalVariableVector_t* lvvi = lvv + i ;
      
      LocalVariableVector_Delete(lvvi) ;
    }
    
    free(lvv) ;
  }
}
