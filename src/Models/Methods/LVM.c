#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "LVM.h"
#include "Message.h"
#include "Tools/Math.h"
#include "Elements.h"



LVM_t* LVM_Create(void)
{
  LVM_t* lvm = (LVM_t*) malloc(sizeof(LVM_t)) ;
  
  if(!lvm) arret("LVM_Create") ;
  
  
  /* Space allocation for the variables */
  {
    size_t sz = LVM_MaxNbOfVariableVectors*sizeof(double*) ;
    double** ptc = (double**) malloc(sz) ;
    
    if(!ptc) arret("LVM_Create") ;
    
    LVM_GetPointerToVariables(lvm) = ptc ;
  }
  {
    int i ;
    size_t sz = LVM_MaxNbOfVariableVectors*LVM_MaxNbOfVariables*sizeof(double) ;
    double* val = (double*) malloc(sz) ;
  
    if(!val) arret("LVM_Create") ;
      
    for(i = 0 ; i < LVM_MaxNbOfVariableVectors ; i++) {
      double* vali = val + i*LVM_MaxNbOfVariables ;
      
      LVM_GetPointerToVariables(lvm)[i] = vali ;
    }
  }
  
  
  /* Space allocation for the variable derivatives */
  {
    size_t sz = LVM_MaxNbOfVariableVectors*sizeof(double*) ;
    double** ptc = (double**) malloc(sz) ;
    
    if(!ptc) arret("LVM_Create") ;
    
    LVM_GetPointerToVariableDerivatives(lvm) = ptc ;
  }
  {
    int i ;
    size_t sz = LVM_MaxNbOfVariableVectors*LVM_MaxNbOfVariables*sizeof(double) ;
    double* val = (double*) malloc(sz) ;
  
    if(!val) arret("LVM_Create") ;
      
    for(i = 0 ; i < LVM_MaxNbOfVariableVectors ; i++) {
      double* vali = val + i*LVM_MaxNbOfVariables ;
      
      LVM_GetPointerToVariableDerivatives(lvm)[i] = vali ;
    }
  }
  
  
  /* Space allocation for the variable fluxes */
  {
    size_t sz = LVM_MaxNbOfVariableVectors*sizeof(double*) ;
    double** ptc = (double**) malloc(sz) ;
    
    if(!ptc) arret("LVM_Create") ;
    
    LVM_GetPointerToVariableFluxes(lvm) = ptc ;
  }
  {
    int i ;
    size_t sz = LVM_MaxNbOfVariableVectors*LVM_MaxNbOfVariableFluxes*sizeof(double) ;
    double* val = (double*) malloc(sz) ;
  
    if(!val) arret("LVM_Create") ;
      
    for(i = 0 ; i < LVM_MaxNbOfVariableVectors ; i++) {
      double* vali = val + i*LVM_MaxNbOfVariableFluxes ;
      
      LVM_GetPointerToVariableFluxes(lvm)[i] = vali ;
    }
  }
  
  
  /* Space allocation for output */
  {
    size_t sz = LVM_MaxSizeOfOutput ;
    double* output = (double*) malloc(sz) ;
    
    if(!output) arret("LVM_Create") ;
    
    LVM_GetOutput(lvm) = output ;
  }
  
  
  /* Space allocation for input */
  {
    size_t sz = LVM_MaxSizeOfInput ;
    double* input = (double*) malloc(sz) ;
    
    if(!input) arret("LVM_Create") ;
    
    LVM_GetInput(lvm)= input ;
  }
  
  
  /* Space allocation for buffer */
  {
    Buffer_t* buf = Buffer_Create(LVM_SizeOfBuffer) ;
    
    LVM_GetBuffer(lvm) = buf ;
  }
  
  return(lvm) ;
}



double* LVM_ComputeVariableDerivatives(LVM_t* lvm,LVM_ComputeSecondaryVariables_t* computesecondaryvariables,double dt,double dxi,int i,int n)
{
  Element_t* el = LVM_GetElement(lvm) ;
  Model_t* model = Element_GetModel(el) ;
  int NbOfVariables = Model_GetNbOfVariables(model) ;
  
  if(NbOfVariables > LVM_MaxNbOfVariables) {
    arret("LVM_ComputeVariableDerivatives") ;
  }
  
  {
    double* x  = LVM_GetVariables(lvm,n) ;
    double* dx = LVM_GetVariableDerivatives(lvm,n) ;
    int j ;
  
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] = x[j] ;
    }
  
    dx[i] += dxi ;
  
    computesecondaryvariables(el,dt,dx) ;
  
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] -= x[j] ;
      dx[j] /= dxi ;
    }

    return(dx) ;
  }
}




double* LVM_ComputeVariableFluxes(LVM_t* lvm,LVM_ComputeFluxes_t* computefluxes,int i,int j)
{
  Element_t* el = LVM_GetElement(lvm) ;
  Model_t* model = Element_GetModel(el) ;
  double* grdij = LVM_GetVariableDerivatives(lvm,i) ;

  {
    int NbOfVariables = Model_GetNbOfVariables(model) ;
  
    if(NbOfVariables > LVM_MaxNbOfVariables) {
      arret("LVM_ComputeVariableFluxes") ;
    }
    
    {
      double* xi  = LVM_GetVariables(lvm,i) ;
      double* xj  = LVM_GetVariables(lvm,j) ;
      int k ;
      
      for(k = 0 ; k < NbOfVariables ; k++)  {
        grdij[k] = xj[k] - xi[k] ;
      }
    }
  }
  
  /* Fluxes */
  {
    int NbOfVariableFluxes = Model_GetNbOfVariableFluxes(model) ;
  
    if(NbOfVariableFluxes > LVM_MaxNbOfVariableFluxes) {
      arret("LVM_ComputeVariableFluxes") ;
    }
    
    {
      double* w = LVM_GetVariableFluxes(lvm,i) ;
      
      computefluxes(el,grdij,w,i,j) ;
    
      return(w) ;
    }
  }
}
