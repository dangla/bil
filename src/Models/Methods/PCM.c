#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "PCM.h"
#include "Message.h"
#include "Tools/Math.h"
#include "Elements.h"


static PCM_t* instancepcm = NULL ;

static PCM_t*  PCM_Create(void) ;



PCM_t* PCM_Create(void)
{
  PCM_t* pcm = (PCM_t*) malloc(sizeof(PCM_t)) ;
  
  if(!pcm) arret("PCM_Create") ;
  
  
  /* Space allocation for the components */
  {
    size_t sz = PCM_MaxNbOfComponentVectors*sizeof(double*) ;
    double** ptc = (double**) malloc(sz) ;
    
    if(!ptc) arret("PCM_Create") ;
    
    PCM_GetPointerToComponents(pcm) = ptc ;
  }
  {
    int i ;
    size_t sz = PCM_MaxNbOfComponentVectors*PCM_MaxNbOfComponents*sizeof(double) ;
    double* val = (double*) malloc(sz) ;
  
    if(!val) arret("PCM_Create") ;
      
    for(i = 0 ; i < PCM_MaxNbOfComponentVectors ; i++) {
      double* vali = val + i*PCM_MaxNbOfComponents ;
      
      PCM_GetPointerToComponents(pcm)[i] = vali ;
    }
  }
  
  
  /* Space allocation for the component derivatives */
  {
    size_t sz = PCM_MaxNbOfComponentVectors*sizeof(double*) ;
    double** ptc = (double**) malloc(sz) ;
    
    if(!ptc) arret("PCM_Create") ;
    
    PCM_GetPointerToComponentDerivatives(pcm) = ptc ;
  }
  {
    int i ;
    size_t sz = PCM_MaxNbOfComponentVectors*PCM_MaxNbOfComponents*sizeof(double) ;
    double* val = (double*) malloc(sz) ;
  
    if(!val) arret("PCM_Create") ;
      
    for(i = 0 ; i < PCM_MaxNbOfComponentVectors ; i++) {
      double* vali = val + i*PCM_MaxNbOfComponents ;
      
      PCM_GetPointerToComponentDerivatives(pcm)[i] = vali ;
    }
  }
  
  
  /* Space allocation for the component fluxes */
  {
    size_t sz = PCM_MaxNbOfComponentVectors*sizeof(double*) ;
    double** ptc = (double**) malloc(sz) ;
    
    if(!ptc) arret("PCM_Create") ;
    
    PCM_GetPointerToComponentFluxes(pcm) = ptc ;
  }
  {
    int i ;
    size_t sz = PCM_MaxNbOfComponentVectors*PCM_MaxNbOfComponentFluxes*sizeof(double) ;
    double* val = (double*) malloc(sz) ;
  
    if(!val) arret("PCM_Create") ;
      
    for(i = 0 ; i < PCM_MaxNbOfComponentVectors ; i++) {
      double* vali = val + i*PCM_MaxNbOfComponentFluxes ;
      
      PCM_GetPointerToComponentFluxes(pcm)[i] = vali ;
    }
  }
  
  
  /* Space allocation for output */
  {
    size_t sz = PCM_MaxSizeOfOutput ;
    double* output = (double*) malloc(sz) ;
    
    if(!output) arret("PCM_Create") ;
    
    PCM_GetOutput(pcm) = output ;
  }
  
  
  /* Space allocation for input */
  {
    size_t sz = PCM_MaxSizeOfInput ;
    double* input = (double*) malloc(sz) ;
    
    if(!input) arret("PCM_Create") ;
    
    PCM_GetInput(pcm)= input ;
  }
  
  
  /* Space allocation for buffer */
  {
    Buffer_t* buf = Buffer_Create(PCM_SizeOfBuffer) ;
    
    PCM_GetBuffer(pcm) = buf ;
  }
  
  return(pcm) ;
}


PCM_t* PCM_GetInstance(Element_t* el,double** u,double** u_n,double* vim_n)
{
  if(!instancepcm) {
    instancepcm = PCM_Create() ;
  }
  
  PCM_GetElement(instancepcm) = el ;
  PCM_GetPointerToNodalUnknowns(instancepcm) = u ;
  PCM_GetPointerToPreviousNodalUnknowns(instancepcm) = u_n ;
  PCM_GetPreviousImplicitTerms(instancepcm) = vim_n ;
  
  PCM_FreeBuffer(instancepcm) ;
  
  return(instancepcm) ;
}



double* PCM_ComputeComponentDerivatives(PCM_t* pcm,PCM_ComputeSecondaryComponents_t* computesecondarycomponents,double dt,double dxi,int i,int n)
{
  Element_t* el = PCM_GetElement(pcm) ;
  Model_t* model = Element_GetModel(el) ;
  int NbOfComponents = Model_GetNbOfVariables(model) ;
  
  if(NbOfComponents > PCM_MaxNbOfComponents) {
    arret("PCM_ComputeComponentDerivatives") ;
  }
  
  {
    double* x  = PCM_GetPointerToComponents(pcm)[n] ;
    double* dx = PCM_GetPointerToComponentDerivatives(pcm)[n] ;
    int j ;
  
    for(j = 0 ; j < NbOfComponents ; j++) {
      dx[j] = x[j] ;
    }
  
    dx[i] += dxi ;
  
    computesecondarycomponents(pcm,dt,dx) ;
  
    for(j = 0 ; j < NbOfComponents ; j++) {
      dx[j] -= x[j] ;
      dx[j] /= dxi ;
    }

    return(dx) ;
  }
}




double* PCM_ComputeComponentFluxes(PCM_t* pcm,PCM_ComputeFluxes_t* computefluxes,int i,int j)
{
  Element_t* el = PCM_GetElement(pcm) ;
  Model_t* model = Element_GetModel(el) ;
  double* grdij = PCM_GetPointerToComponentDerivatives(pcm)[i] ;

  {
    int NbOfComponents = Model_GetNbOfVariables(model) ;
  
    if(NbOfComponents > PCM_MaxNbOfComponents) {
      arret("PCM_ComputeComponentFluxes") ;
    }
    
    {
      double* xi  = PCM_GetPointerToComponents(pcm)[i] ;
      double* xj  = PCM_GetPointerToComponents(pcm)[j] ;
      int k ;
      
      for(k = 0 ; k < NbOfComponents ; k++)  {
        grdij[k] = xj[k] - xi[k] ;
      }
    }
  }
  
  /* Fluxes */
  {
    int NbOfComponentFluxes = Model_GetNbOfVariableFluxes(model) ;
  
    if(NbOfComponentFluxes > PCM_MaxNbOfComponentFluxes) {
      arret("PCM_ComputeComponentFluxes") ;
    }
    
    {
      double* w = computefluxes(pcm,grdij,i,j) ;
    
      return(w) ;
    }
  }
}
