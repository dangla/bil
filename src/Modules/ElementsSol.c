#include <stdio.h>
#include "ElementsSol.h"
#include "Message.h"
#include "GenericData.h"


static void   (ElementSol_Initialize)(ElementSol_t*) ;


/* Extern functions */

ElementsSol_t*   (ElementsSol_Create)(Mesh_t* mesh)
{
  ElementsSol_t* elementssol = (ElementsSol_t*) malloc(sizeof(ElementsSol_t)) ;
  
  if(!elementssol) arret("ElementsSol_Create") ;
  
  
  {
    int NbOfElements = Mesh_GetNbOfElements(mesh) ;
    ElementSol_t* elementsol = (ElementSol_t*) malloc(NbOfElements*sizeof(ElementSol_t)) ;
  
    if(!elementsol) {
      arret("ElementsSol_Create (1) : impossible d\'allouer la memoire") ;
    }

    ElementsSol_GetElementSol(elementssol)   = elementsol ;
    ElementsSol_GetNbOfElements(elementssol) = NbOfElements ;
    
    
    /* Initialization */
    {
      int i ;
      
      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_t* elementsol_i = elementsol + i ;
        
        ElementSol_Initialize(elementsol_i) ;
      }
    }
  
    ElementsSol_GetNbOfImplicitTerms(elementssol) = 0 ;
    ElementsSol_GetNbOfExplicitTerms(elementssol) = 0 ;
    ElementsSol_GetNbOfConstantTerms(elementssol) = 0 ;
  }
  
  return(elementssol) ;
}



void (ElementsSol_AllocateMemoryForImplicitTerms)(ElementsSol_t* elementssol)
/**  Allocate the memory space for the implicit terms
 *   assuming that the numbers of these terms have been initialized
 *   in "elementssol".
 */
{

  /* Compute the total nb of implicit terms */
  {
    int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
    ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
    int    i ;
    int    n_vi = 0 ;
  
    for(i = 0 ; i < NbOfElements ; i++) {
      ElementSol_t* elementsol_i = elementsol + i ;
      
      n_vi += ElementSol_GetNbOfImplicitTerms(elementsol_i) ;
    }
  
    ElementsSol_GetNbOfImplicitTerms(elementssol) = n_vi ;
  }
  
  
  /* Allocation of space for implicit terms */
  {
    int  n_vi = ElementsSol_GetNbOfImplicitTerms(elementssol) ;
  
    {
      int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
      ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
      int    i ;
      double* vi = (double*) calloc(n_vi,sizeof(double)) ;
      
      if(!vi) {
        arret("ElementsSol_AllocateMemoryForImplicitTerms") ;
      }
  
  
      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_t* elementsol_i = elementsol + i ;
        int ni = ElementSol_GetNbOfImplicitTerms(elementsol_i) ;
        GenericData_t* gdat = ElementSol_GetImplicitGenericData(elementsol_i) ;
        
        GenericData_Initialize(gdat,ni,vi,double) ;
          
        vi += ni ;
      }
    }
  }
}



void (ElementsSol_AllocateMemoryForExplicitTerms)(ElementsSol_t* elementssol)
/**  Allocate the memory space for the explicit terms
 *   assuming that the numbers of these terms have been initialized
 *   in "elementssol".
 */
{

  /* Compute the total nb of explicit terms */
  {
    int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
    ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
    int    i ;
    int    n_ve = 0 ;
  
    for(i = 0 ; i < NbOfElements ; i++) {
      ElementSol_t* elementsol_i = elementsol + i ;
      
      n_ve += ElementSol_GetNbOfExplicitTerms(elementsol_i) ;
    }
  
    ElementsSol_GetNbOfExplicitTerms(elementssol) = n_ve ;
  }
  
  
  /* Allocation of space for explicit terms */
  {
    int  n_ve = ElementsSol_GetNbOfExplicitTerms(elementssol) ;
  
    {
      int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
      ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
      int    i ;
      double* ve = (double*) calloc(n_ve,sizeof(double)) ;
      
      if(!ve) {
        arret("ElementsSol_AllocateMemoryForExplicitTerms") ;
      }

  
      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_t* elementsol_i = elementsol + i ;
        int ne = ElementSol_GetNbOfExplicitTerms(elementsol_i) ;
        GenericData_t* gdat = ElementSol_GetExplicitGenericData(elementsol_i) ;
        
        GenericData_Initialize(gdat,ne,ve,double) ;
          
        ve += ne ;
      }
    }
  }
}



void (ElementsSol_AllocateMemoryForConstantTerms)(ElementsSol_t* elementssol)
/**  Allocate the memory space for the constant terms
 *   assuming that the numbers of these terms have been initialized
 *   in "elementssol".
 */
{

  /* Compute the total nb of constant terms */
  {
    int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
    ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
    int    i ;
    int    n_v0 = 0 ;
  
    for(i = 0 ; i < NbOfElements ; i++) {
      ElementSol_t* elementsol_i = elementsol + i ;
      
      n_v0 += ElementSol_GetNbOfConstantTerms(elementsol_i) ;
    }
  
    ElementsSol_GetNbOfConstantTerms(elementssol) = n_v0 ;
  }
  
  
  /* Allocation of space for constant terms */
  {
    int  n_v0 = ElementsSol_GetNbOfConstantTerms(elementssol) ;
  
    {
      int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
      ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
      int    i ;
      double* v0 = (double*) calloc(n_v0,sizeof(double)) ;
      
      if(!v0) {
        arret("ElementsSol_AllocateMemoryForConstantTerms") ;
      }
  
  
      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_t* elementsol_i = elementsol + i ;
        int n0 = ElementSol_GetNbOfConstantTerms(elementsol_i) ;
        GenericData_t* gdat = ElementSol_GetConstantGenericData(elementsol_i) ;
        
        GenericData_Initialize(gdat,n0,v0,double) ;
          
        v0 += n0 ;
      }
    }
  }
}



void   (ElementsSol_ShareConstantTermsFrom)(ElementsSol_t* elementssol0,ElementsSol_t* elementssol)
/**  Initialize the constant terms on the basis of "elementssol0".
 */
{

  /* Compute the nb of constant terms */
  {
    int  n_v0 = ElementsSol_GetNbOfConstantTerms(elementssol0) ;
  
    ElementsSol_GetNbOfConstantTerms(elementssol) = n_v0 ;
  }
  
  
  /* Initialize the constant terms */
  {
    {
      int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
      ElementSol_t* elementsol0 = ElementsSol_GetElementSol(elementssol0) ;
      ElementSol_t* elementsol  = ElementsSol_GetElementSol(elementssol) ;
      int    i ;

      /* Elements of elementssol and elementssol0 share the same constant terms */
      for(i = 0 ; i < NbOfElements ; i++) {
        int n0 = ElementSol_GetNbOfConstantTerms(elementsol0 + i) ;
        void* v0 = ElementSol_GetConstantTerm(elementsol0 + i) ;
        GenericData_t* gdat = ElementSol_GetConstantGenericData(elementsol + i) ;
        
        GenericData_Initialize(gdat,n0,v0,double) ;
      }

    }
  }
}



void ElementsSol_Copy(ElementsSol_t* elementssol_dest,ElementsSol_t* elementssol_src)
/** Copy the (im/ex)plicit and constant terms 
 *  from elementssol_src to elementssol_dest */
{
  unsigned int nelts = ElementsSol_GetNbOfElements(elementssol_src) ;
  ElementSol_t* elementsol_s = ElementsSol_GetElementSol(elementssol_src) ;
  ElementSol_t* elementsol_d = ElementsSol_GetElementSol(elementssol_dest) ;
  
  ElementsSol_GetNbOfElements(elementssol_dest) = nelts ;
  
    /* Implicit terms */
    {
      unsigned int ni = ElementsSol_GetNbOfImplicitTerms(elementssol_src) ;
      double* vi_s = (double*) ElementSol_GetImplicitTerm(elementsol_s) ;
      double* vi_d = (double*) ElementSol_GetImplicitTerm(elementsol_d) ;
      unsigned int i ;
      
      ElementsSol_GetNbOfImplicitTerms(elementssol_dest) = ni ;
    
      for(i = 0 ; i < ni ; i++) {
        vi_d[i] = vi_s[i] ;
      }
    }

    /* Explicit terms */
    {
      unsigned int ne = ElementsSol_GetNbOfExplicitTerms(elementssol_src) ;
      double* ve_s = (double*) ElementSol_GetExplicitTerm(elementsol_s) ;
      double* ve_d = (double*) ElementSol_GetExplicitTerm(elementsol_d) ;
      unsigned int i ;
      
      ElementsSol_GetNbOfExplicitTerms(elementssol_dest) = ne ;
      
      for(i = 0 ; i < ne ; i++) {
        ve_d[i] = ve_s[i] ;
      }
    }

    /* Constant terms */
    {
      unsigned int nc = ElementsSol_GetNbOfConstantTerms(elementssol_src) ;
      double* vc_s = (double*) ElementSol_GetConstantTerm(elementsol_s) ;
      double* vc_d = (double*) ElementSol_GetConstantTerm(elementsol_d) ;
      unsigned int i ;
      
      ElementsSol_GetNbOfConstantTerms(elementssol_dest) = nc ;
      
      for(i = 0 ; i < nc ; i++) {
        vc_d[i] = vc_s[i] ;
      }
    }
}




void   (ElementSol_Initialize)(ElementSol_t* elementsol)
{
  ElementSol_GetImplicitGenericData(elementsol) = GenericData_New() ;
  ElementSol_GetExplicitGenericData(elementsol) = GenericData_New() ;
  ElementSol_GetConstantGenericData(elementsol) = GenericData_New() ;
}



ElementSol_t* (ElementSol_GetDeepElementSol)(ElementSol_t* elementsol,unsigned int depth)
{
  while(depth--) elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
  return(elementsol) ;
}
