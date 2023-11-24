#include <stdio.h>
#include "ElementsSol.h"
#include "Message.h"
#include "GenericData.h"
#include "Mry.h"


/* Extern functions */

ElementsSol_t*   (ElementsSol_Create)(Mesh_t* mesh)
{
  ElementsSol_t* elementssol = (ElementsSol_t*) Mry_New(ElementsSol_t) ;
  
  
  {
    int NbOfElements = Mesh_GetNbOfElements(mesh) ;
    ElementSol_t* elementsol = (ElementSol_t*) Mry_New(ElementSol_t[NbOfElements]) ;

    ElementsSol_GetElementSol(elementssol)   = elementsol ;
    ElementsSol_GetNbOfElements(elementssol) = NbOfElements ;
    
    
    /* Initialization */
    {
      int i ;
      
      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_t* elementsol_i = ElementSol_New() ;
        
        elementsol[i] = elementsol_i[0] ;
        free(elementsol_i) ;
      }
    }
  }
  
  return(elementssol) ;
}



void (ElementsSol_Delete)(void* self)
{
  ElementsSol_t* elementssol = (ElementsSol_t*) self ;
  
  {
    int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
    ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
    
    if(elementsol) {
      int i ;
      
      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_t* elementsol_i = elementsol + i ;
      
        ElementSol_Delete(elementsol_i) ;
      }
    
      free(elementsol) ;
      ElementsSol_GetElementSol(elementssol) = NULL ;
    }
  }
}



void (ElementsSol_AllocateMemoryForImplicitTerms)(ElementsSol_t* elementssol)
/**  Allocate the memory space for the implicit terms
 *   assuming that the numbers of these terms have been initialized
 *   in "elementssol".
 */
{
  int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
  ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
  
  {
    int    i ;

    for(i = 0 ; i < NbOfElements ; i++) {
      ElementSol_AllocateMemoryForImplicitTerms(elementsol + i) ;
    }
  }
}



void (ElementsSol_AllocateMemoryForExplicitTerms)(ElementsSol_t* elementssol)
/**  Allocate the memory space for the explicit terms
 *   assuming that the numbers of these terms have been initialized
 *   in "elementssol".
 */
{
  int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
  ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
  
    {
      int    i ;

      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_AllocateMemoryForExplicitTerms(elementsol + i) ;
      }
    }
}



void (ElementsSol_AllocateMemoryForConstantTerms)(ElementsSol_t* elementssol)
/**  Allocate the memory space for the constant terms
 *   assuming that the numbers of these terms have been initialized
 *   in "elementssol".
 */
{
  int NbOfElements = ElementsSol_GetNbOfElements(elementssol) ;
  ElementSol_t* elementsol = ElementsSol_GetElementSol(elementssol) ;
  
    {
      int    i ;

      for(i = 0 ; i < NbOfElements ; i++) {
        ElementSol_AllocateMemoryForConstantTerms(elementsol + i) ;
      }
    }
}



void (ElementsSol_Copy)(ElementsSol_t* elementssol_dest,ElementsSol_t* elementssol_src)
/** Copy the (im/ex)plicit and constant terms 
 *  from elementssol_src to elementssol_dest */
{
  unsigned int nelts = ElementsSol_GetNbOfElements(elementssol_src) ;
  ElementSol_t* elementsol_s = ElementsSol_GetElementSol(elementssol_src) ;
  ElementSol_t* elementsol_d = ElementsSol_GetElementSol(elementssol_dest) ;
  
  ElementsSol_GetNbOfElements(elementssol_dest) = nelts ;
  
    /* Implicit terms */
    {
      int ie ;
        
      for(ie = 0 ; ie < nelts ; ie++) {
        ElementSol_t* elementsoli_s = elementsol_s + ie ;
        ElementSol_t* elementsoli_d = elementsol_d + ie ;
        double* vi_s = (double*) ElementSol_GetImplicitTerm(elementsoli_s) ;
        double* vi_d = (double*) ElementSol_GetImplicitTerm(elementsoli_d) ;
        int nvi = ElementSol_GetNbOfImplicitTerms(elementsoli_s) ;
        unsigned int i ;
        
        if(vi_d != vi_s) {
          for(i = 0 ; i < nvi ; i++) {
            vi_d[i] = vi_s[i] ;
          }
        }
      }
    }

    /* Explicit terms */
    {
      int ie ;
        
      for(ie = 0 ; ie < nelts ; ie++) {
        ElementSol_t* elementsoli_s = elementsol_s + ie ;
        ElementSol_t* elementsoli_d = elementsol_d + ie ;
        double* ve_s = (double*) ElementSol_GetExplicitTerm(elementsoli_s) ;
        double* ve_d = (double*) ElementSol_GetExplicitTerm(elementsoli_d) ;
        int nve = ElementSol_GetNbOfExplicitTerms(elementsoli_s) ;
        unsigned int i ;
      
        if(ve_d != ve_s) {
          for(i = 0 ; i < nve ; i++) {
            ve_d[i] = ve_s[i] ;
          }
        }
      }
    }

    /* Constant terms */
    {
      int ie ;
        
      for(ie = 0 ; ie < nelts ; ie++) {
        ElementSol_t* elementsoli_s = elementsol_s + ie ;
        ElementSol_t* elementsoli_d = elementsol_d + ie ;
        double* vc_s = (double*) ElementSol_GetConstantTerm(elementsoli_s) ;
        double* vc_d = (double*) ElementSol_GetConstantTerm(elementsoli_d) ;
        int nvc = ElementSol_GetNbOfConstantTerms(elementsoli_s) ;
        unsigned int i ;
        
        if(vc_d != vc_s) {
          for(i = 0 ; i < nvc ; i++) {
            vc_d[i] = vc_s[i] ;
          }
        }
      }
    }
}
