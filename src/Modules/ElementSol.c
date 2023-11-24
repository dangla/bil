#include <stdio.h>
#include "ElementSol.h"
#include "Message.h"
#include "Mry.h"
#include "GenericData.h"



/* Extern functions */

ElementSol_t*   (ElementSol_New)(void)
{
  ElementSol_t* elementsol = (ElementSol_t*) Mry_New(ElementSol_t) ;
  
  #if 1
  ElementSol_GetImplicitGenericData(elementsol) = GenericData_New() ;
  ElementSol_GetExplicitGenericData(elementsol) = GenericData_New() ;
  ElementSol_GetConstantGenericData(elementsol) = GenericData_New() ;
  #endif
  
  #if 0
  ElementSol_GetImplicitGenericData(elementsol) = NULL ;
  ElementSol_GetExplicitGenericData(elementsol) = NULL ;
  ElementSol_GetConstantGenericData(elementsol) = NULL ;
  #endif
  
  //ElementSol_GetPreviousElementSol(elementsol)  = NULL ;
  //ElementSol_GetNextElementSol(elementsol)      = NULL ;
  
  return(elementsol) ;
}



void   (ElementSol_Delete)(void* self)
{
  ElementSol_t*   elementsol = (ElementSol_t*) self ;
  
  {
    GenericData_t* gdat = ElementSol_GetImplicitGenericData(elementsol) ;
    
    if(gdat) {
      GenericData_Delete(gdat) ;
      free(gdat) ;
      ElementSol_GetImplicitGenericData(elementsol) = NULL ;
    }
  }
  
  {
    GenericData_t* gdat = ElementSol_GetExplicitGenericData(elementsol) ;
    
    if(gdat) {
      GenericData_Delete(gdat) ;
      free(gdat) ;
      ElementSol_GetExplicitGenericData(elementsol) = NULL ;
    }
  }
  
  {
    GenericData_t* gdat = ElementSol_GetConstantGenericData(elementsol) ;
    
    if(gdat) {
      GenericData_Delete(gdat) ;
      free(gdat) ;
      ElementSol_GetConstantGenericData(elementsol) = NULL ;
    }
  }
}



#if 0
void (ElementSol_AllocateMemoryForImplicitTerms)(ElementSol_t* elementsol)
/**  Allocate the memory space for the implicit terms
 */
{
  int ni = ElementSol_GetNbOfImplicitTerms(elementsol) ;
  double* vi = (double*) Mry_New(double[ni]) ;
  GenericData_t* gdat = GenericData_Create(ni,vi,double,"implicit terms") ;
  
  free(GenericData_GetName(ElementSol_GetImplicitGenericData(elementsol))) ;
  ElementSol_GetImplicitGenericData(elementsol)[0] = gdat[0] ;
  free(gdat) ;
}



void (ElementSol_AllocateMemoryForExplicitTerms)(ElementSol_t* elementsol)
/**  Allocate the memory space for the explicit terms
 */
{
  int ne = ElementSol_GetNbOfExplicitTerms(elementsol) ;
  double* ve = (double*) Mry_New(double[ne]) ;
  GenericData_t* gdat = GenericData_Create(ne,ve,double,"explicit terms") ;
  
  free(GenericData_GetName(ElementSol_GetExplicitGenericData(elementsol))) ;
  ElementSol_GetExplicitGenericData(elementsol)[0] = gdat[0] ;
  free(gdat) ;
}



void (ElementSol_AllocateMemoryForConstantTerms)(ElementSol_t* elementsol)
/**  Allocate the memory space for the constant terms
 */
{
  int n0 = ElementSol_GetNbOfConstantTerms(elementsol) ;
  double* v0 = (double*) Mry_New(double[n0]) ;
  GenericData_t* gdat = GenericData_Create(n0,v0,double,"constant terms") ;
  
  free(GenericData_GetName(ElementSol_GetConstantGenericData(elementsol))) ;
  ElementSol_GetConstantGenericData(elementsol)[0] = gdat[0] ;
  free(gdat) ;
}
#endif



#if 1
void (ElementSol_AllocateMemoryForImplicitTerms)(ElementSol_t* elementsol)
/**  Allocate the memory space for the implicit terms
 */
{
  GenericData_t* gdat = ElementSol_GetImplicitGenericData(elementsol) ;
  int ni = ElementSol_GetNbOfImplicitTerms(elementsol) ;
  double* vi = (double*) Mry_New(double[ni]) ;
        
  GenericData_Initialize(gdat,ni,vi,double,"implicit terms") ;
}



void (ElementSol_AllocateMemoryForExplicitTerms)(ElementSol_t* elementsol)
/**  Allocate the memory space for the explicit terms
 */
{
  GenericData_t* gdat = ElementSol_GetExplicitGenericData(elementsol) ;
  int ne = ElementSol_GetNbOfExplicitTerms(elementsol) ;
  double* ve = (double*) Mry_New(double[ne]) ;
        
  GenericData_Initialize(gdat,ne,ve,double,"explicit terms") ;
}



void (ElementSol_AllocateMemoryForConstantTerms)(ElementSol_t* elementsol)
/**  Allocate the memory space for the constant terms
 */
{
  GenericData_t* gdat = ElementSol_GetConstantGenericData(elementsol) ;
  int n0 = ElementSol_GetNbOfConstantTerms(elementsol) ;
  double* v0 = (double*) Mry_New(double[n0]) ;
        
  GenericData_Initialize(gdat,n0,v0,double,"constant terms") ;
}
#endif




void (ElementSol_Copy)(ElementSol_t* elementsol_d,ElementSol_t* elementsol_s)
/** Copy the (im/ex)plicit and constant terms 
 *  from elementsol_src to elementsol_dest */
{
  
  /* Implicit terms */
  {
    double* vi_s = (double*) ElementSol_GetImplicitTerm(elementsol_s) ;
    double* vi_d = (double*) ElementSol_GetImplicitTerm(elementsol_d) ;
    int nvi = ElementSol_GetNbOfImplicitTerms(elementsol_s) ;
    unsigned int i ;
        
    if(vi_d != vi_s) {
      for(i = 0 ; i < nvi ; i++) {
        vi_d[i] = vi_s[i] ;
      }
    }
  }

  /* Explicit terms */
  {
    double* ve_s = (double*) ElementSol_GetExplicitTerm(elementsol_s) ;
    double* ve_d = (double*) ElementSol_GetExplicitTerm(elementsol_d) ;
    int nve = ElementSol_GetNbOfExplicitTerms(elementsol_s) ;
    unsigned int i ;
      
    if(ve_d != ve_s) {
      for(i = 0 ; i < nve ; i++) {
        ve_d[i] = ve_s[i] ;
      }
    }
  }

  /* Constant terms */
  {
    double* vc_s = (double*) ElementSol_GetConstantTerm(elementsol_s) ;
    double* vc_d = (double*) ElementSol_GetConstantTerm(elementsol_d) ;
    int nvc = ElementSol_GetNbOfConstantTerms(elementsol_s) ;
    unsigned int i ;
        
    if(vc_d != vc_s) {
      for(i = 0 ; i < nvc ; i++) {
        vc_d[i] = vc_s[i] ;
      }
    }
  }
}
