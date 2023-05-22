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
