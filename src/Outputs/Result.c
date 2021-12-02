#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "Result.h"
#include "View.h"




Result_t* (Result_Create)(void)
{
  Result_t* result = (Result_t*) Mry_New(Result_t) ;
  
  /* Allocate memory for the values */
  {
    double* v = (double*) Mry_New(double[9]) ;
    
    Result_GetValue(result) = v ; ;
  }
  
  /* Views */
  {
    View_t* view = View_Create() ;
    
    Result_GetView(result) = view ;
    /* Should be eliminated (used in old models) */
    Result_GetNameOfValue(result) = Result_GetNameOfView(result) ;
  }
  
  return(result) ;
}



void (Result_Delete)(void* self)
{
  Result_t* result = (Result_t*) self ;
  
  free(Result_GetValue(result)) ;
  
  {
    View_t* view = Result_GetView(result) ;
  
    View_Delete(view) ;
    free(view) ;
  }
}



void (Result_Store)(Result_t *r,double *v,const char* name,int n)
{
  double* value = Result_GetValue(r) ;
  int i ;
  
  strcpy(Result_GetNameOfView(r),name) ;
  
  Result_GetNbOfValues(r) = n ;
  Result_GetNbOfComponents(r) = n ;
  
  for(i = 0 ; i < n ; i++) {
    value[i] = v[i] ;
  }
}



void (Result_SetValuesToZero)(Result_t *r)
{
  double* v = Result_GetValue(r) ;
  int i ;
  
  for(i = 0 ; i < 9 ; i++) {
    v[i] = 0. ;
  }
}
