#include <string.h>
#include "Message.h"
#include "Result.h"
#include "View.h"




Result_t* Result_Create(int n)
{
  Result_t* result = (Result_t*) malloc(n*sizeof(Result_t)) ;
  
  if(!result) {
    arret("Result_Create(1)") ;
  }
  
  /* Allocate memory for the values */
  {
    double* v = (double*) malloc(n*9*sizeof(double)) ;
    int i ;
    
    if(!v) {
      arret("Result_Create(2)") ;
    }
    
    for(i = 0 ; i < n ; i++) {
      Result_GetValue(result + i)   = v + 9*i ; ;
    }
  }
  
  /* Views */
  {
    View_t* view = View_Create(n) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Result_GetView(result + i) = view + i ;
      /* Should be eliminated (used in old models) */
      Result_GetNameOfValue(result + i) = Result_GetNameOfView(result + i) ;
    }
  }
  
  return(result) ;
}



void Result_Delete(Result_t** result)
{
  View_t* view = Result_GetView(*result) ;
  
  View_Delete(&view) ;
  
  free(Result_GetValue(*result)) ;
  free(*result) ;
  *result = NULL ;
}



void Result_Store(Result_t *r,double *v,const char* name,int n)
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



void Result_SetValuesToZero(Result_t *r)
{
  double* v = Result_GetValue(r) ;
  int i ;
  
  for(i = 0 ; i < 9 ; i++) {
    v[i] = 0. ;
  }
}
