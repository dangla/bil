#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Message.h"
#include "Periodicity.h"



Periodicity_t* (Periodicity_New)(const int n)
{
  Periodicity_t* periodicity = (Periodicity_t*) malloc(n*sizeof(Periodicity_t)) ;
    
  assert(periodicity) ;
   
    
  if(n > 0) {
    double* vector = (double*) malloc(n*3*sizeof(double)) ;
    int i ;
    
    assert(vector) ;
    
    for(i = 0 ; i < n ; i++) {
      Periodicity_GetPeriodVector(periodicity + i) = vector + 3*i ;
    }
  }
  
  return(periodicity) ;
}




void Periodicity_Delete(void* self)
{
  Periodicity_t** pperiodicity = (Periodicity_t**) self ;
  Periodicity_t*   periodicity = *pperiodicity ;
  
  free(Periodicity_GetPeriodVector(periodicity)) ;
  free(periodicity) ;
  
  *pperiodicity = NULL ;
}
