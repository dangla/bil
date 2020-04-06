#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Mry.h"
#include "Message.h"
#include "Periodicity.h"



Periodicity_t* (Periodicity_New)(const int n)
{
  Periodicity_t* periodicity = (Periodicity_t*) Mry_New(Periodicity_t[n]) ;
   
    
  if(n > 0) {
    double* vector = (double*) Mry_New(double[n*3]) ;
    int i ;
    
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
  
  //*pperiodicity = NULL ;
}
