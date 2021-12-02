#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Mry.h"
#include "Message.h"
#include "Periodicity.h"



Periodicity_t* (Periodicity_New)(void)
{
  Periodicity_t* periodicity = (Periodicity_t*) Mry_New(Periodicity_t) ;

  {
    double* vector = (double*) Mry_New(double[3]) ;
    
    Periodicity_GetPeriodVector(periodicity) = vector ;
  }
  
  return(periodicity) ;
}



void (Periodicity_Delete)(void* self)
{
  Periodicity_t* periodicity = (Periodicity_t*) self ;
  
  free(Periodicity_GetPeriodVector(periodicity)) ;
}
