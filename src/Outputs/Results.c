#include <string.h>
#include "Message.h"
#include "Results.h"




Results_t* Results_Create(int n)
{
  Results_t* results = (Results_t*) malloc(sizeof(Results_t)) ;
  
  if(!results) {
    arret("Results_Create") ;
  }
  
  Results_GetNbOfResults(results) = n ;
  
  Results_GetResult(results) = Result_Create(n) ;
  
  return(results);
}



void Results_Delete(void* self)
{
  Results_t** presults = (Results_t**) self ;
  Results_t*   results =  *presults;
  Result_t* result = Results_GetResult(results) ;
  
  Result_Delete(&result) ;
  free(results) ;
  *presults = NULL ;
}
