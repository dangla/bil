#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "Results.h"




Results_t* (Results_Create)(int n)
{
  Results_t* results = (Results_t*) Mry_New(Results_t) ;
  
  Results_GetNbOfResults(results) = n ;
  
  {
    Result_t* result = (Result_t*) Mry_New(Result_t[n]) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Result_t* rs = Result_Create() ;
  
      result[i] = rs[0] ;
      free(rs) ;
    }
  
    Results_GetResult(results) = result ;
  }
  
  return(results);
}



void (Results_Delete)(void* self)
{
  Results_t* results = (Results_t*) self ;
  
  {
    int n = Results_GetNbOfResults(results) ;
    Result_t* result = Results_GetResult(results) ;
    int i ;

    for(i = 0 ; i < n ; i++) {
      Result_Delete(result + i) ;
    }

    free(result) ;
  }
}
