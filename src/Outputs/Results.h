#ifndef RESULTS_H
#define RESULTS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Results_s      ; typedef struct Results_s      Results_t ;


extern Results_t* Results_Create(int) ;
extern void       Results_Delete(Results_t**) ;


#define Results_GetNbOfResults(results)      ((results)->nbofresults)
#define Results_GetResult(results)           ((results)->result)


#include "Result.h"

struct Results_s {            /* Results */
  unsigned int nbofresults ;  /* nb of results */
  Result_t *result ;          /* result */
} ;

#endif
