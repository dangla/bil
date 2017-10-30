#ifndef RESULTS_H
#define RESULTS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Results_s      ; typedef struct Results_s      Results_t ;
/*     1. Results_t attributes */
struct Result_s       ; typedef struct Result_s       Result_t ;


/* 1. Results_t */
extern Results_t* Results_Create(int) ;
extern void       Results_Delete(Results_t**) ;


#define Results_GetNbOfResults(results)      ((results)->nbofresults)
#define Results_GetResult(results)           ((results)->result)



struct Results_s {            /* Results */
  unsigned int nbofresults ;  /* nb of results */
  Result_t *result ;          /* result */
} ;



/* 2. Result_t */
extern Result_t* Result_GetInstances(int) ;
extern void      Result_Store(Result_t*,double*,const char*,int) ;
extern void      Result_SetValuesToZero(Result_t*) ;



#define Result_GetValue(result)           ((result)->v)
#define Result_GetView(result)            ((result)->view)
/* These 2 next attributes should be eliminated (used in old models only) */
#define Result_GetNbOfValues(result)      ((result)->n)
#define Result_GetNameOfValue(result)     ((result)->text)


#define Result_GetNbOfComponents(result)  (View_GetNbOfComponents(Result_GetView(result)))
#define Result_GetNameOfView(result)      (View_GetNameOfView(Result_GetView(result)))


#include "Views.h"

struct Result_s {             /* Result */
  double* v ;                 /* Values */
  short int n ;               /* Nb of values (1,3,9) */
  char*   text ;              /* Name of the result */
  View_t* view ;              /* View (scalar, vector, tensor) */
} ;

/* Old notations which I try to eliminate little by little */
#define resu_t    Result_t

#endif
