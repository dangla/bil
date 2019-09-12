#ifndef BCOND_H
#define BCOND_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct BCond_s        ; typedef struct BCond_s        BCond_t ;


#include "DataFile.h"


extern BCond_t*  (BCond_New)     (void) ;
extern BCond_t*  (BCond_Create)  (void) ;
extern void      (BCond_Scan)    (BCond_t*,DataFile_t*) ;


#define BCond_MaxLengthOfKeyWord        (30)

#define BCond_GetRegionIndex(BC)         ((BC)->reg)
#define BCond_GetNameOfUnknown(BC)       ((BC)->inc)
#define BCond_GetNameOfEquation(BC)      ((BC)->eqn)
#define BCond_GetFunction(BC)            ((BC)->fn)
#define BCond_GetField(BC)               ((BC)->ch)
#define BCond_GetEntryInDataFile(BC)     ((BC)->data)
#define BCond_GetFunctionIndex(BC)       ((BC)->fctindex)
#define BCond_GetFieldIndex(BC)          ((BC)->fldindex)
#define BCond_GetFunctions(BC)           ((BC)->functions)
#define BCond_GetFields(BC)              ((BC)->fields)



#include "Functions.h"
#include "Fields.h"

struct BCond_s {              /* Boundary condition */
  char*  data ;               /* Entry point in the string of data */
  int    reg ;                /* Region index */
  char*  inc ;                /* Name of unknown to be eliminated */
  char*  eqn ;                /* Name of equation to be eliminated */
  int    fctindex ;           /* Time function index */
  int    fldindex ;           /* Field index */
  Function_t* fn ;            /* Time function */
  Field_t* ch ;               /* Field */
  Functions_t* functions ;    /* Time functions */
  Fields_t* fields ;          /* Fields */
} ;


/* Old notations which I try to eliminate little by little */
#define cond_t    BCond_t

#endif
