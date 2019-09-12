#ifndef ICOND_H
#define ICOND_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct ICond_s        ; typedef struct ICond_s        ICond_t ;


extern ICond_t* (ICond_Create)(int) ;



#define ICond_MaxLengthOfKeyWord        (30)
#define ICond_MaxLengthOfFileName       (60)

#define ICond_GetRegionIndex(IC)             ((IC)->reg)
#define ICond_GetNameOfUnknown(IC)           ((IC)->inc)
#define ICond_GetFunction(IC)                ((IC)->fn)
#define ICond_GetField(IC)                   ((IC)->ch)
#define ICond_GetFileNameOfNodalValues(IC)   ((IC)->file)
#define ICond_GetEntryInDataFile(IC)         ((IC)->data)
#define ICond_GetFunctionIndex(IC)           ((IC)->fctindex)
#define ICond_GetFieldIndex(IC)              ((IC)->fldindex)
#define ICond_GetFunctions(IC)               ((IC)->functions)
#define ICond_GetFields(IC)                  ((IC)->fields)



#include "Functions.h"
#include "Fields.h"


struct ICond_s {              /* Initial condition */
  char*  data ;               /* Entry point in the string of data */
  int    reg ;                /* Region index */
  char*  inc ;                /* Name of unknown */
  int    fctindex ;           /* Time function index */
  int    fldindex ;           /* Field index */
  Function_t* fn ;            /* Time function */
  Field_t* ch ;               /* Field */
  char*  file ;               /* Name of file containing nodal values */
  Functions_t* functions ;    /* Time functions */
  Fields_t* fields ;          /* Fields */
} ;

#endif
