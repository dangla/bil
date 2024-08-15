#ifndef LOAD_H
#define LOAD_H

/* vacuous declarations and typedef names */

/* class-like structure and attributes */
struct Load_s         ; typedef struct Load_s         Load_t ;



#include "DataFile.h"

extern Load_t* (Load_New)    (void) ;
extern void    (Load_Delete) (void*) ;
extern void    (Load_Scan)   (Load_t*,DataFile_t*) ;



#include "Region.h"

#define Load_MaxLengthOfRegionName      Region_MaxLengthOfRegionName
#define Load_MaxLengthOfKeyWord               (30)


#define Load_GetRegionTag(LOAD)          ((LOAD)->RegionTag)
#define Load_GetRegionName(LOAD)         ((LOAD)->RegionName)
#define Load_GetType(LOAD)               ((LOAD)->Type)
#define Load_GetNameOfEquation(LOAD)     ((LOAD)->NameOfEquation)
#define Load_GetFunctionIndex(LOAD)      ((LOAD)->FunctionIndex)
#define Load_GetFieldIndex(LOAD)         ((LOAD)->FieldIndex)
#define Load_GetFunction(LOAD)           ((LOAD)->Function)
#define Load_GetField(LOAD)              ((LOAD)->Field)
#define Load_GetFunctions(LOAD)          ((LOAD)->Functions)
#define Load_GetFields(LOAD)             ((LOAD)->Fields)


#define Load_TypeIs(LOAD,TYPE)           (!strcmp(Load_GetType(LOAD),TYPE))


#include "Functions.h"
#include "Fields.h"


struct Load_s {
  Function_t*  Function ;
  Field_t*     Field ;
  Functions_t* Functions ;
  Fields_t*    Fields ;
  int    RegionTag ;
  char*  RegionName ;
  char*  Type ;
  char*  NameOfEquation ;
  int    FunctionIndex ;
  int    FieldIndex ;
} ;


/* Old notations which I try to eliminate little by little */
#define char_t    Load_t

#endif
