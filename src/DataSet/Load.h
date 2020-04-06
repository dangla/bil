#ifndef LOAD_H
#define LOAD_H

/* vacuous declarations and typedef names */

/* class-like structure and attributes */
struct Load_s         ; typedef struct Load_s         Load_t ;



#include "DataFile.h"

extern Load_t* (Load_New)(void) ;
extern void    (Load_Scan)(Load_t*,DataFile_t*) ;



#define Load_MaxLengthOfKeyWord               (30)

#define Load_GetRegionIndex(LOAD)        ((LOAD)->reg)
#define Load_GetType(LOAD)               ((LOAD)->t)
#define Load_GetNameOfEquation(LOAD)     ((LOAD)->eqn)
#define Load_GetFunctionIndex(LOAD)      ((LOAD)->fctindex)
#define Load_GetFieldIndex(LOAD)         ((LOAD)->fldindex)
#define Load_GetFunction(LOAD)           ((LOAD)->fn)
#define Load_GetField(LOAD)              ((LOAD)->ch)
#define Load_GetFunctions(LOAD)          ((LOAD)->functions)
#define Load_GetFields(LOAD)             ((LOAD)->fields)


#define Load_TypeIs(LOAD,TYPE)           (!strcmp(Load_GetType(LOAD),TYPE))


#include "Functions.h"
#include "Fields.h"


struct Load_s {               /* chargement */
  int    reg ;                /* numero de la region */
  char*   t ;                 /* type de chargement */
  char*   eqn ;               /* nom de l'equation */
  int    fctindex ;           /* Time function index */
  int    fldindex ;           /* Field index */
  Function_t* fn ;            /* fonction du temps */
  Field_t* ch ;               /* champ */
  Functions_t* functions ;    /* Functions */
  Fields_t* fields ;          /* Fields */
} ;


/* Old notations which I try to eliminate little by little */
#define char_t    Load_t

#endif
