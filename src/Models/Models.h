/*
** File:
** $Id: Models.h$
**
**
** Purpose:
** Unit specification for models.
**
*/
#ifndef MODELS_H
#define MODELS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* Vacuous declarations and typedef names */

/* class-like structure "Models_t" */
struct Models_s       ; typedef struct Models_s       Models_t ;


#include <stdio.h>


/* Declaration of Macros, Methods and Structures */

/* 1. Models_t */
#include "Model.h"
#include "Geometry.h"
#include "DataFile.h"


extern Models_t* (Models_New)(const int) ;
extern Models_t* (Models_Create)(DataFile_t*,Geometry_t*) ;
extern void      (Models_Delete)(void*) ;
extern void      (Models_Print)(char*,FILE *) ;
extern Model_t*  (Models_FindModel)(Models_t*,const char*) ;
extern int       (Models_FindModelIndex)(Models_t*,const char*) ;
extern Model_t*  (Models_FindOrAppendModel)(Models_t*,const char*,Geometry_t*,DataFile_t*) ;


#include "ListOfModels.h"

#define Models_NbOfModels               (ListOfModels_Nb)
#define Models_ListOfNames               ListOfModels_Names
#define Models_ListOfSetModelProp        ListOfModels_Methods(_SetModelProp)


#define Models_GetMaxNbOfModels(MODS) ((MODS)->maxn_model)
#define Models_GetNbOfModels(MODS)    ((MODS)->n_model)
#define Models_GetModel(MODS)         ((MODS)->model)



struct Models_s {             /* Models */
  unsigned int maxn_model ;   /* Maximun nb of models */
  unsigned int n_model ;      /* Nb of models */
  Model_t* model ;            /* Point to the first model */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
