/*
** File:
** $Id: Model.h$
**
**
** Purpose:
** Unit specification for model.
**
*/
#ifndef MODEL_H
#define MODEL_H


/* Vacuous declarations and typedef names */

/* class-like structure "Model_t" */
struct Model_s        ; typedef struct Model_s        Model_t ;


/*  Typedef names of Methods */
#include <stdio.h>

typedef int    (Model_SetModelProp_t)        (Model_t*) ;
typedef int    (Model_ComputePropertyIndex_t)(const char*) ;
typedef int    (Model_PrintModelProp_t)      (Model_t*,FILE*) ;

#include "Element.h"

typedef int    (Model_ComputeInitialState_t) (Element_t*,double) ;
typedef int    (Model_ComputeExplicitTerms_t)(Element_t*,double) ;
typedef int    (Model_ComputeImplicitTerms_t)(Element_t*,double,double) ;
typedef int    (Model_ComputeMatrix_t)       (Element_t*,double,double,double*) ;
typedef int    (Model_ComputeResidu_t)       (Element_t*,double,double,double*) ;

typedef void*  (Model_ComputeVariables_t)(Element_t*,void*,void*,void*,const double,const double,const int) ;
//typedef double*   (Model_ComputeVariableFluxes_t)(Element_t* el,double** u,double t,double dt,int i,int j,...) ;
typedef void   (Model_ComputeSecondaryVariables_t)(Element_t*,const double,const double,double*) ;
typedef void*  (Model_ComputeVariableDerivatives_t)(Element_t*,const double,const double,void*,const double,const int) ;

#include "IntFcts.h"
#include "ShapeFcts.h"

typedef int    (Model_DefineElementProperties_t)(Element_t*,IntFcts_t*,ShapeFcts_t*) ;

#include "Load.h"

typedef int    (Model_ComputeLoads_t)(Element_t*,double,double,Load_t*,double*) ;

#include "Result.h"

typedef int    (Model_ComputeOutputs_t)(Element_t*,double,double*,Result_t*) ;

#include "Material.h"
#include "DataFile.h"

typedef int    (Model_ReadMaterialProperties_t)(Material_t*,DataFile_t*) ;



/* 2. Model_t */
#include "Geometry.h"
#include "DataFile.h"

extern Model_t*  (Model_New)       (void) ;
extern void      (Model_Delete)    (void*) ;
extern Model_t*  (Model_Initialize)(Model_t*,const char*,Geometry_t*,DataFile_t*) ;
extern double*   (Model_ComputeVariableDerivatives)(Element_t*,double,double,double,int,int) ;
extern void      (Model_Scan)(Model_t*,DataFile_t*,Geometry_t*) ;


#include "Views.h"

#define Model_MaxLengthOfKeyWord        (30)
#define Model_MaxNbOfEquations          (10)
#define Model_MaxLengthOfShortTitle     (80)
#define Model_MaxLengthOfAuthorNames    (80)
#define Model_MaxNbOfViews              (Views_MaxNbOfViews)

#define Model_MaxNbOfVariables          (100)
#define Model_MaxNbOfVariableFluxes     (Model_MaxNbOfVariables)


/* Accessors */
#define Model_GetCodeNameOfModel(MOD)      ((MOD)->codename)
#define Model_GetNbOfEquations(MOD)        ((MOD)->nbofequations)
#define Model_GetNameOfEquation(MOD)       ((MOD)->nameofequations)
#define Model_GetNameOfUnknown(MOD)        ((MOD)->nameofunknowns)
#define Model_GetSequentialIndexOfUnknown(MOD)        ((MOD)->sequentialindex)
#define Model_GetGeometry(MOD)             ((MOD)->geometry)
#define Model_GetDataFile(MOD)             ((MOD)->datafile)
#define Model_GetShortTitle(MOD)           ((MOD)->shorttitle)
#define Model_GetNameOfAuthors(MOD)        ((MOD)->authors)
#define Model_GetNumericalMethod(MOD)      ((MOD)->numericalmethod)
#define Model_GetObjectiveValue(MOD)       ((MOD)->obval)
#define Model_GetViews(MOD)                ((MOD)->views)
#define Model_GetLocalVariableVectors(MOD) ((MOD)->localvariable)
#define Model_GetLocalFluxVectors(MOD)     ((MOD)->localflux)
//#define Model_GetNbOfVariables(MOD)        ((MOD)->nbofvariables)
//#define Model_GetNbOfVariableFluxes(MOD)   ((MOD)->nbofvariablefluxes)


#define Model_GetSetModelProp(MOD)            ((MOD)->setmodelprop)
#define Model_GetReadMaterialProperties(MOD)  ((MOD)->readmatprop)
#define Model_GetPrintModelProp(MOD)          ((MOD)->printmodelprop)
#define Model_GetDefineElementProperties(MOD) ((MOD)->defineelementprop)
#define Model_GetComputeInitialState(MOD)     ((MOD)->computeinitialstate)
#define Model_GetComputeExplicitTerms(MOD)    ((MOD)->computeexplicitterms)
#define Model_GetComputeImplicitTerms(MOD)    ((MOD)->computeimplicitterms)
#define Model_GetComputeMatrix(MOD)           ((MOD)->computematrix)
#define Model_GetComputeResidu(MOD)           ((MOD)->computeresidu)
#define Model_GetComputeLoads(MOD)            ((MOD)->computeloads)
#define Model_GetComputeOutputs(MOD)          ((MOD)->computeoutputs)
#define Model_GetComputePropertyIndex(MOD)    ((MOD)->computepropertyindex)

#define Model_GetComputeSecondaryVariables(MOD)  ((MOD)->computesecondaryvariables)



/* Copy operations */
#define Model_CopyNameOfEquation(MOD,index,name) \
        (strcpy(Model_GetNameOfEquation(MOD)[index],name))

#define Model_CopyNameOfUnknown(MOD,index,name) \
        (strcpy(Model_GetNameOfUnknown(MOD)[index],name))

#define Model_CopyCodeNameOfModel(MOD,codename) \
        (strcpy(Model_GetCodeNameOfModel(MOD),codename))

#define Model_CopyShortTitle(MOD,title) \
        (strcpy(Model_GetShortTitle(MOD),title))

#define Model_CopyNameOfAuthors(MOD,authors) \
        (strcpy(Model_GetNameOfAuthors(MOD),authors))



/* Dimension */
#define Model_GetDimension(MOD) \
        Geometry_GetDimension(Model_GetGeometry(MOD))


/* Short hands */
#define Model_SetModelProp(MOD) \
        Model_GetSetModelProp(MOD)(MOD)

#define Model_PrintModelProp(MOD,file) \
        Model_GetPrintModelProp(MOD)(MOD,file)



#include "LocalVariableVectors.h"

/* Operations on local variables */
#define Model_GetNbOfVariables(MOD) \
        LocalVariableVectors_GetNbOfVariables(Model_GetLocalVariableVectors(MOD))

#define Model_GetVariable(MOD,n) \
        LocalVariableVectors_GetVariable(Model_GetLocalVariableVectors(MOD),n)

#define Model_GetVariableDerivative(MOD,n) \
        LocalVariableVectors_GetVariableDerivative(Model_GetLocalVariableVectors(MOD),n)


#define Model_GetNbOfVariableFluxes(MOD) \
        LocalVariableVectors_GetNbOfVariables(Model_GetLocalFluxVectors(MOD))

#define Model_GetVariableFlux(MOD,n) \
        LocalVariableVectors_GetVariable(Model_GetLocalFluxVectors(MOD),n)

#define Model_GetVariableFluxDerivative(MOD,n) \
        LocalVariableVectors_GetVariableDerivative(Model_GetLocalFluxVectors(MOD),n)









#include "ObVal.h"


struct Model_s {              /* model */
  Model_SetModelProp_t*             setmodelprop ;
  Model_ReadMaterialProperties_t*   readmatprop ;
  Model_PrintModelProp_t*           printmodelprop ;
  Model_DefineElementProperties_t*  defineelementprop ;
  Model_ComputeInitialState_t*      computeinitialstate ;
  Model_ComputeExplicitTerms_t*     computeexplicitterms ;
  Model_ComputeImplicitTerms_t*     computeimplicitterms ;
  Model_ComputeMatrix_t*            computematrix ;
  Model_ComputeResidu_t*            computeresidu ;
  Model_ComputeLoads_t*             computeloads ;
  Model_ComputeOutputs_t*           computeoutputs ;
  Model_ComputePropertyIndex_t*     computepropertyindex ;
  
  Model_ComputeSecondaryVariables_t* computesecondaryvariables ;
  
  char*   codename ;          /* code name of the model */
  char*   shorttitle ;        /* Short title of the model */
  char*   authors ;           /* Authors of the model */
  
  unsigned short nbofequations ;    /* Number of equations */
  char**   nameofequations ;  /* Names of equations */
  char**   nameofunknowns ;   /* Names of unknowns */
  int*     sequentialindex ;  /* Sequential indexes of unknowns/equations */
  
  Geometry_t* geometry  ;     /* Geometry of the problem being dealt with */
  DataFile_t* datafile ;      /* Datafile being dealt with */
  
  void*    numericalmethod ;  /* Numerical method */
  ObVal_t* obval ;            /* Objective values of unknowns */
  Views_t* views ;            /* Views */
  
  LocalVariableVectors_t* localvariable ;
  LocalVariableVectors_t* localflux ;
  
  //unsigned int nbofvariables ;
  //unsigned int nbofvariablefluxes ;
} ;


/* Old notations which should be eliminated */
#define MAX_EQUATIONS     Model_MaxNbOfEquations

#endif
