#ifndef PREDEFINEDMODELMETHODS_H
#define PREDEFINEDMODELMETHODS_H


#include "Utils.h"

#define ComputeInitialState(...) \
        Utils_CAT_NARG(ComputeInitialState,__VA_ARGS__)(__VA_ARGS__)
        
#define ComputeInitialState1(A) \
        ComputeInitialState2(A,double t)
        

#include "IntFcts.h"
#include "ShapeFcts.h"

#define DefineElementProp(...) \
        Utils_CAT_NARG(DefineElementProp,__VA_ARGS__)(__VA_ARGS__)
        
#define DefineElementProp1(A) \
        DefineElementProp3(A,IntFcts_t* intfcts,ShapeFcts_t* shapefcts)
        
#define DefineElementProp2(A,B) \
        DefineElementProp3(A,B,ShapeFcts_t* shapefcts)


#define PrintModelProp(...) \
        Utils_CAT_NARG(PrintModelProp,__VA_ARGS__)(__VA_ARGS__)

#define PrintModelProp1(A) \
        PrintModelProp2(A,FILE* NULL)


#include "Model.h"

static Model_SetModelProperties_t       SetModelProp ;
static Model_ReadMaterialProperties_t   ReadMatProp ;
static Model_PrintModelProperties_t     PrintModelProp2 ;
static Model_DefineElementProperties_t  DefineElementProp3 ;
static Model_ComputeInitialState_t      ComputeInitialState2 ;
static Model_ComputeExplicitTerms_t     ComputeExplicitTerms ;
static Model_ComputeImplicitTerms_t     ComputeImplicitTerms ;
static Model_ComputeLoads_t             ComputeLoads ;
static Model_ComputeMatrix_t            ComputeMatrix ;
static Model_ComputeResidu_t            ComputeResidu ;
static Model_ComputeOutputs_t           ComputeOutputs ;
static Model_ComputePropertyIndex_t     ComputePropertyIndex ;

#include "BaseName.h"

#define BaseName_SetModelProp  BaseName(_SetModelProp)

extern Model_SetModelProperties_t BaseName_SetModelProp ;

int BaseName_SetModelProp(Model_t* model)
{
  Model_GetSetModelProperties(model) = &SetModelProp ;
  Model_GetReadMaterialProperties(model) = &ReadMatProp ;
  Model_GetPrintModelProperties(model) = &PrintModelProp2 ;
  Model_GetDefineElementProperties(model) = &DefineElementProp3 ;
  Model_GetComputeInitialState(model) = &ComputeInitialState2 ;
  Model_GetComputeExplicitTerms(model) = &ComputeExplicitTerms ;
  Model_GetComputeImplicitTerms(model) = &ComputeImplicitTerms ;
  Model_GetComputeLoads(model) = &ComputeLoads ;
  Model_GetComputeMatrix(model) = &ComputeMatrix ;
  Model_GetComputeResidu(model) = &ComputeResidu ;
  Model_GetComputeOutputs(model) = &ComputeOutputs ;
  
  /* ComputePropertyIndex must be defined in the model file */
  #if defined (ComputePropertyIndex) || defined (COMPUTEPROPERTYINDEX)
  Model_GetComputePropertyIndex(model) = &ComputePropertyIndex ;
  #endif
  
  #ifdef TITLE
  Model_CopyShortTitle(model,TITLE) ;
  #else
  Model_CopyShortTitle(model,"") ;
  #endif
  
  #ifdef AUTHORS
  Model_CopyNameOfAuthors(model,AUTHORS) ;
  #else
  Model_CopyNameOfAuthors(model,"Unknowns") ;
  #endif
  
  return(SetModelProp(model));
}


#define PrintModelChar    PrintModelProp

#endif
