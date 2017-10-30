#ifndef PREDEFINEDMETHODS_H
#define PREDEFINEDMETHODS_H


#include "Utils.h"

#define ComputeInitialState(...) \
        Utils_CAT_NARG(ComputeInitialState_,__VA_ARGS__)(__VA_ARGS__)
        
#define ComputeInitialState_1(A) \
        ComputeInitialState_2(A,double t)

#define ComputeInitialState_2 \
        ComputeInitialState_
        


#include "Models.h"

static Model_SetModelProp_t             SetModelProp ;
static Model_ReadMaterialProperties_t   ReadMatProp ;
static Model_PrintModelProp_t           PrintModelProp ;
static Model_DefineElementProperties_t  DefineElementProp ;
static Model_ComputeInitialState_t      ComputeInitialState_ ;
static Model_ComputeExplicitTerms_t     ComputeExplicitTerms ;
static Model_ComputeImplicitTerms_t     ComputeImplicitTerms ;
static Model_ComputeLoads_t             ComputeLoads ;
static Model_ComputeMatrix_t            ComputeMatrix ;
static Model_ComputeResidu_t            ComputeResidu ;
static Model_ComputeOutputs_t           ComputeOutputs ;

#include "BaseName_SetModelProp.h"

extern Model_SetModelProp_t BaseName_SetModelProp ;

int BaseName_SetModelProp(Model_t* model)
{
  Model_GetReadMaterialProperties(model) = ReadMatProp ;
  Model_GetPrintModelProp(model) = PrintModelProp ;
  Model_GetDefineElementProperties(model) = DefineElementProp ;
  Model_GetComputeInitialState(model) = ComputeInitialState_ ;
  Model_GetComputeExplicitTerms(model) = ComputeExplicitTerms ;
  Model_GetComputeImplicitTerms(model) = ComputeImplicitTerms ;
  Model_GetComputeLoads(model) = ComputeLoads ;
  Model_GetComputeMatrix(model) = ComputeMatrix ;
  Model_GetComputeResidu(model) = ComputeResidu ;
  Model_GetComputeOutputs(model) = ComputeOutputs ;
  
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
