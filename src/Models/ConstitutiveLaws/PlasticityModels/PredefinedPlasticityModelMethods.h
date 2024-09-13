#ifndef PREDEFINEDPLASTICITYMODELMETHODS_H
#define PREDEFINEDPLASTICITYMODELMETHODS_H

#include "autodiff"

template<typename T>
T* (YieldFunction)(Plasticity_t*,const T*,const T*);

template<typename T>
T* (FlowRules)(Plasticity_t*,const T*,const T*);


static Plasticity_SetModelProp_t                     SetModelProp ;
static Plasticity_ComputeTangentStiffnessTensor_t    ComputeTangentStiffnessTensor ;
static Plasticity_ReturnMapping_t                    ReturnMapping ;

static Plasticity_YieldFunction_t                    YieldFunction ;
static Plasticity_FlowRules_t                        FlowRules ;
#ifdef HAVE_AUTODIFF
static Plasticity_YieldFunctionDual_t                YieldFunction ;
static Plasticity_FlowRulesDual_t                    FlowRules ;
#endif

static Plasticity_SetParameters_t                    SetParameters ;


#include "BaseName.h"

#define BaseName_SetPlasticityModelProp  BaseName(_SetPlasticityModelProp)

extern Plasticity_SetModelProp_t BaseName_SetPlasticityModelProp ;

int BaseName_SetPlasticityModelProp(Plasticity_t* plasty)
{
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = ComputeTangentStiffnessTensor ;
    Plasticity_GetReturnMapping(plasty)                 = ReturnMapping ;
    Plasticity_GetYieldFunction(plasty)                 = YieldFunction<double> ;
    Plasticity_GetFlowRules(plasty)                     = FlowRules<double> ;
    #ifdef HAVE_AUTODIFF
    Plasticity_GetYieldFunctionDual(plasty)             = YieldFunction<real> ;
    Plasticity_GetFlowRulesDual(plasty)                 = FlowRules<real> ;
    #endif
    Plasticity_GetSetParameters(plasty)                 = SetParameters ;
  
  return(SetModelProp(plasty));
}

#endif
