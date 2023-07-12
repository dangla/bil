#ifndef DAMAGE_H
#define DAMAGE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Damage_s     ; typedef struct Damage_s     Damage_t ;


/* Typedef names of methods */
typedef double  (Damage_ComputeTangentStiffnessTensor_t)(Damage_t*,const double*,const double*,const double*) ;
typedef double  (Damage_ReturnMapping_t)(Damage_t*,double*,double*,double*) ;
typedef void    (Damage_SetParameters_t)(Damage_t*,...) ;
typedef void    (Damage_SetModelProp_t) (Damage_t*) ;


/* 1. Damage_t */

extern Damage_t*      (Damage_Create)(void) ;
extern void           (Damage_Delete)(void*) ;
extern void           (Damage_Initialize)                  (Damage_t*) ;
//extern void           (Damage_SetParameters)               (Damage_t*,...) ;
extern double         (Damage_UpdateTangentStiffnessTensor)(Damage_t*,double*) ;
extern void           (Damage_PrintTangentStiffnessTensor) (Damage_t*) ;
extern void           (Damage_CopyTangentStiffnessTensor)  (Damage_t*,double*) ;
extern void           (Damage_CopyDamagedStiffnessTensor)  (Damage_t*,double*) ;
//extern Damage_ComputeFunctionGradients_t     Damage_Criterion ;
//extern Damage_ReturnMapping_t                Damage_ReturnMapping ;



/* Accessors */
#define Damage_GetCodeNameOfModel(D)            ((D)->codenameofmodel)
#define Damage_GetYieldFunctionGradient(D)      ((D)->dfsds)
#define Damage_GetPotentialFunctionGradient(D)  ((D)->dgsds)
#define Damage_GetHardeningVariable(D)          ((D)->hardv)
#define Damage_GetHardeningModulus(D)           ((D)->hardm)
#define Damage_GetCriterionValue(D)             ((D)->criterion)
#define Damage_GetFjiCijkl(D)                   ((D)->fc)
#define Damage_GetCijklGlk(D)                   ((D)->cg)
#define Damage_GetFjiCijklGlk(D)                ((D)->fcg)
#define Damage_GetDamagedStiffnessTensor(D)     ((D)->cdamaged)
#define Damage_GetTangentStiffnessTensor(D)     ((D)->ctangent)
#define Damage_GetSetParameters(D)              ((D)->setparameters)
#define Damage_GetSetModelProp(D)               ((D)->setmodelprop)
#define Damage_GetElasticity(D)                 ((D)->elasty)
#define Damage_GetParameter(D)                  ((D)->parameter)
#define Damage_GetComputeTangentStiffnessTensor(D)   ((D)->computetangentstiffnesstensor)
#define Damage_GetReturnMapping(D)              ((D)->returnmapping)


#define Damage_MaxLengthOfKeyWord             (100)
#define Damage_MaxNbOfParameters              (6)
#define Damage_MaxNbOfHardeningVariables      (2)



#include "Arg.h"

#define Damage_SetParameters(D,...) \
        do {\
          if(Damage_MaxNbOfParameters < Arg_NARG(__VA_ARGS__)) {\
            Message_RuntimeError("Damage_SetParameters: too many parameters") ;\
          }\
          Damage_GetSetParameters(D)(D,__VA_ARGS__) ;\
        } while(0)


#include "Elasticity.h"
#include "GenericData.h"
#include "Math_.h"



        
/* Shorthands */
#define Damage_ComputeTangentStiffnessTensor(D,...) \
        Damage_GetComputeTangentStiffnessTensor(D)(D,__VA_ARGS__)
        
#define Damage_ReturnMapping(D,...) \
        Damage_GetReturnMapping(D)(D,__VA_ARGS__)

#define Damage_ComputeElasticTensor(D,...) \
        Elasticity_ComputeStiffnessTensor(Damage_GetElasticity(D),__VA_ARGS__)
        
#define Damage_CopyElasticTensor(D,...) \
        Elasticity_CopyStiffnessTensor(Damage_GetElasticity(D),__VA_ARGS__)


#define Damage_ComputeFunctionGradients \
        Damage_ComputeTangentStiffnessTensor
        

#define Damage_SetModelProp(D) \
        Damage_GetSetModelProp(D)(D)
        


/* Implementation */
#define Damage_CopyCodeName(D,STR) \
        memcpy(Damage_GetCodeNameOfModel(D),STR,MAX(strlen(STR),Damage_MaxLengthOfKeyWord))

#define Damage_Is(D,STR) \
        (!strcmp(Damage_GetCodeNameOfModel(D),STR))


#define Damage_SetTo(D,MOD) \
        do { \
          Damage_CopyCodeName(D,Utils_STR(MOD)) ; \
          Damage_Initialize(D) ; \
          Damage_SetModelProp(D) ; \
        } while(0)




struct Damage_s {
  char*   codenameofmodel ;
  double* dfsds ;  /** Yield function gradient */
  double* dgsds ;  /** Potential function gradient */
  double* hardv ;  /** Hardening variable */
  double* hardm ;  /** Hardening modulus */
  double  criterion ;     /** Value of the yield function criterion */
  double* fc ;     /** fc(k,l) = dfsds(j,i) * C(i,j,k,l) */
  double* cg ;     /** cg(i,j) = C(i,j,k,l) * dgsds(l,k) */
  double  fcg ;    /** fcg = hardm + dfsds(j,i) * C(i,j,k,l) * dgsds(l,k)) */
  double* cdamaged ;  /** Damaged stiffness tensor */
  double* ctangent ;  /** Tangent stiffness tensor */
  double* parameter ;
  GenericData_t* genericdata ;  /** Damage properties */
  Elasticity_t* elasty ;
  Damage_ComputeTangentStiffnessTensor_t*  computetangentstiffnesstensor ;
  Damage_ReturnMapping_t*             returnmapping ;
  Damage_SetParameters_t*             setparameters ;
  Damage_SetModelProp_t*              setmodelprop ;
} ;

#endif
