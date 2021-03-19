#ifndef DAMAGE_H
#define DAMAGE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Damage_s     ; typedef struct Damage_s     Damage_t ;


/* Typedef names of methods */
typedef double  (Damage_ComputeFunctionGradients_t)(Damage_t*,const double*,const double*,const double*) ;
typedef double  (Damage_Criterion_t)(Damage_t*,const double*,const double*) ;
typedef double  (Damage_ReturnMapping_t)(Damage_t*,double*,double*,double*) ;


/* 1. Damage_t */

extern Damage_t*      (Damage_Create)(void) ;
extern void           (Damage_Delete)(void*) ;
extern void           (Damage_Initialize)                  (Damage_t*) ;
extern void           (Damage_SetParameters)               (Damage_t*,...) ;
extern void           (Damage_SetParameter)                (Damage_t*,const char*,double) ;
extern double         (Damage_UpdateTangentStiffnessTensor)(Damage_t*,double*) ;
extern void           (Damage_PrintStiffnessTensor)        (Damage_t*) ;
extern void           (Damage_CopyStiffnessTensor)         (Damage_t*,double*) ;
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
#define Damage_GetTangentStiffnessTensor(D)     ((D)->cijkl)
#define Damage_GetStiffnessTensor(D)            ((D)->cijkl)
#define Damage_GetElasticity(D)                 ((D)->elasty)
#define Damage_GetParameter(D)                  ((D)->parameter)
#define Damage_GetComputeFunctionGradients(D)   ((D)->computefunctiongradients)
#define Damage_GetReturnMapping(D)              ((D)->returnmapping)


#define Damage_MaxLengthOfKeyWord             (100)
#define Damage_MaxNbOfParameters              (6)
#define Damage_MaxNbOfHardeningVariables      (2)


#include "Elasticity.h"
#include "GenericData.h"
#include "Tools/Math.h"





/* Mazars
 * ------ */
#define Damage_IsMazars(D) \
        Damage_Is(D,"Mazars")

#define Damage_SetToMazars(D) \
        do { \
          Damage_CopyCodeName(D,"Mazars") ; \
          Damage_Initialize(D) ; \
        } while(0)
        
#define Damage_GetStrainAtUniaxialTensileStrength(D) \
        Damage_GetParameter(D)[0]
        
#define Damage_GetA_t(D) \
        Damage_GetParameter(D)[1]
        
#define Damage_GetB_t(D) \
        Damage_GetParameter(D)[2]
        
#define Damage_GetA_c(D) \
        Damage_GetParameter(D)[3]
        
#define Damage_GetB_c(D) \
        Damage_GetParameter(D)[4]





/* Marigo-Jirasek
 * -------------- */
#define Damage_IsMarigoJirasek(D) \
        Damage_Is(D,"MarigoJirasek")

#define Damage_SetToMarigoJirasek(D) \
        do { \
          Damage_CopyCodeName(D,"MarigoJirasek") ; \
          Damage_Initialize(D) ; \
        } while(0)
        
#define Damage_GetCriticalEnergyReleaseRate(D) \
        Damage_GetParameter(D)[0]
        
#define Damage_GetMaximumEnergyReleaseRate(D) \
        Damage_GetParameter(D)[1]
        
#define Damage_GetUniaxialTensileStrength(D) \
        Damage_GetParameter(D)[2]
        
#define Damage_GetFractureEnergy(D) \
        Damage_GetParameter(D)[3]
        
#define Damage_GetCrackBandWidth(D) \
        Damage_GetParameter(D)[4]






/* Basic
 * ----- */
#define Damage_IsBasic(D) \
        Damage_Is(D,"Basic")

#define Damage_SetToBasic(D) \
        do { \
          Damage_CopyCodeName(D,"Basic") ; \
          Damage_Initialize(D) ; \
        } while(0)
        
#define Damage_GetEnergyThreshold(D) \
        Damage_GetParameter(D)[0]
        
#define Damage_GetA_t(D) \
        Damage_GetParameter(D)[1]
        
#define Damage_GetB_t(D) \
        Damage_GetParameter(D)[2]
        
#define Damage_GetA_c(D) \
        Damage_GetParameter(D)[3]
        
#define Damage_GetB_c(D) \
        Damage_GetParameter(D)[4]
        
        
        
/* Shorthands */
#define Damage_ComputeFunctionGradients(D,...) \
        Damage_GetComputeFunctionGradients(D)(D,__VA_ARGS__)
        
#define Damage_ReturnMapping(D,...) \
        Damage_GetReturnMapping(D)(D,__VA_ARGS__)

#define Damage_ComputeElasticTensor(D,...) \
        Elasticity_ComputeStiffnessTensor(Damage_GetElasticity(D),__VA_ARGS__)
        
#define Damage_CopyElasticTensor(D,...) \
        Elasticity_CopyStiffnessTensor(Damage_GetElasticity(D),__VA_ARGS__)
        


/* Implementation */
#define Damage_CopyCodeName(D,TYP) \
        memcpy(Damage_GetCodeNameOfModel(D),TYP,MAX(strlen(TYP),Damage_MaxLengthOfKeyWord))

#define Damage_Is(D,TYP) \
        (!strcmp(Damage_GetCodeNameOfModel(D),TYP))




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
  double* cijkl ;  /** Stiffness tensor of the damaged material */
  double* parameter ;
  GenericData_t* genericdata ;  /** Damage properties */
  Elasticity_t* elasty ;
  Damage_ComputeFunctionGradients_t*  computefunctiongradients ;
  Damage_ReturnMapping_t*             returnmapping ;
} ;

#endif
