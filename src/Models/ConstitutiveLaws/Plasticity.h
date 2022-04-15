#ifndef PLASTICITY_H
#define PLASTICITY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Plasticity_s     ; typedef struct Plasticity_s     Plasticity_t ;


/* Typedef names of methods */
typedef double  (Plasticity_ComputeTangentStiffnessTensor_t)(Plasticity_t*,const double*,const double*,const double) ;
typedef double  (Plasticity_ComputeFunctionGradients_t)(Plasticity_t*,const double*,const double*) ;
typedef double  (Plasticity_Criterion_t)(Plasticity_t*,const double*,const double*) ;
typedef double  (Plasticity_ReturnMapping_t)(Plasticity_t*,double*,double*,double*) ;


/* 1. Plasticity_t */

extern Plasticity_t*  (Plasticity_Create)(void) ;
extern void           (Plasticity_Delete)(void*) ;
extern void           (Plasticity_Initialize)                  (Plasticity_t*) ;
extern void           (Plasticity_SetParameters)               (Plasticity_t*,...) ;
extern void           (Plasticity_SetParameter)                (Plasticity_t*,const char*,double) ;
extern double         (Plasticity_UpdateElastoplasticTensor)   (Plasticity_t*,double*) ;
extern void           (Plasticity_PrintTangentStiffnessTensor) (Plasticity_t*) ;
extern void           (Plasticity_CopyTangentStiffnessTensor)  (Plasticity_t*,double*) ;
//extern Plasticity_ComputeFunctionGradients_t     Plasticity_Criterion ;
//extern Plasticity_ReturnMapping_t                Plasticity_ReturnMapping ;



/* Accessors */
#define Plasticity_GetCodeNameOfModel(PL)            ((PL)->codenameofmodel)
#define Plasticity_GetYieldFunctionGradient(PL)      ((PL)->dfsds)
#define Plasticity_GetPotentialFunctionGradient(PL)  ((PL)->dgsds)
#define Plasticity_GetHardeningVariable(PL)          ((PL)->hardv)
#define Plasticity_GetHardeningModulus(PL)           ((PL)->hardm)
#define Plasticity_GetCriterionValue(PL)             ((PL)->criterion)
#define Plasticity_GetFjiCijkl(PL)                   ((PL)->fc)
#define Plasticity_GetCijklGlk(PL)                   ((PL)->cg)
#define Plasticity_GetFjiCijklGlk(PL)                ((PL)->fcg)
#define Plasticity_GetTangentStiffnessTensor(PL)     ((PL)->cijkl)
#define Plasticity_GetElasticity(PL)                 ((PL)->elasty)
#define Plasticity_GetParameter(PL)                  ((PL)->parameter)
#define Plasticity_GetComputeFunctionGradients(PL)   ((PL)->computefunctiongradients)
#define Plasticity_GetComputeTangentStiffnessTensor(PL)   ((PL)->computetangentstiffnesstensor)
#define Plasticity_GetReturnMapping(PL)              ((PL)->returnmapping)
#define Plasticity_GetPlasticMultiplier(PL)          ((PL)->lambda)


#define Plasticity_MaxLengthOfKeyWord             (100)
#define Plasticity_MaxNbOfParameters              (6)
#define Plasticity_MaxNbOfHardeningVariables      (2)
#define Plasticity_MaxNbOfPlasticMultiplier       (2) // Not used


#include "Elasticity.h"
#include "GenericData.h"
#include "Math_.h"
#include "Utils.h"




/* Drucker-Prager
 * -------------- */
#define Plasticity_IsDruckerPrager(PL) \
        Plasticity_Is(PL,"Drucker-Prager")

#define Plasticity_SetToDruckerPrager(PL) \
        do { \
          Plasticity_CopyCodeName(PL,"Drucker-Prager") ; \
          Plasticity_Initialize(PL) ; \
        } while(0)
        
#define Plasticity_GetFrictionAngle(PL) \
        Plasticity_GetParameter(PL)[0]

#define Plasticity_GetDilatancyAngle(PL) \
        Plasticity_GetParameter(PL)[1]
        
#define Plasticity_GetCohesion(PL) \
        Plasticity_GetParameter(PL)[2]
        
/* Cam-clay
 * -------- */
#define Plasticity_IsCamClay(PL) \
        Plasticity_Is(PL,"Cam-clay")

#define Plasticity_SetToCamClay(PL) \
        do { \
          Plasticity_CopyCodeName(PL,"Cam-clay") ; \
          Plasticity_Initialize(PL) ; \
        } while(0)
        
#define Plasticity_GetSlopeSwellingLine(PL) \
        Plasticity_GetParameter(PL)[0]

#define Plasticity_GetSlopeVirginConsolidationLine(PL) \
        Plasticity_GetParameter(PL)[1]
        
#define Plasticity_GetSlopeCriticalStateLine(PL) \
        Plasticity_GetParameter(PL)[2]
        
#define Plasticity_GetInitialPreconsolidationPressure(PL) \
        Plasticity_GetParameter(PL)[3]
        
#define Plasticity_GetInitialVoidRatio(PL) \
        Plasticity_GetParameter(PL)[4]
        
/* Cam-clayEp
 * ---------- */
#define Plasticity_IsCamClayEp(PL) \
        Plasticity_Is(PL,"Cam-clayEp")

#define Plasticity_SetToCamClayEp(PL) \
        do { \
          Plasticity_CopyCodeName(PL,"Cam-clayEp") ; \
          Plasticity_Initialize(PL) ; \
        } while(0)
        
/* Cam-clayOffset
 * -------------- */
#define Plasticity_IsCamClayOffset(PL) \
        Plasticity_Is(PL,"Cam-clayOffset")

#define Plasticity_SetToCamClayOffset(PL) \
        do { \
          Plasticity_CopyCodeName(PL,"Cam-clayOffset") ; \
          Plasticity_Initialize(PL) ; \
        } while(0)
        
        
        
/* Shorthands */
#define Plasticity_ComputeTangentStiffnessTensor(...) \
        Utils_CAT_NARG(Plasticity_ComputeTangentStiffnessTensor,__VA_ARGS__)(__VA_ARGS__)
        
#define Plasticity_ComputeTangentStiffnessTensor3(...) \
        Plasticity_ComputeTangentStiffnessTensor4(__VA_ARGS__,0)
        
#define Plasticity_ComputeTangentStiffnessTensor4(PL,...) \
        Plasticity_GetComputeTangentStiffnessTensor(PL)(PL,__VA_ARGS__)
        
#define Plasticity_ComputeFunctionGradients \
        Plasticity_ComputeTangentStiffnessTensor3
        
/*
#define Plasticity_ComputeFunctionGradients(PL,...) \
        Plasticity_GetComputeFunctionGradients(PL)(PL,__VA_ARGS__)
*/


#define Plasticity_ReturnMapping(PL,...) \
        Plasticity_GetReturnMapping(PL)(PL,__VA_ARGS__)

#define Plasticity_ComputeElasticTensor(PL,...) \
        Elasticity_ComputeStiffnessTensor(Plasticity_GetElasticity(PL),__VA_ARGS__)
        
#define Plasticity_CopyElasticTensor(PL,...) \
        Elasticity_CopyStiffnessTensor(Plasticity_GetElasticity(PL),__VA_ARGS__)
        


/* Implementation */
#define Plasticity_CopyCodeName(PL,TYP) \
        memcpy(Plasticity_GetCodeNameOfModel(PL),TYP,MAX(strlen(TYP),Plasticity_MaxLengthOfKeyWord))

#define Plasticity_Is(PL,TYP) \
        (!strcmp(Plasticity_GetCodeNameOfModel(PL),TYP))




struct Plasticity_s {
  char*   codenameofmodel ;
  double* dfsds ;  /** Yield function gradient */
  double* dgsds ;  /** Potential function gradient */
  double* hardv ;  /** Hardening variable */
  double* hardm ;  /** Hardening modulus */
  double  criterion ;     /** Value of the yield function criterion */
  double* fc ;     /** fc(k,l) = dfsds(j,i) * C(i,j,k,l) */
  double* cg ;     /** cg(i,j) = C(i,j,k,l) * dgsds(l,k) */
  double  fcg ;    /** fcg = hardm + dfsds(j,i) * C(i,j,k,l) * dgsds(l,k)) */
  double* cijkl ;  /** Tangent stiffness tensor */
  double  lambda ; /** Plastic multiplier */
  double* parameter ;
  GenericData_t* genericdata ;  /** Plastic properties */
  Elasticity_t* elasty ;
  //Plasticity_ComputeFunctionGradients_t*  computefunctiongradients ;
  Plasticity_ComputeTangentStiffnessTensor_t* computetangentstiffnesstensor ;
  Plasticity_ReturnMapping_t*             returnmapping ;
} ;

#endif
