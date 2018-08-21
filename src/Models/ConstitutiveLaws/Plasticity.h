#ifndef PLASTICITY_H
#define PLASTICITY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Plasticity_s     ; typedef struct Plasticity_s     Plasticity_t ;


/* 1. Plasticity_t */

extern Plasticity_t*  (Plasticity_Create)(void) ;
extern double         (Plasticity_UpdateElastoplasticTensor)(Plasticity_t*,double*) ;


#define Plasticity_MaxLengthOfKeyWord      (30)


/* Accessors */
#define Plasticity_GetCodeNameOfModel(PL)            ((PL)->codenameofmodel)
#define Plasticity_GetYieldFunctionGradient(PL)      ((PL)->dfsds)
#define Plasticity_GetPotentialFunctionGradient(PL)  ((PL)->dgsds)
#define Plasticity_GetHardeningModulus(PL)           ((PL)->hm)
#define Plasticity_GetCriterionValue(PL)             ((PL)->criterion)
#define Plasticity_GetFjiCijkl(PL)                   ((PL)->fc)
#define Plasticity_GetCijklGlk(PL)                   ((PL)->cg)
#define Plasticity_GetTangentStiffnessTensor(PL)     ((PL)->cijkl)
#define Plasticity_GetElasticity(PL)                 ((PL)->elasty)


#include "Elasticity.h"
#include "GenericData.h"

struct Plasticity_s {
  char*   codenameofmodel ;
  double* dfsds ;  /** Yield function gradient */
  double* dgsds ;  /** Potential function gradient */
  double* hm ;     /** Hardening modulus */
  double  criterion ;     /** Value of the yield function criterion */
  double* fc ;     /** fc(k,l) = dfsds(j,i) * C(i,j,k,l) */
  double* cg ;     /** cg(i,j) = C(i,j,k,l) * dgsds(l,k) */
  double  fcg ;    /** fcg = dfsds(j,i) * C(i,j,k,l) * dgsds(l,k)) */
  double* cijkl ;  /** Tangent stiffness tensor */
  GenericData_t* genericdata ;  /** Plastic properties */
  Elasticity_t* elasty ;
} ;

#endif
