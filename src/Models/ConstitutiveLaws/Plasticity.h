#ifndef PLASTICITY_H
#define PLASTICITY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Plasticity_s     ; typedef struct Plasticity_s     Plasticity_t ;


/* Typedef names of methods */
typedef double* (Plasticity_ComputeTangentStiffnessTensor_t)(Plasticity_t*,const double*,const double*,const double*) ;
typedef double* (Plasticity_ReturnMapping_t)(Plasticity_t*,double*,double*,double*) ;
typedef double* (Plasticity_YieldFunction_t)(Plasticity_t*,const double*,const double*) ;
typedef double* (Plasticity_FlowRules_t)    (Plasticity_t*,const double*,const double*) ;
typedef void    (Plasticity_SetParameters_t)(Plasticity_t*,...) ;
typedef void    (Plasticity_SetModelProp_t) (Plasticity_t*) ;


/* 1. Plasticity_t */

extern Plasticity_t*  (Plasticity_Create)(void) ;
extern void           (Plasticity_Delete)(void*) ;
extern void           (Plasticity_Initialize)                  (Plasticity_t*) ;
//extern void           (Plasticity_SetParameters)               (Plasticity_t*,...) ;
//extern void           (Plasticity_SetParameter)                (Plasticity_t*,const char*,double) ;
extern double         (Plasticity_UpdateElastoplasticTensor)   (Plasticity_t*,double*) ;
extern void           (Plasticity_PrintTangentStiffnessTensor) (Plasticity_t*) ;
extern void           (Plasticity_CopyTangentStiffnessTensor)  (Plasticity_t*,double*) ;
extern Plasticity_ComputeTangentStiffnessTensor_t (Plasticity_GenericTangentStiffnessTensor) ;
extern Plasticity_ComputeTangentStiffnessTensor_t (Plasticity_GenericTangentStiffnessTensor1) ;
extern Plasticity_ReturnMapping_t                 (Plasticity_GenericReturnMapping) ;



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
#define Plasticity_GetComputeTangentStiffnessTensor(PL)   ((PL)->computetangentstiffnesstensor)
#define Plasticity_GetReturnMapping(PL)              ((PL)->returnmapping)
#define Plasticity_GetYieldFunction(PL)              ((PL)->yieldfunction)
#define Plasticity_GetFlowRules(PL)                  ((PL)->flowrules)
#define Plasticity_GetSetParameters(PL)              ((PL)->setparameters)
#define Plasticity_GetSetModelProp(PL)               ((PL)->setmodelprop)
#define Plasticity_GetPlasticMultiplier(PL)          ((PL)->lambda)
#define Plasticity_GetBuffers(PL)                    ((PL)->buffers)
#define Plasticity_GetNbOfCriteria(PL)               ((PL)->ncriteria)
#define Plasticity_GetNbOfHardeningVariables(PL)     ((PL)->nhardv)
#define Plasticity_GetNbOfNonHardeningVariables(PL)  ((PL)->nnhardv)
#define Plasticity_GetTypicalSmallIncrementOfHardeningVariable(PL)   ((PL)->dhardv)
#define Plasticity_GetTypicalSmallIncrementOfStress(PL)              ((PL)->dstress)
#define Plasticity_GetGenericData(PL)                ((PL)->genericdata)
#define Plasticity_GetCurves(PL)                     ((PL)->curves)


/* Buffer */
#define Plasticity_GetBuffer(PL) \
        Buffers_GetBufferOfCurrentThread(Plasticity_GetBuffers(PL))


#define Plasticity_MaxLengthOfKeyWord             (100)
#define Plasticity_MaxNbOfParameters              (20)
#define Plasticity_MaxNbOfHardeningVariables      (5)
#define Plasticity_MaxNbOfCriteria                (2)
#define Plasticity_MaxNbOfCurves                  (4)

#define Plasticity_SizeOfBuffer \
        (1000*sizeof(double))



#include "Utils.h"
#include "Arg.h"

#define Plasticity_SetParameters(PL,...) \
        do {\
          if(Plasticity_MaxNbOfParameters < Arg_NARG(__VA_ARGS__)) {\
            Message_RuntimeError("Plasticity_SetParameters: too many parameters") ;\
          }\
          Plasticity_GetSetParameters(PL)(PL,__VA_ARGS__) ;\
        } while(0)
        
        
        
/* Shorthands */
#define Plasticity_ComputeTangentStiffnessTensor(...) \
        Utils_CAT_NARG(Plasticity_ComputeTangentStiffnessTensor,__VA_ARGS__)(__VA_ARGS__)
        
#define Plasticity_ComputeTangentStiffnessTensor3(...) \
        Plasticity_ComputeTangentStiffnessTensor4(__VA_ARGS__,NULL)

        
/* We use a C extension provided by GNU C:
 * A compound statement enclosed in parentheses may appear 
 * as an expression in GNU C.
 * (https://gcc.gnu.org/onlinedocs/gcc/Statement-Exprs.html#Statement-Exprs) */
#if 1
#define Plasticity_ComputeTangentStiffnessTensor4(PL,...) \
        ({ \
          Plasticity_ComputeTangentStiffnessTensor_t* Plasticity_ct =  Plasticity_GetComputeTangentStiffnessTensor(PL); \
          (Plasticity_ct) ? Plasticity_ct(PL,__VA_ARGS__) : \
          Plasticity_GenericTangentStiffnessTensor(PL,__VA_ARGS__); \
        })
#endif
        
#define Plasticity_ComputeFunctionGradients \
        Plasticity_ComputeTangentStiffnessTensor3


#if 1
#define Plasticity_ReturnMapping(PL,...) \
        ({ \
          Plasticity_ReturnMapping_t* Plasticity_rm = Plasticity_GetReturnMapping(PL); \
          (Plasticity_rm) ? Plasticity_rm(PL,__VA_ARGS__) : \
          Plasticity_GenericReturnMapping(PL,__VA_ARGS__) ; \
        })
#endif


#define Plasticity_ComputeElasticTensor(PL,...) \
        Elasticity_ComputeStiffnessTensor(Plasticity_GetElasticity(PL),__VA_ARGS__)
        
#define Plasticity_CopyElasticStiffnessTensor(PL,...) \
        Elasticity_CopyStiffnessTensor(Plasticity_GetElasticity(PL),__VA_ARGS__)
        
#define Plasticity_CopyElasticTensor \
        Plasticity_CopyElasticStiffnessTensor
        
#define Plasticity_CopyElasticComplianceTensor(PL,...) \
        Elasticity_CopyComplianceTensor(Plasticity_GetElasticity(PL),__VA_ARGS__)


#define Plasticity_YieldFunction(PL,...) \
        Plasticity_GetYieldFunction(PL)(PL,__VA_ARGS__)

#define Plasticity_FlowRules(PL,...) \
        Plasticity_GetFlowRules(PL)(PL,__VA_ARGS__)

#define Plasticity_SetModelProp(PL) \
        Plasticity_GetSetModelProp(PL)(PL)
        


/* Implementation */
#define Plasticity_CopyCodeName(PL,STR) \
        memcpy(Plasticity_GetCodeNameOfModel(PL),STR,MAX(strlen(STR),Plasticity_MaxLengthOfKeyWord))

#define Plasticity_Is(PL,STR) \
        (!strcmp(Plasticity_GetCodeNameOfModel(PL),STR))

#define Plasticity_SetTo(PL,MOD) \
        do { \
          Plasticity_CopyCodeName(PL,Utils_STR(MOD)) ; \
          Plasticity_Initialize(PL) ; \
          Plasticity_SetModelProp(PL) ; \
        } while(0)
          


/* Operations on buffer */
#define Plasticity_AllocateInBuffer(PL,sz) \
        Buffer_Allocate(Plasticity_GetBuffer(PL),(sz))
        
#define Plasticity_FreeBuffer(PL)  \
        Buffer_Free(Plasticity_GetBuffer(PL))
        
#define Plasticity_FreeBufferFrom(PL,p) \
        Buffer_FreeFrom(Plasticity_GetBuffer(PL),(char*) (p))



/* GenericData */
#define Plasticity_AppendGenericData(PL,GD) \
        do { \
          if(Plasticity_GetGenericData(PL)) { \
            GenericData_Append(Plasticity_GetGenericData(PL),GD) ; \
          } else { \
            Plasticity_GetGenericData(PL) = GD ; \
          } \
        } while(0)
        
#define Plasticity_FindGenericData(PL,...) \
        GenericData_Find(Plasticity_GetGenericData(PL),__VA_ARGS__)
        
#define Plasticity_FindData(PL,...) \
        GenericData_FindData(Plasticity_GetGenericData(PL),__VA_ARGS__)
        
#define Plasticity_FindNbOfData(PL,...) \
        GenericData_FindNbOfData(Plasticity_GetGenericData(PL),__VA_ARGS__)

#define Plasticity_AppendData(PL,...) \
        Plasticity_AppendGenericData(PL,GenericData_Create(__VA_ARGS__))





#include "Elasticity.h"
#include "GenericData.h"
#include "Curves.h"
#include "Buffers.h"

struct Plasticity_s {
  char*   codenameofmodel ;
  double* dfsds ;  /** Yield function gradient */
  double* dgsds ;  /** Potential function gradient */
  double* hardv ;  /** Hardening variable */
  double* dhardv ; /** Typical small increment of hardening variable */
  double  dstress ; /** Typical small increment of stress */
  int     ncriteria ;  /** Nb of criteria */
  int     nhardv ;  /** Nb of hardening variables */
  int     nnhardv ;  /** Nb of non hardening variables */
  double* hardm ;  /** Hardening modulus */
  double* criterion ; /** Value of the yield function criterion */
  double* fc ;     /** fc(k,l) = dfsds(j,i) * C(i,j,k,l) */
  double* cg ;     /** cg(i,j) = C(i,j,k,l) * dgsds(l,k) */
  double* fcg ;    /** fcg = hardm + dfsds(j,i) * C(i,j,k,l) * dgsds(l,k)) */
  double* cijkl ;  /** Tangent stiffness tensor */
  double* lambda ; /** Plastic multiplier */
  double* parameter ;
  GenericData_t* genericdata ;  /** Plastic properties */
  Curves_t* curves ;          /**< Curves */
  Elasticity_t* elasty ;
  Plasticity_ComputeTangentStiffnessTensor_t* computetangentstiffnesstensor ;
  Plasticity_ReturnMapping_t*                 returnmapping ;
  Plasticity_YieldFunction_t*                 yieldfunction ;
  Plasticity_FlowRules_t*                     flowrules ;
  Plasticity_SetParameters_t*                 setparameters ;
  Plasticity_SetModelProp_t*                  setmodelprop ;
  Buffers_t*  buffers ;
} ;

#endif
