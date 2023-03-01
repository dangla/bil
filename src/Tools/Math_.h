#ifndef MATH_H
#define MATH_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Math_s         ; typedef struct Math_s         Math_t ;


/* Constants */
#define Math_Pi           (3.14159265358979323846)
#define Math_Ln10         (2.3025850929940456840179914546843642076011)

extern void     (Math_Delete)(void*) ;
extern double   (Math_ComputeFirstStressInvariant)(const double*) ;
extern double   (Math_ComputeSecondStressInvariant)(const double*) ;
extern double   (Math_ComputeSecondDeviatoricStressInvariant)(const double*) ;
extern double*  (Math_SolveByGaussEliminationWithPartialPivoting)(double*,double*,int) ;
extern double*  (Math_SolveByGaussJordanElimination)(double*,double*,int,int) ;
extern int      (Math_ComputePolynomialEquationRoots)(double*,int) ;
extern int      (Math_PolishPolynomialEquationRoot)(double*,int,double*,double,int) ;
extern double   (Math_EvaluateExpression)(char*) ;
extern double   (Math_EvaluateExpressions)(char*,char*) ;
extern double   (Math_Compute3x3MatrixDeterminant)(const double*) ;
extern double*  (Math_Inverse3x3Matrix)(const double*) ;
extern double*  (Math_ComputePrincipalStresses)(const double*) ;
extern double*  (Math_ComputeRealEigenvaluesAndEigenvectorsOf3x3Matrix)(double*,const char) ;
extern void     (Math_PrintStiffnessTensor)(const double*) ;
extern void     (Math_PrintStressTensor)(const double*) ;
extern void     (Math_PrintMatrix)(const double*,const int) ;
extern double*  (Math_ComputeDeviatoricStress)(const double*) ;
extern double*  (Math_SolveByGaussEliminationKIJ)(double*,double*,int,int*) ;
extern double*  (Math_SolveByGaussEliminationJIK)(double*,double*,int,int*) ;





#define Math_SolveByGaussEliminationWithoutPartialPivoting(...) \
        Math_SolveByGaussEliminationKIJ(__VA_ARGS__,NULL)

#define Math_SolveByGaussElimination(...) \
        Utils_CAT_NARG(Math_SolveByGaussElimination,__VA_ARGS__)(__VA_ARGS__)

/* Implementation */
#define Math_SolveByGaussElimination3 \
        Math_SolveByGaussEliminationWithPartialPivoting
        
#define Math_SolveByGaussElimination4 \
        Math_SolveByGaussEliminationKIJ


#include "BilExtraLibs.h"

#ifdef LAPACKLIB
#if defined(__cplusplus)
  extern "C" {
#endif

extern int dgeev_(const char*,const char*,int*,double*,int*,double*, double*,double*,int*,double*,int*,double*,int*,int*);

#if defined(__cplusplus)
  }
#endif
#endif


#define Math_ComputeSecondDeviatoricStrainInvariant \
        Math_ComputeSecondDeviatoricStressInvariant
        
#define Math_ComputeDeviatoricStrain \
        Math_ComputeDeviatoricStress
        
#define Math_ComputeMatrixInverse(a,n) \
        Math_SolveByGaussJordanElimination(a,NULL,n,0)
        
#define Math_InvertMatrix(A,N)  \
        Math_SolveByGaussJordanElimination(A,NULL,N,0)

#define Math_Max(a,b)          (((a) > (b)) ? (a) : (b))
#define Math_Min(a,b)          (((a) < (b)) ? (a) : (b))

#define Math_Swap(a,b,Type_t) \
        do { \
          Type_t tmp = (a) ; \
          (a) = (b) ; \
          (b) = tmp ; \
        } while(0)
        
#define Math_SwapDouble(a,b)   Math_Swap(a,b,double)
#define Math_SwapInt(a,b)      Math_Swap(a,b,int)

#define Math_Sign(a)           ((a < 0) ? -1 : (a > 0))


#ifndef M_PI
#define M_PI       Math_Pi
#endif

#ifndef MAX
#define MAX        Math_Max
#endif

#ifndef MIN
#define MIN        Math_Min
#endif



#define Math_MaxNbOfMatrices                      (4)
#define Math_MaxSizeOfMatrix                      (9*sizeof(double))
#define Math_SizeOfBuffer                         (Math_MaxNbOfMatrices*Math_MaxSizeOfMatrix)


#define Math_GetOutputs(M)             ((M)->outputs)
#define Math_GetBuffers(M)             ((M)->buffers)
#define Math_GetDelete(M)              ((M)->Delete)



#define Math_GetBuffer(M) \
        Buffers_GetBufferOfCurrentThread(Math_GetBuffers(M))

/* Operations on buffer */
#define Math_AllocateInBuffer(M,sz) \
        (Buffer_Allocate(Math_GetBuffer(M),(sz)))
        
#define Math_FreeBuffer(M) \
        (Buffer_Free(Math_GetBuffer(M)))
        
#define Math_FreeBufferFrom(M,p) \
        (Buffer_FreeFrom(Math_GetBuffer(M),(char*) (p)))



#include "Buffers.h"
#include "GenericObject.h"

struct Math_s {          /* some math methods */
  void* outputs ;
  Buffers_t*  buffers ;         /* Buffer */
  GenericObject_Delete_t* Delete ;
} ;



/* Old notations that I try to eliminate little by little */
#define j2       Math_ComputeSecondDeviatoricStressInvariant
//#define gausse   Math_SolveByGaussElimination

#endif
