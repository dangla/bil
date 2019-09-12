#ifndef MATH_H
#define MATH_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Math_s         ; typedef struct Math_s         Math_t ;


/* Constants */
#define Math_Pi           (3.14159265358979323846)
#define Math_Ln10         (2.3025850929940456840179914546843642076011)

extern double   (Math_ComputeSecondDeviatoricStressInvariant)(const double*) ;
extern double*  (Math_SolveByGaussElimination)(double*,double*,int) ;
extern int      (Math_ComputePolynomialEquationRoots)(double*,int) ;
extern int      (Math_PolishPolynomialEquationRoot)(double*,int,double*,double,int) ;
extern double   (Math_EvaluateExpression)(char*) ;
extern double   (Math_EvaluateExpressions)(char*,char*) ;


#define Math_ComputeSecondDeviatoricStrainInvariant \
        Math_ComputeSecondDeviatoricStressInvariant
        
#define Math_ComputeMatrixInverse(a,n) \
        Math_SolveByGaussJordanElimination(a,NULL,n,0)

#define Math_Max(a,b)          (((a) > (b)) ? (a) : (b))
#define Math_Min(a,b)          (((a) < (b)) ? (a) : (b))
#define Math_Swap(a,b,Type_t)  {Type_t tmp = (a) ; (a) = (b) ; (b) = tmp ;}
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


#define Math_GetOutputs(math)             ((math)->outputs)


struct Math_s {          /* some math methods */
  void* outputs ;
} ;



/* Old notations that I try to eliminate little by little */
#define j2       Math_ComputeSecondDeviatoricStressInvariant
//#define gausse   Math_SolveByGaussElimination

#endif
