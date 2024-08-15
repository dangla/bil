#ifndef INTFCT_H
#define INTFCT_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure "DataSet_t and attributes */
struct IntFct_s       ; typedef struct IntFct_s       IntFct_t ;


extern IntFct_t*  (IntFct_New)(void) ;
extern IntFct_t*  (IntFct_Create)(int,int,const char*) ;
extern void       (IntFct_Delete)(void*) ;
//extern void       (IntFct_ComputeIsoShapeFctInActualSpace)(int,int,double**,double*,int,double*,double*) ;
extern int        (IntFct_ComputeFunctionIndexAtPointOfReferenceFrame)(IntFct_t*,double*) ;
extern double     (IntFct_InterpolateAtPoint)(IntFct_t*,double*,int,int) ;


#define IntFct_MaxLengthOfKeyWord          (30)

#define IntFct_MaxNbOfIntPoints            (20)

#include "ShapeFct.h"

#define IntFct_MaxNbOfFunctions            (ShapeFct_MaxNbOfNodes)


#define IntFct_GetType(IFCT)             ((IFCT)->type)
#define IntFct_GetNbOfNodes(IFCT)        ((IFCT)->nn)
#define IntFct_GetNbOfFunctions(IFCT)    ((IFCT)->nn)
#define IntFct_GetNbOfPoints(IFCT)       ((IFCT)->np)
#define IntFct_GetWeight(IFCT)           ((IFCT)->w)
#define IntFct_GetFunction(IFCT)         ((IFCT)->h)
#define IntFct_GetFunctionGradient(IFCT) ((IFCT)->dh)
#define IntFct_GetDimension(IFCT)        ((IFCT)->dim)
#define IntFct_GetPointCoordinates(IFCT) ((IFCT)->a)
//#define IntFct_GetComputeValuesAtPoint(IFCT)        ((IFCT)->computevaluesatpoint)



#define IntFct_GetCoordinatesAtPoint(IFCT,p) \
        (IntFct_GetPointCoordinates(IFCT) + (p)*3)
//        (IntFct_GetPointCoordinates(IFCT) + (p)*IntFct_GetDimension(IFCT))
        
#define IntFct_GetFunctionAtPoint(IFCT,p) \
        (IntFct_GetFunction(IFCT) + (p)*IntFct_MaxNbOfFunctions)
//        (IntFct_GetFunction(IFCT) + (p)*IntFct_GetNbOfFunctions(IFCT))

#define IntFct_GetFunctionGradientAtPoint(IFCT,p) \
        (IntFct_GetFunctionGradient(IFCT) + \
        (p)*3*IntFct_MaxNbOfFunctions)
//        (p)*IntFct_GetDimension(IFCT)*IntFct_GetNbOfFunctions(IFCT))
        
#define IntFct_TypeIs(IFCT,type) \
        (!strcmp(IntFct_GetType(IFCT),type))
        
#define IntFct_SetType(IFCT,type) \
        (strcpy(IntFct_GetType(IFCT),type))




/*  Typedef names of Methods */
//typedef double* IntFct_ComputeValuesAtPoint_t(IntFct_t*,double*) ;


struct IntFct_s {             /* Interpolation function */
  char*   type ;              /* Type of the function */
  unsigned short int nn ;     /* Number of functions/nodes */
  unsigned short int np ;     /* Number of integration points */
  double* a ;                 /* Reference coordinates of integration points */
  double* w ;                 /* Weights */
  double* h ;                 /* Values of interpolation functions */
  double* dh ;                /* Values of function gradients */
  unsigned short int dim ;    /* sous-dimension (0,1,2,3) */
  //IntFct_ComputeValuesAtPoint_t* computevaluesatpoint ;
} ;


/* Old notations which should be eliminated */
#define inte_t    IntFct_t
#define MAX_PGAUSS                    IntFct_MaxNbOfIntPoints
//#define fint_abs                      IntFct_ComputeIsoShapeFctInActualSpace


#ifdef __CPLUSPLUS
}
#endif
#endif
