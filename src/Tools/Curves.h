#ifndef CURVES_H
#define CURVES_H



/* vacuous declarations and typedef names */

/* class-like structures */
struct Curves_s       ; typedef struct Curves_s       Curves_t ;
/*   1. Curves_t attributes */
struct Curve_s        ; typedef struct Curve_s        Curve_t ;


extern Curves_t* (Curves_Create)(unsigned int) ;
extern void      (Curves_Delete)(Curves_t**) ;
extern int       (Curves_ReadCurves)(Curves_t*,const char*) ;
//extern int       (Curves_WriteCurves1)(char*) ;
//extern int       (Curves_WriteCurves2)(char*) ;
extern int       (Curves_Append)(Curves_t*,Curve_t*) ;
extern int       (Curves_FindCurveIndex)(Curves_t*,const char*) ;
//#define Curves_WriteCurves Curves_WriteCurves2


#define Curves_MaxLengthOfTextLine      (500)
#define Curves_SizeOfBuffer             (Curves_MaxLengthOfTextLine*sizeof(char))


#define Curves_GetNbOfAllocatedCurves(curves)    ((curves)->n_allocatedcurves)
#define Curves_GetNbOfCurves(curves)             ((curves)->n_cb)
#define Curves_GetCurve(curves)                  ((curves)->cb)
#define Curves_GetBuffer(curves)                 ((curves)->buffer)


#define Curves_AllocateInBuffer(curves,sz) \
        (Buffer_Allocate(Curves_GetBuffer(curves),(sz)))
        
#define Curves_FreeBuffer(curves) \
        (Buffer_Free(Curves_GetBuffer(curves)))
        
#define Curves_CreateDerivative(curves,cv) \
        (Curves_Append(curves,Curve_CreateDerivative(cv)))
        
#define Curves_CreateIntegral(curves,cv) \
        (Curves_Append(curves,Curve_CreateIntegral(cv)))
        
#define Curves_CreateInverse(curves,cv,sc) \
        (Curves_Append(curves,Curve_CreateInverse(cv,sc)))

#define Curves_CannotAppendCurves(curves,i) \
        (Curves_GetNbOfCurves(curves) + i > Curves_GetNbOfAllocatedCurves(curves))


extern Curve_t* (Curve_Create)(unsigned int) ;
extern void     (Curve_Delete)(Curve_t**) ;
extern Curve_t* (Curve_CreateDerivative)(Curve_t*) ;
extern Curve_t* (Curve_CreateIntegral)(Curve_t*) ;
extern Curve_t* (Curve_CreateInverse)(Curve_t*,const char) ;
extern double*  (Curve_CreateSamplingOfX)(Curve_t*) ;
extern double   (Curve_ComputeValue)(Curve_t*,double) ;
extern double   (Curve_ComputeDerivative)(Curve_t*,double) ;
extern double   (Curve_ComputeIntegral)(Curve_t*,double) ;
extern char*    (Curve_PrintInFile)(Curve_t*) ;


#define Curve_MaxLengthOfKeyWord        (30)
#define Curve_MaxLengthOfFileName       (60)
#define Curve_MaxLengthOfTextLine       (500)
#define Curve_MaxLengthOfCurveName      (30)


#define Curve_GetNbOfPoints(curve)      ((curve)->n)
#define Curve_GetXRange(curve)          ((curve)->a)
#define Curve_GetYValue(curve)          ((curve)->f)
#define Curve_GetScaleType(curve)       ((curve)->echelle)
#define Curve_GetNameOfXAxis(curve)     ((curve)->xname)
#define Curve_GetNameOfYAxis(curve)     ((curve)->yname)



#include "Buffer.h"

struct Curves_s {             /* Curves */
  unsigned int n_allocatedcurves ;         /* Nb of allocated curves */
  unsigned int n_cb ;         /* nb of curves */
  Curve_t *cb ;               /* curves */
  Buffer_t   *buffer ;        /* Buffer */
} ;

struct Curve_s {              /* courbe */
  char*  xname ;              /* Name of the x-axis */
  char*  yname ;              /* Name of the y-axis */
  char   echelle ;            /* echelle = n(ormale) ou l(ogarithmique) */
  int    n ;                  /* nombre de points */
  double* a ;                 /* abscisses */
  double* f ;                 /* valeurs f(a) */
} ;



/* Old notations which should be eliminated */
#define crbe_t                 Curve_t
#define courbe(a,b)            Curve_ComputeValue(&(b),(a))
#define dcourbe(a,b)           Curve_ComputeDerivative(&(b),(a))
#define lit_courbe(mat,b)      Curves_ReadCurves(Material_GetCurves(mat),(b))
#define ecrit_courbe           Curves_WriteCurves

#endif
