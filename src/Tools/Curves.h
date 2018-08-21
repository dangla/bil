#ifndef CURVES_H
#define CURVES_H



/* vacuous declarations and typedef names */

/* class-like structures */
struct Curves_s       ; typedef struct Curves_s       Curves_t ;


#include "Curve.h"

extern Curves_t* (Curves_Create)(unsigned int) ;
extern void      (Curves_Delete)(void*) ;
extern int       (Curves_ReadCurves)(Curves_t*,const char*) ;
//extern int       (Curves_WriteCurves1)(char*) ;
//extern int       (Curves_WriteCurves2)(char*) ;
extern int       (Curves_Append)(Curves_t*,Curve_t*) ;
extern int       (Curves_FindCurveIndex)(Curves_t*,const char*) ;
extern Curve_t*  (Curves_FindCurve)(Curves_t*,const char*) ;
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



#include "Buffer.h"


struct Curves_s {             /* Curves */
  unsigned int n_allocatedcurves ;         /* Nb of allocated curves */
  unsigned int n_cb ;         /* nb of curves */
  Curve_t *cb ;               /* curves */
  Buffer_t   *buffer ;        /* Buffer */
} ;



/* Old notations which should be eliminated */
#define lit_courbe(mat,b)      Curves_ReadCurves(Material_GetCurves(mat),(b))
#define ecrit_courbe           Curves_WriteCurves

#endif
