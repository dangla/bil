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
extern int       (Curves_CreateDerivative)(Curves_t*,Curve_t*) ;
extern int       (Curves_CreateIntegral)(Curves_t*,Curve_t*) ;
extern int       (Curves_CreateInverse)(Curves_t*,Curve_t*,const char) ;


#define Curves_MaxLengthOfTextLine      (500)
#define Curves_SizeOfBuffer             (Curves_MaxLengthOfTextLine*sizeof(char))


#define Curves_GetNbOfAllocatedCurves(CVS)    ((CVS)->n_allocatedcurves)
#define Curves_GetNbOfCurves(CVS)             ((CVS)->n_cb)
#define Curves_GetCurve(CVS)                  ((CVS)->cb)
#define Curves_GetBuffer(CVS)                 ((CVS)->buffer)


#define Curves_AllocateInBuffer(CVS,SZ) \
        (Buffer_Allocate(Curves_GetBuffer(CVS),(SZ)))
        
#define Curves_FreeBuffer(CVS) \
        (Buffer_Free(Curves_GetBuffer(CVS)))



/* Below also works */
#if 0
#define Curves_CreateDerivative(CVS,CV) \
        ({ \
          int Curves_i ; \
          do { \
            Curve_t* dcv = Curve_CreateDerivative(CV) ; \
            Curves_i = Curves_Append(CVS,dcv) ; \
            free(dcv) ; \
          } while(0) ; \
          Curves_i ; \
        })
        
#define Curves_CreateIntegral(CVS,CV) \
        ({ \
          int Curves_i ; \
          do { \
            Curve_t* icv = Curve_CreateIntegral(CV) ; \
            Curves_i = Curves_Append(CVS,icv) ; \
            free(icv) ; \
          } while(0) ; \
          Curves_i ; \
        })
        
#define Curves_CreateInverse(CVS,CV,SCALE) \
        ({ \
          int Curves_i ; \
          do { \
            Curve_t* icv = Curve_CreateInverse(CV,SCALE) ; \
            Curves_i = Curves_Append(CVS,icv) ; \
            free(icv) ; \
          } while(0) ; \
          Curves_i ; \
        })
#endif



#define Curves_CannotAppendCurves(CVS,i) \
        (Curves_GetNbOfCurves(CVS) + i > Curves_GetNbOfAllocatedCurves(CVS))



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
