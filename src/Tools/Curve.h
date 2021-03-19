#ifndef CURVE_H
#define CURVE_H



/* vacuous declarations and typedef names */

/* class-like structures */
struct Curve_s        ; typedef struct Curve_s        Curve_t ;



extern Curve_t* (Curve_Create)(unsigned int) ;
extern void     (Curve_Delete)(void*) ;
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



/* Set the names of x-axis and y-axis */
#define Curve_SetNameOfXAxis(CV,NAME) \
        strcpy(Curve_GetNameOfXAxis(CV),NAME)

#define Curve_SetNameOfYAxis(CV,NAME) \
        strcpy(Curve_GetNameOfYAxis(CV),NAME)



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

#endif
