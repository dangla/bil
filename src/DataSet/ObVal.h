#ifndef OBVAL_H
#define OBVAL_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct ObVal_s        ; typedef struct ObVal_s        ObVal_t ;



#include "DataFile.h"

extern ObVal_t*  (ObVal_New)    (void) ;
extern void      (ObVal_Delete) (void*) ;
extern void      (ObVal_Scan)   (ObVal_t*,DataFile_t*) ;



#define ObVal_MaxLengthOfKeyWord        (30)


#define ObVal_GetType(OV)             ((OV)->type)
#define ObVal_GetNameOfUnknown(OV)    ((OV)->inc)
#define ObVal_GetValue(OV)            ((OV)->val)
#define ObVal_GetRelaxationFactor(OV) ((OV)->relaxfactor)



#define ObVal_IsRelativeValue(OV) \
        (ObVal_GetType(OV) == 'r')
        
#define ObVal_SetTypeToRelative(OV) \
        (ObVal_GetType(OV) = 'r')


#define ObVal_IsAbsoluteValue(OV) \
        (ObVal_GetType(OV) == 'a')
        
#define ObVal_SetTypeToAbsolute(OV) \
        (ObVal_GetType(OV) = 'a')


#include <math.h>

#define ObVal_GetAbsoluteValue(OV,U) \
        (ObVal_GetValue(OV) * ((ObVal_IsAbsoluteValue(OV)) ? 1 : fabs(U)))

#define ObVal_GetRelativeValue(OV,U) \
        (ObVal_GetValue(OV) / ((ObVal_IsRelativeValue(OV)) ? 1 : fabs(U)))




struct ObVal_s {              /* Objective variation */
  char    type ;              /* Type = a(bsolute) or r(elative) */
  char*   inc ;               /* Name of the unknown */
  double  val ;               /* Objective variation */
  double  relaxfactor ;       /* Relaxation factor */
} ;

#endif
