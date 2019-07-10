#ifndef ELASTICITY_H
#define ELASTICITY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Elasticity_s     ; typedef struct Elasticity_s     Elasticity_t ;


/* 1. Elasticity_t */

extern Elasticity_t*  (Elasticity_Create)(void) ;
extern void           (Elasticity_Delete)(void*) ;
extern void           (Elasticity_SetParameters)           (Elasticity_t*,...) ;
extern void           (Elasticity_SetParameter)            (Elasticity_t*,const char*,double) ;
extern double*        (Elasticity_ComputeStiffnessTensor)  (Elasticity_t*,double*) ;
extern void           (Elasticity_PrintStiffnessTensor)    (Elasticity_t*) ;


/* Accessors */
#define Elasticity_GetStiffnessTensor(EL)      ((EL)->cijkl)
#define Elasticity_GetType(EL)                 ((EL)->type)
#define Elasticity_GetParameter(EL)            ((EL)->parameter)


#define Elasticity_MaxLengthOfKeyWord     (100)
#define Elasticity_MaxNbOfParameters      (8)


#include "Tools/Math.h"


/* Isotropy */
#define Elasticity_IsIsotropic(EL) \
        (!strcmp(Elasticity_GetType(EL),"isotropy"))

#define Elasticity_SetToIsotropy(EL) \
        Elasticity_CopyType(EL,"isotropy")
        
#define Elasticity_GetYoungModulus(EL) \
        Elasticity_GetParameter(EL)[0]

#define Elasticity_GetPoissonRatio(EL) \
        Elasticity_GetParameter(EL)[1]
        
#define Elasticity_GetBulkModulus(EL) \
        Elasticity_GetParameter(EL)[2]

#define Elasticity_GetShearModulus(EL) \
        Elasticity_GetParameter(EL)[3]


/* Transversely isotropy */
#define Elasticity_IsTransverselyIsotropic(EL) \
        (!strcmp(Elasticity_GetType(EL),"transiso"))
        
#define Elasticity_SetToTransverselyIsotropy(EL) \
        Elasticity_CopyType(EL,"transiso")
        
#define Elasticity_GetYoungModulus3(EL) \
        Elasticity_GetParameter(EL)[4]

#define Elasticity_GetPoissonRatio3(EL) \
        Elasticity_GetParameter(EL)[5]
        
#define Elasticity_GetShearModulus3(EL) \
        Elasticity_GetParameter(EL)[6]
        
#define Elasticity_GetAxis3(EL) \
        Elasticity_GetParameter(EL)[7]
        
        

/* Implementation */
#define Elasticity_CopyType(EL,typ) \
        memcpy(Elasticity_GetType(EL),typ,MAX(strlen(typ),Elasticity_MaxLengthOfKeyWord))



struct Elasticity_s {
  double* cijkl ;
  char*   type ;
  double* parameter ;
} ;

#endif
