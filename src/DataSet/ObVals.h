#ifndef OBVALS_H
#define OBVALS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct ObVals_s       ; typedef struct ObVals_s       ObVals_t ;
/*   1. ObVals_t  attributes */
struct ObVal_s        ; typedef struct ObVal_s        ObVal_t ;



/* 1. ObVals_t 
 * -----------*/
#include "DataFile.h"
#include "Mesh.h"

extern ObVals_t*  ObVals_Create(DataFile_t*,Mesh_t*,Materials_t*) ;

#define ObVals_GetNbOfObVals(OVS)    ((OVS)->n_obj)
#define ObVals_GetObVal(OVS)         ((OVS)->obj)




/* 2. ObVal_t 
 * ----------*/
#define ObVal_MaxLengthOfKeyWord        (30)

#define ObVal_GetType(OV)             ((OV)->type)
#define ObVal_GetNameOfUnknown(OV)    ((OV)->inc)
#define ObVal_GetValue(OV)            ((OV)->val)
#define ObVal_GetRelaxationFactor(OV) ((OV)->relaxfactor)



#define ObVal_IsRelativeValue(OV) \
        (ObVal_GetType(OV) == 'r')

#define ObVal_IsAbsoluteValue(OV) \
        (ObVal_GetType(OV) == 'a')



struct ObVals_s {             /* objective variations */
  unsigned int n_obj ;        /* nb */
  ObVal_t* obj ;              /* objective variation */
} ;

struct ObVal_s {              /* Objective variation */
  char    type ;              /* Type = a(bsolute) or r(elative) */
  char*   inc ;               /* Name of the unknown */
  double  val ;               /* Objective variation */
  double  relaxfactor ;       /* Relaxation factor */
} ;

#endif
