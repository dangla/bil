#ifndef OBVALS_H
#define OBVALS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct ObVals_s       ; typedef struct ObVals_s       ObVals_t ;



#include "DataFile.h"
#include "Mesh.h"
#include "Materials.h"


extern ObVals_t*  (ObVals_New)(const int) ;
extern ObVals_t*  (ObVals_Create)(DataFile_t*,Mesh_t*,Materials_t*) ;
extern int        (ObVals_FindObValIndex)(ObVals_t*,char*) ;



#define ObVals_GetNbOfObVals(OVS)    ((OVS)->n_obj)
#define ObVals_GetObVal(OVS)         ((OVS)->obj)



#include "ObVal.h"


struct ObVals_s {             /* objective variations */
  unsigned int n_obj ;        /* nb */
  ObVal_t* obj ;              /* objective variation */
} ;


#endif
