#ifndef GENERICDATA_H
#define GENERICDATA_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct GenericData_s     ; typedef struct GenericData_s     GenericData_t ;



/* GenericData_t
 * -------------*/


#include <stdio.h>
#include "TypeId.h"

extern GenericData_t* (GenericData_New)          (void) ;
extern void           (GenericData_Delete)       (GenericData_t**) ;
extern GenericData_t* (GenericData_Create_)      (int,size_t,TypeId_t) ;
extern void           (GenericData_Initialize_)  (GenericData_t*,int,void*,TypeId_t) ;
extern void           (GenericData_InsertBefore) (GenericData_t*,GenericData_t*) ;
extern void           (GenericData_InsertAfter)  (GenericData_t*,GenericData_t*) ;
extern GenericData_t* (GenericData_First)        (GenericData_t*) ;
extern GenericData_t* (GenericData_Last)         (GenericData_t*) ;



#define GenericData_GetTypeId(GD)                ((GD)->typ)
#define GenericData_GetNbOfData(GD)              ((GD)->n)
#define GenericData_GetData(GD)                  ((GD)->data)
#define GenericData_GetNextGenericData(GD)       ((GD)->next)
#define GenericData_GetPreviousGenericData(GD)   ((GD)->prev)



        
#define GenericData_Create(N,T) \
        GenericData_Create_(N,sizeof(T),TypeId_Create(T))
        
#define GenericData_Initialize(A,B,C,T) \
        GenericData_Initialize_(A,B,C,TypeId_Create(T))
        
        
#define GenericData_Append(A,B) \
        GenericData_InsertAfter(GenericData_Last(A),B)



/* Typedef names of methods */
//typedef void* GenericData_Allocate(int,TypedId_t) ;




/* Generic data */
struct GenericData_s {
  TypeId_t typ ;                /* The type id of data */
  int n ;                       /* Nb of data */
  void* data ;                  /* The data */
  GenericData_t* prev ;         /* Previous generic data */
  GenericData_t* next ;         /* Next generic data */
} ;


#endif
