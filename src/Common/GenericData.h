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
extern void           (GenericData_Delete)       (void*) ;
extern GenericData_t* (GenericData_Create_)      (int,void*,TypeId_t,char const*) ;
extern void           (GenericData_Initialize_)  (GenericData_t*,int,void*,TypeId_t,char const*) ;
//extern void           (GenericData_InsertBefore) (GenericData_t*,GenericData_t*) ;
//extern void           (GenericData_InsertAfter)  (GenericData_t*,GenericData_t*) ;
extern GenericData_t* (GenericData_Append)       (GenericData_t*,GenericData_t*) ;
extern GenericData_t* (GenericData_First)        (GenericData_t*) ;
extern GenericData_t* (GenericData_Last)         (GenericData_t*) ;
extern GenericData_t* (GenericData_Find_)        (GenericData_t*,TypeId_t,char const*) ;




#define GenericData_MaxLengthOfKeyWord           (30)



#define GenericData_GetTypeId(GD)                ((GD)->typ)
#define GenericData_GetName(GD)                  ((GD)->name)
#define GenericData_GetNbOfData(GD)              ((GD)->n)
#define GenericData_GetData(GD)                  ((GD)->data)
#define GenericData_GetNextGenericData(GD)       ((GD)->next)
#define GenericData_GetPreviousGenericData(GD)   ((GD)->prev)
#define GenericData_GetDelete(GD)                ((GD)->Delete)



        
#define GenericData_Create(A,B,T,S) \
        GenericData_Create_(A,B,TypeId_Create(T),S)
        
#define GenericData_Initialize(GD,A,B,T,S) \
        GenericData_Initialize_(GD,A,B,TypeId_Create(T),S)
        
#define GenericData_Find(GD,T,N) \
        GenericData_Find_(GD,TypeId_Create(T),N)
        
        
#define GenericData_Is(GD,I,S) \
        ((GenericData_GetTypeId(GD) == I) && !strcmp(GenericData_GetName(GD),S))
        
        
#define GenericData_FindData(GD,T,N) \
        GenericData_GetData(GenericData_Find(GD,T,N))
        
        
#define GenericData_Merge(A,B) \
        GenericData_Append(A,GenericData_First(B))



/* Typedef names of methods */
//typedef void* GenericData_Allocate(int,TypedId_t) ;



#include <GenericObject.h>

/* Generic data */
struct GenericData_s {
  TypeId_t typ ;                /* The type id of data */
  char*    name ;               /* Name of the data */
  int n ;                       /* Nb of data */
  void* data ;                  /* The data */
  GenericData_t* prev ;         /* Previous generic data */
  GenericData_t* next ;         /* Next generic data */
  GenericObject_Delete_t* Delete ;
} ;


#endif
