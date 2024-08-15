#ifndef GENERICDATA_H
#define GENERICDATA_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct GenericData_s     ; typedef struct GenericData_s     GenericData_t ;



/* GenericData_t
 * -------------*/


#include <stdio.h>
#include "TypeId.h"

extern GenericData_t* (GenericData_New)          (void) ;
extern void           (GenericData_Delete)       (void*) ;
extern GenericData_t* (GenericData_Create_)      (int,void*,TypeId_t*,const char*) ;
extern void           (GenericData_Initialize_)  (GenericData_t*,int,void*,TypeId_t*,const char*) ;
//extern void           (GenericData_InsertBefore) (GenericData_t*,GenericData_t*) ;
//extern void           (GenericData_InsertAfter)  (GenericData_t*,GenericData_t*) ;
extern GenericData_t* (GenericData_Append)       (GenericData_t*,GenericData_t*) ;
extern GenericData_t* (GenericData_First)        (GenericData_t*) ;
extern GenericData_t* (GenericData_Last)         (GenericData_t*) ;
extern GenericData_t* (GenericData_Find_)        (GenericData_t*,const int,const char*) ;




#define GenericData_MaxLengthOfKeyWord           (30)



#define GenericData_GetTypeId(GD)                ((GD)->typ)
#define GenericData_GetName(GD)                  ((GD)->name)
#define GenericData_GetNbOfData(GD)              ((GD)->n)
#define GenericData_GetData(GD)                  ((GD)->data)
#define GenericData_GetNextGenericData(GD)       ((GD)->next)
#define GenericData_GetPreviousGenericData(GD)   ((GD)->prev)
#define GenericData_GetDelete(GD)                ((GD)->Delete)



#define GenericData_GetSize(GD) \
        TypeId_GetSize(GenericData_GetTypeId(GD))

#define GenericData_GetIdNumber(GD) \
        TypeId_GetIdNumber(GenericData_GetTypeId(GD))



        
#define GenericData_Create(N,DATA,T,NAME) \
        GenericData_Create_(N,DATA,TypeId_Create(T),NAME)
        
#define GenericData_Initialize(GD,N,DATA,T,NAME) \
        GenericData_Initialize_(GD,N,DATA,TypeId_Create(T),NAME)
        
#define GenericData_Find(GD,T,NAME) \
        GenericData_Find_(GD,TypeId_IdNumber(T),NAME)
        
        
#define GenericData_Is(GD,ID,NAME) \
        ((GenericData_GetIdNumber(GD) == ID) && !strcmp(GenericData_GetName(GD),NAME))
        
        
#define GenericData_FindData(GD,T,NAME) \
        ((T*) GenericData_GetData_(GenericData_Find(GD,T,NAME)))
        
        
#define GenericData_FindNbOfData(GD,T,NAME) \
        GenericData_GetNbOfData_(GenericData_Find(GD,T,NAME))
        
        
#define GenericData_Merge(A,B) \
        GenericData_Append(A,GenericData_First(B))





/* Implementation */
#define GenericData_GetNbOfData_(GD) \
        ((GD) ? GenericData_GetNbOfData(GD) : 0)

#define GenericData_GetData_(GD) \
        ((GD) ? GenericData_GetData(GD) : NULL)



/* Typedef names of methods */
//typedef void* GenericData_Allocate(int,TypedId_t) ;



#include "GenericObject.h"

/* Generic data */
struct GenericData_s {
  TypeId_t* typ ;               /* The type id of data */
  char*    name ;               /* Name of the data */
  int n ;                       /* Nb of data */
  void* data ;                  /* The data */
  GenericData_t* prev ;         /* Previous generic data */
  GenericData_t* next ;         /* Next generic data */
  GenericObject_Delete_t* Delete ;
} ;



#ifdef __CPLUSPLUS
}
#endif
#endif
