#ifndef SESSION_H
#define SESSION_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structures */
struct Session_s       ; typedef struct Session_s       Session_t ;



#include "GenericData.h"
#include "TypeId.h"

extern Session_t*     (Session_GetCurrentInstance)(void) ;
extern Session_t*     (Session_New)(void) ;
extern void           (Session_Delete)(void) ;
extern Session_t*     (Session_Open)(void) ;
extern Session_t*     (Session_Close)(void) ;
//extern GenericData_t* (Session_FindGenericData)(Session_t*,TypeId_t*,char const*) ;
//extern void           (Session_AddGenericData)(Session_t*,TypeId_t*,char const*) ;


#define Session_GetIndex(SS)                ((SS)->index)
#define Session_GetGenericData(SS)          ((SS)->gdat)
#define Session_GetPreviousSession(SS)      ((SS)->prev)


  
/* The current instance */
#define Session_CurrentInstance \
        Session_GetCurrentInstance()


/* The generic data of the current instance */
#define Session_GenericData \
        Session_GetGenericData(Session_CurrentInstance)


/* Find and add generic data in the current instance */
#define Session_FindGenericData(T,N) \
        GenericData_Find(Session_GenericData,T,N)
        
        
#define Session_AddGenericData(GD) \
        do { \
          Session_GenericData = GenericData_Append(Session_GenericData,GD); \
        } while(0)
        
/*
#define Session_MergeGenericData(SS,GD) \
        (Session_GetGenericData(SS) = GenericData_Merge(Session_GetGenericData(SS),GD))
*/




struct Session_s {
  int            index ;
  GenericData_t* gdat ;
  Session_t*     prev ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
