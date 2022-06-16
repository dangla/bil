#ifndef STACK_H
#define STACK_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Stack_s     ; typedef struct Stack_s     Stack_t ;



/* Stack_t
 * -------*/


#include <stdio.h>
#include "TypeId.h"

extern Stack_t*       (Stack_New)          (void) ;
extern void           (Stack_Delete)       (void*) ;
extern Stack_t*       (Stack_Push_)        (Stack_t*,int,void*,TypeId_t*) ;
extern Stack_t*       (Stack_Pop)          (Stack_t*) ;



#define Stack_GetTypeId(STK)                ((STK)->typ)
#define Stack_GetNbOfData(STK)              ((STK)->n)
#define Stack_GetData(STK)                  ((STK)->data)
#define Stack_GetBelow(STK)                 ((STK)->below)



        
#define Stack_Push(N,T) \
        Stack_Push_(N,sizeof(T),TypeId_Create(T))


#define Stack_IsEmpty(STK) \
        (!STK)




/* Stack */
struct Stack_s {
  TypeId_t* typ ; 
  int n ;
  void* data ;
  Stack_t* below ;
} ;


#endif
