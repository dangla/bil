#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "TypeId.h"
#include "Stack.h"




/* Global functions */

Stack_t* (Stack_New)(void)
{
  Stack_t* stack = (Stack_t*) malloc(sizeof(Stack_t)) ;
  
  if(!stack) {
    arret("Stack_New(1)") ;
  }
  
  {
    Stack_GetTypeId(stack) = TypeId_Create(undefined) ;
    Stack_GetNbOfData(stack) = 0 ;
    Stack_GetData(stack)  = NULL ;
    Stack_GetBelow(stack) = NULL ;
  }
  
  return(stack) ;
}



void (Stack_Delete)(void* self)
{
  Stack_t* stack = (Stack_t*) self ;
  
  if(stack) {
    Stack_t* below = Stack_GetBelow(stack) ;
    void* data = Stack_GetData(stack) ;
      
    free(data) ;
    stack = below ;
  }
}



Stack_t* (Stack_Push_)(Stack_t* head,int n,void* data,TypeId_t typ)
{
  Stack_t* stack = Stack_New() ;
  
  {
    Stack_GetTypeId(stack)   = typ ;
    Stack_GetNbOfData(stack) = n ;
    Stack_GetData(stack)     = data ;
    Stack_GetBelow(stack)    = head ;
  }
  
  return(stack) ;
}



Stack_t* (Stack_Pop)(Stack_t* head)
{
  if(head) {
    Stack_t* below = Stack_GetBelow(head) ;
    void* data = Stack_GetData(head) ;
      
    free(data) ;
    return(below) ;
  }
  
  return(head) ;
}
