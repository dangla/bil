#include <stdio.h>
#include <stdlib.h>

#include "TypeId.h"
#include "Mry.h"




TypeId_t* (TypeId_Create_)(const int id,const size_t sz)
{
  TypeId_t* typeid = (TypeId_t*) Mry_New(TypeId_t) ;
  
  TypeId_GetIdNumber(typeid) = id ;
  TypeId_GetSize(typeid) = sz ;
  
  return(typeid) ;
}




void  (TypeId_Delete)(void* self)
{
  TypeId_t* typeid = (TypeId_t*) self ;
  
  if(typeid) {
  }

}

