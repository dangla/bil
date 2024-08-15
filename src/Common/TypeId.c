#include <stdio.h>
#include <stdlib.h>

#include "TypeId.h"
#include "Mry.h"




TypeId_t* (TypeId_Create_)(const TypeId_e id,const size_t sz)
{
  TypeId_t* tid = (TypeId_t*) Mry_New(TypeId_t) ;
  
  TypeId_GetIdNumber(tid) = id ;
  TypeId_GetSize(tid) = sz ;
  
  return(tid) ;
}




void  (TypeId_Delete)(void* self)
{
  TypeId_t* tid = (TypeId_t*) self ;
  
  if(tid) {
  }

}

