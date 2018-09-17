#if 0

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "Class.h"
#include "GenericObject.h"





/* Global functions */

void* (GenericObject_Ctor)(const void* class0, ...)
{
  const Class_t* class1 = (const Class_t*) class0 ;
  void* p  = calloc(1,Class_GetSize(class1)) ;
  
  assert(p) ;
  
  /* Remind (Programming languages - C, ISO/IEC 9899:TC2 p. 103): 
   * A pointer to a structure object, suitably converted, points to its
   * initial member. So if TYP is the structure type pointed to by p,
   * i.e. obj = (TYP*) p, and typ that of the first member of obj,
   * then         (typ*) obj  = &((typ) obj->first_member)
   * and so    * ((typ*) obj) =   (typ) obj->first_member.
   * We force a conversion of p which treats the beginning of the object 
   * as a pointer to Class_t and set the argument class as the value 
   * of this pointer. */

  * (const Class_t**) p = class1 ;
  
  if(Class_GetCtor(class1)) {
    va_list ap ;
    
    va_start(ap,class0) ;
    p = Class_GetCtor(class1)(p,&ap) ;
    va_end(ap) ;
  }
  
  return(p) ;
}


/* Example of programming of "Obj_Ctor" in the body file Obj.c
void* Obj_Ctor(void* _self,va_list* app)
{
  const char* str = va_arg(*app,const char*) ;

  {
    Obj_t* self = _self ;

    Obj_Initialize(self,str) ;

    return(self) ;
  }
}
*/



void* (GenericObject_New_)(const void* class0,const int n)
{
  const Class_t* class1 = (const Class_t*) class0 ;
  void* p  = calloc(n,Class_GetSize(class1)) ;
  
  assert(p) ;
  
  * (const Class_t**) p = class1 ;
  
  return(p) ;
}



void (GenericObject_Dtor)(void* self)
{
  const Class_t** cp = (const Class_t**) self ;
  
  if(self && (*cp) && Class_GetDtor(*cp)) {
    self = Class_GetDtor(*cp)(self) ;
  }

  free(self);
}


/* Example of programming of "Obj_Dtor" in the body file Obj.c
void* Obj_Dtor(void* _self)
{
  Obj_t* self = _self ;

  Obj_Delete(&self) ;

  return(self) ;
}
*/



size_t (GenericObject_SizeOf)(const void* self)
{
  const Class_t* const* cp = (const Class_t* const*) self ;
  
  assert(self && (*cp)) ;
  
  return(Class_GetSize(*cp)) ;
}



const void* (GenericObject_ClassOf)(const void* self)
{
  const Class_t* const* cp = (const Class_t* const*) self ;
  
  assert(self && (*cp)) ;
  
  return(*cp) ;
}

#endif
