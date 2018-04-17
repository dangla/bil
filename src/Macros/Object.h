#ifndef OBJECT_H
#define OBJECT_H

#include <stdlib.h>
#include <stdarg.h>


/* T stands or a type */
extern void* Object_New_(const int = 1) ;
template <typename T> inline void* Object_New_(const int n)
{
  T* A = new T[n] ;
  //T* A = (T*) calloc(n,sizeof(T)) ;
  if(A) assert(A) ;
  return(A) ;
}
#define Object_New(T, ...) \
        Object_New_<T>(__VA_ARGS__)


/* OBJ stands for a pointer to an object. 
 * (OBJ)->delete should point to a function which is expected 
 *  to free the space allocated for the creation of OBJ. */
#define Object_Delete(OBJ) \
        (*(OBJ))->Delete(OBJ)


#if 0
#include "Utils.h"

#define Object_Type(OBJ) \
        Utils_CAT(OBJ,_t)
#endif
        
#endif
