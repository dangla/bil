#ifndef MRY_H
#define MRY_H


#include <stdlib.h>

extern void*     (Mry_Allocate)(size_t) ;
extern void*     (Mry_AllocateZeroed)(size_t,size_t) ;
extern void*     (Mry_Realloc)(void*,size_t) ;
extern void      (Mry_Free)(void*) ;


#include <stdarg.h>
#include "Utils.h"

#define Mry_New(...) \
        Utils_CAT_NARG(Mry_New_,__VA_ARGS__)(__VA_ARGS__)


/* Mry_New(T[N]) is also valid */
#define Mry_New_1(T) \
        Mry_AllocateZeroed(1,sizeof(T))

#define Mry_New_2(T,N) \
        Mry_AllocateZeroed((size_t) (N),sizeof(T))
        

#endif
