#ifndef MRY_H
#define MRY_H


#include <stdlib.h>

extern void*     (Mry_Malloc)(size_t) ;
extern void*     (Mry_Calloc)(size_t,size_t) ;
extern void*     (Mry_Realloc)(void*,size_t) ;
extern void      (Mry_Free)(void*) ;



#include "Utils.h"

#define Mry_New(...) \
        Utils_CAT_NARG(Mry_New_,__VA_ARGS__)(__VA_ARGS__)
        
#define Mry_New_1(T) \
        Mry_New_2(T,1)
        
#define Mry_New_2(T,N) \
        (T*) Mry_Calloc(N,sizeof(T))
        

#endif
