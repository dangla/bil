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
        Utils_CAT_NARG(Mry_New,__VA_ARGS__)(__VA_ARGS__)


/* Mry_New(T[N]) is also valid */
#define Mry_New1(T) \
        Mry_AllocateZeroed(1,sizeof(T))

#define Mry_New2(T,N) \
        Mry_AllocateZeroed((size_t) (N),sizeof(T))


#define Mry_Create(T,N,CREATE) \
        ({ \
          T* Mry_v = (T*) Mry_New(T,N) ; \
          do { \
            int Mry_i ; \
            for(Mry_i = 0 ; Mry_i < N ; Mry_i++) { \
              T* Mry_v0 = CREATE ; \
              Mry_v[Mry_i] = Mry_v0[0] ; \
              free(Mry_v0) ; \
            } \
          } while(0); \
          Mry_v ; \
        })


#define Mry_Delete(OBJ,N,DELETE) \
        do { \
          if(OBJ) { \
            int Mry_i ; \
            for(Mry_i = 0 ; Mry_i < N ; Mry_i++) { \
              DELETE(OBJ + Mry_i) ; \
            } \
          } \
        } while(0)
        

#endif
