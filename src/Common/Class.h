#ifndef CLASS_H
#define CLASS_H

#endif
#if 0

/* 
 * After:
 * A. Schreiner, Object oriented programming with ANSI-C, 2011.
 */

/* vacuous declarations and typedef names */

/* class-like structure */
struct Class_s     ; typedef struct Class_s     Class_t ;



/* Class_t
 * -------*/

#include <stdarg.h>
#include <stdio.h>

typedef void*    (Class_Ctor_t)  (void*,va_list*) ;
typedef void*    (Class_Dtor_t)  (void*) ;
//typedef void*    (Class_New_t)   (void*,const int) ;





#define Class_GetSize(CL)                ((CL)->size)
#define Class_GetCtor(CL)                ((CL)->ctor)
#define Class_GetDtor(CL)                ((CL)->dtor)



#define Class_Ctor(A,...) \
        Class_GetCtor((Class_t*) A)(A,__VA_ARGS__)
        
#define Class_Dtor(A) \
        Class_GetDtor((Class_t*) A)(A)

#define Class_New(A,...) \
        Class_GetNew((Class_t*) A)(A,__VA_ARGS__)
        

#include "Utils.h"

/* This is how to initialize any object of type OB_t */
#define Class_Instance(OB) \
        {sizeof(Utils_CAT(OB,_t)), \
         Utils_CAT(OB,_Ctor), \
         Utils_CAT(OB,_Dtor)}


/* To be included in the body file "OB.c" */
#define Class_Define(OB) \
        static const Class_t Utils_CAT(OB,_ClassInstance_) = Class_Instance(OB) ; \
        const void* Utils_CAT(OB,_ClassInstance) = &Utils_CAT(OB,_ClassInstance_) ;
/*
        void* (Utils_CAT(OB,_Dtor))(void* self) \
        { \
          Utils_CAT(OB,_t)* obj = (Utils_CAT(OB,_t)*) self ; \
          Utils_CAT(OB,_Delete)(&obj) ; \
          return(obj) ; \
        } \
        void* (Utils_CAT(OB,_New))(const int n) \
        { \
          Utils_CAT(OB,_t)* obj = (Utils_CAT(OB,_t)*) malloc(n*sizeof(Utils_CAT(OB,_t))) ; \
          assert(obj) ; \
          return(obj) ; \
        }
*/





/* To be included in the header file "OB.h" */
#define Class_Declare(OB) \
        extern const void* Utils_CAT(OB,_ClassInstance) ;
/*
        extern Class_New_t  Utils_CAT(OB,_New) ; \
*/





/* Class */
struct Class_s {
  size_t size ;
  Class_Ctor_t*    ctor ;
  Class_Dtor_t*    dtor ;
  //Class_Clone_t*   clone ;
  //Class_Differ_t*  differ ;
} ;

#endif
