#ifndef GENERICOBJECT_H
#define GENERICOBJECT_H

/* 
 * After:
 * A. Schreiner, Object oriented programming with ANSI-C, 2011.
 */

/* vacuous declarations and typedef names */

/* class-like structure */
//struct GenericObject_s     ; typedef const struct GenericObject_s     GenericObject_t ;



/*  Typedef names of Methods */
typedef void (GenericObject_Delete_t)(void*) ;



/* GenericObject_t
 * ---------------*/
 

/* OBJ stands for a pointer to an object. 
 * (OBJ)->delete should point to a function which is expected 
 *  to free the space allocated for the creation of OBJ. */
#define GenericObject_Delete(OBJ) \
        (*(OBJ))->Delete(OBJ)



#endif
#if 0

template <typename T> extern void* GenericObject_New_(const int = 1) ;
template <typename T> inline void* GenericObject_New_(const int n)
{
  T* A = new T[n] ;
  //T* A = (T*) calloc(n,sizeof(T)) ;
  if(A) assert(A) ;
  return(A) ;
}

#define GenericObject_New(T, ...) \
        GenericObject_New_<T>(__VA_ARGS__)


#endif
#if 0

#include "Utils.h"
#include "TypeId.h"


#define GenericObject_Type(GO) \
        Utils_CAT(GO,_t)


//#define GenericObject_Delete(GO,...) \
//        TypeId_Delete(TypeId_Create(GenericObject_Type(GO)),__VA_ARGS__)

 
 
 
extern void*            (GenericObject_Ctor)   (const void*, ...) ;
extern void             (GenericObject_Dtor)   (void*) ;
extern size_t           (GenericObject_SizeOf) (const void*) ;
extern void*            (GenericObject_New_)   (const void*,const int) ;
extern const void*      (GenericObject_ClassOf)(const void*) ;







#define GenericObject_ClassInstance(GO) \
        Utils_CAT(GO,_ClassInstance)


#define GenericObject_Get(GO,A) \
        Utils_CAT3(GO,_Get,A)


#define GenericObject_SizeOf(GO,A) \
        ((size_t) sizeof(GenericObject_Get(GO,A)((GenericObject_Type(GO)*)0)))


#define GenericObject_OffSetOf(GO,A) \
        ((size_t) offsetof(GenericObject_Type(GO),GenericObject_Get(GO,A)((GenericObject_Type(GO)*)0)))


#define GenericObject_Create(GO,...) \
        GenericObject_Ctor(GenericObject_ClassInstance(GO),__VA_ARGS__)


#define GenericObject_New(GO,...) \
        GenericObject_New_(GenericObject_ClassInstance(GO),__VA_ARGS__)


#define GenericObject_Delete(...) \
        GenericObject_Dtor(__VA_ARGS__)


#define GenericObject_Size(GO) \
        ((size_t) sizeof(GenericObject_Type(GO)))
        
        
        

/* Declaration and definition of OB_GetInstance */
#define GenericObject_DeclareGetInstance(GO,...) \
        void*  (Utils_CAT(GO,_GetInstance))(__VA_ARGS__)


#define GenericObject_DefineGetInstance(GO,...) \
        GenericObject_DeclareGetInstance(GO,__VA_ARGS__) \
        { \
          GenericData_t* gdat = Session_FindGenericData(GenericObject_Type(GO),Utils_STR(GO)) ; \
          if(!gdat) { \
            GenericObject_Type(GO)* obj = Utils_CAT(GO,_Create)(__VA_ARGS__) ; \
            gdat = GenericData_Create(1,obj,GenericObject_Type(GO),Utils_STR(GO)) ; \
            Session_AddGenericData(gdat) ; \
            assert(gdat == Session_FindGenericData(GenericObject_Type(GO),Utils_STR(GO))) ; \
          } \
          return((GenericObject_Type(GO)*) GenericData_GetData(gdat)) ; \
        }


#endif
