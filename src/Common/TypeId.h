#ifndef TYPEID_H
#define TYPEID_H


/* Type identifiers */
enum TypeId_e {
  /* From C */
  TypeId_undefined,
  TypeId_unsigned_char,
  TypeId_char,
  TypeId_unsigned_int,
  TypeId_short_int,
  TypeId_int,
  TypeId_unsigned_long,
  TypeId_long_int,
  TypeId_float,
  TypeId_double,
  TypeId_long_double,
  /* From Bil */
  TypeId_Solution_t,
  TypeId_Solutions_t,
  TypeId_last
} ;

typedef enum TypeId_e     TypeId_t ;
        

#include "Utils.h"
#include "Tuple.h"
#include "Algos.h"


/* T  stands for a real type */
/* ID stands for a TypeId */


/* Create a TypeId */
#define TypeId_Create(T) \
        Tuple_SEQ(TypeId_Create_(T))


/* Test the type */
#define TypeId_Is(ID,T) \
        (ID == TypeId_Create(T))
        
#define TypeId_SetTo(ID,T) \
        do {ID = TypeId_Create(T) ;} while(0)


/* Implementation */

#define TypeId_Create_(...) \
        (Utils_CAT(TypeId_,__VA_ARGS__))
        
#define TypeId_long \
        TypeId_long_(
        
#define TypeId_long_(...) \
        Utils_CAT(TypeId_long_,__VA_ARGS__))
        
#define TypeId_short \
        TypeId_short_(
        
#define TypeId_short_(...) \
        Utils_CAT(TypeId_short_,__VA_ARGS__))
        
#define TypeId_unsigned \
        TypeId_unsigned_(
        
#define TypeId_unsigned_(...) \
        Utils_CAT(TypeId_unsigned_,__VA_ARGS__))
        
        


#if 0
/* Qualifiers must be used with a comma: TypeId_GetSizeOf(long,double) */
#define TypeId_GetSizeOf(...) \
        (sizeof(Algos_SEP(Tuple_TUPLE(__VA_ARGS__))))

#define TypeId_Create(...) \
        Utils_CAT(TypeId,Algos_FOLDL(Tuple_TUPLE(__VA_ARGS__),,TypeId_CAT))
        
#define TypeId_CAT(a,b) \
        Utils_CAT3(a,_,b)
#endif


#endif
