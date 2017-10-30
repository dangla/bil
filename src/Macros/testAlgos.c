#include "Algos.h"
#include "Codes.h"


#define T1 (A,b,C,d,E,f,G,h,I,j)

#define T2 (1,2,3,4,5,6,7,8,9,10,11,12)

//#define T2 (K,L,M,N,O,P,Q,R,S,T,U,V)

#define T3 (aA,bB,cC,dD)

#define Fa(a)        Utils_CAT(t_,Codes_ISLOWER(a))

#define Ga(a)        Utils_CAT(_,a)

#define Fab(a,b)     Utils_CAT3(a,_,b)

#define Fabc(a,b,c)  Utils_CAT3(a,b,c)

#define PF(a)        Codes_ISUPPER(a)

#define T4 double

#define T5 long,double

#define T6 long double


/* test Algos_MAP */
#define TEST1 \
        Algos_FOLDL(T1,,Fab)

#define TEST2 \
        Algos_FOLDR(T1,_P,Fab)

#define TEST3 \
        Algos_SCANL(T1,U,Fab)

#define TEST4 \
        Algos_SCANR(T2,U,Fab)

#define TEST5 \
        Algos_FILTER(T1,PF)

#define TEST6 \
        Algos_MAP(T1,Ga)

#define TEST7 \
        Algos_MAP2(T1,T2,Fab)

#define TEST8 \
        Algos_MAP3(T1,T2,T3,Fabc)

#define TEST9 \
        Algos_SEP(T1) = \
        Algos_SEPWITH(T1, )

#define TEST10 \
        Algos_SEPWITH(T3,@)
        
#define TEST11(...) \
        Utils_CAT(TypeId,Algos_FOLDL(Tuple_TUPLE(__VA_ARGS__),,Fab))
        
#define TEST12(...) \
        TypeId_Create(__VA_ARGS__)
        
#define TypeId_Create(...) \
        Tuple_SEQ(TypeId_Create_(__VA_ARGS__))
        
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



#define DEBUG 1


#if 0
#include <stdio.h>
int main(void)
{
  printf("TEST1  = %s\n",Utils_STR(TEST1)) ;
  printf("TEST2  = %s\n",Utils_STR(TEST2)) ;
  printf("TEST3  = %s\n",Utils_STR(TEST3)) ;
  printf("TEST4  = %s\n",Utils_STR(TEST4)) ;
  printf("TEST5  = %s\n",Utils_STR(TEST5)) ;
  printf("TEST6  = %s\n",Utils_STR(TEST6)) ;
  printf("TEST7  = %s\n",Utils_STR(TEST7)) ;
  printf("TEST8  = %s\n",Utils_STR(TEST8)) ;
  printf("TEST9  = %s\n",Utils_STR(TEST9)) ;
  printf("TEST10 = %s\n",Utils_STR(TEST10)) ;
}
#endif

#if 0
static char test1[] = Utils_STR(TEST1) ;
static char test2[] = Utils_STR(TEST2) ;
static char test3[] = Utils_STR(TEST3) ;
static char test4[] = Utils_STR(TEST4) ;
static char test5[] = Utils_STR(TEST5) ;
static char test6[] = Utils_STR(TEST6) ;
static char test7[] = Utils_STR(TEST7) ;
static char test8[] = Utils_STR(TEST8) ;
static char test9[] = Utils_STR(TEST9) ;
static char test10[] = Utils_STR(TEST10) ;
#endif

#if 0
"Algos_FOLDL(T1,P_,Fab)" = {
  TEST1
} ;
"TEST2" = {
  TEST2
} ;

"TEST3" = {
  TEST3
} ;

"TEST4" = {
  TEST4
} ;

"TEST5" = {
  TEST5
} ;

"TEST6" = {
  TEST6
} ;

"TEST7" = {
  TEST7
} ;

"TEST8" = {
  TEST8
} ;

"TEST9" = {
  TEST9
} ;

"TEST10" = {
  TEST10
} ;

"TEST11" = {
  TEST11(T5)
} ;

"TEST12" = {
  TEST12(unsigned int)
} ;
#endif
