#ifndef TUPLE_H
#define TUPLE_H

/*
 * Preprocessor tuple: adapted from 
 * http://cern.ch/laurent.deniau/cos.html
 */


/* NOTE-USER: definitions
   - a sequence is a list of comma separated tokens, e.g. a,b,c
   - a tuple is a sequence within parenthesis, e.g. (a,b,c)
     tuples are limited to 63 elements (C99 limits)
*/
 
#include "Arg.h"
#include "Utils.h"

/**
 * Primitives (split and rev could be implemented using eval)
 */

// length of tuple, i.e Tuple_LEN(()) is 1
#define Tuple_LEN(T) \
        Arg_NARG T

// expand tuple to sequence (e.g curry)
#define Tuple_SEQ(T) \
        Utils_ID T

// build tuple from sequence (e.g uncurry)
#define Tuple_TUPLE(...) \
        (__VA_ARGS__)

// eval function F() up to n times on tuple (run in O(n))
#define Tuple_EVAL(n,T,F) \
        Utils_CAT_(Tuple_EVAL_,n)(T,F)

// split a tuple at position n into a tuple of 2 tuples (T1,T2)
#define Tuple_SPLIT(n,T) \
        Utils_CAT_(Tuple_SPLIT_,n) T

// reverse tuple T
#define Tuple_REV(T) \
        Utils_CAT(Tuple_REV_,Arg_NARG T) T

/***********************************************************
 * Getters (assume that element exists)
 */

// elements accessors
#define Tuple_1ST(T)  Utils_PAIR(Utils_ARG1,(Utils_ID T,))
#define Tuple_2ND(T)  Utils_PAIR(Utils_ARG2,(Utils_ID T,))
#define Tuple_3RD(T)  Utils_PAIR(Utils_ARG3,(Utils_ID T,))
#define Tuple_4TH(T)  Utils_PAIR(Utils_ARG4,(Utils_ID T,))
#define Tuple_5TH(T)  Utils_PAIR(Utils_ARG5,(Utils_ID T,))
#define Tuple_6TH(T)  Utils_PAIR(Utils_ARG6,(Utils_ID T,))
#define Tuple_7TH(T)  Utils_PAIR(Utils_ARG7,(Utils_ID T,))
#define Tuple_8TH(T)  Utils_PAIR(Utils_ARG8,(Utils_ID T,))
#define Tuple_9TH(T)  Utils_PAIR(Utils_ARG9,(Utils_ID T,))
#define Tuple_LAST(T) Tuple_1ST(Tuple_REV(T))

// elements indexed accessor in range 1..length(T)
#define Tuple_ELEM(n,T) \
        Tuple_1ST(Tuple_2ND(Tuple_SPLIT(n,(,Utils_ID T))))

// take first n <= length(T) elements of tuple
#define Tuple_TAKE(n,T) \
        Tuple_1ST(Tuple_SPLIT(n,(Utils_ID T,)))

// take last n <= length(T) elements of tuple
#define Tuple_RTAKE(n,T) \
        Tuple_REV(Tuple_TAKE(n,Tuple_REV(T)))

// drop first n < length(T) elements of tuple
#define Tuple_DROP(n,T) \
        Tuple_2ND(Tuple_SPLIT(n,T))

// drop last n < length(T) elements of tuple
#define Tuple_RDROP(n,T) \
        Tuple_REV(Tuple_DROP(n,Tuple_REV(T)))

/***********************************************************
 * Setters
 */

// concat element a in front of tuple
#define Tuple_CONS(a,T) \
        (a,Utils_ID T)

// concat element a at the end of tuple (reverse CONS)
#define Tuple_RCONS(T,a) \
        (Utils_ID T,a)

// concat tuples T1 and T2
#define Tuple_CONCAT(T1,T2) \
        (Utils_ID T1,Utils_ID T2)

/***********************************************************
 * Sequences
 */

// sequence of n elements from duplication of va_args elements
#define Tuple_DUPSEQ(n,...) \
        Tuple_SEQ(Tuple_TAKE(n,(Arg_DUPSEQ_N(__VA_ARGS__))))

// numbered sequence of n elements from duplication of va_args elements
#define Tuple_NUMSEQ(n,...) \
        Tuple_SEQ(Tuple_TAKE(n,(Arg_NUMSEQ_N(__VA_ARGS__))))

// reverse numbered sequence of n elements from duplication of va_args elements
#define Tuple_REVSEQ(n,...) \
        Tuple_SEQ(Tuple_TAKE(n,(Arg_REVSEQ_N(__VA_ARGS__))))

/***********************************************************
 * Implementation
 */

// ----- EVAL -----

#define Tuple_EVAL_0( T,F) T
#define Tuple_EVAL_1( T,F) F(T)
#define Tuple_EVAL_2( T,F) F(F(T))
#define Tuple_EVAL_3( T,F) F(F(F(T)))
#define Tuple_EVAL_4( T,F) F(F(F(F(T))))
#define Tuple_EVAL_5( T,F) F(F(F(F(F(T)))))
#define Tuple_EVAL_6( T,F) F(F(F(F(F(F(T))))))
#define Tuple_EVAL_7( T,F) F(F(F(F(F(F(F(T)))))))
#define Tuple_EVAL_8( T,F) F(F(F(F(F(F(F(F(T))))))))
#define Tuple_EVAL_9( T,F) F(F(F(F(F(F(F(F(F(T)))))))))
#define Tuple_EVAL_10(T,F) F(F(F(F(F(F(F(F(F(F(T))))))))))
#define Tuple_EVAL_11(T,F) F(F(F(F(F(F(F(F(F(F(F(T)))))))))))
#define Tuple_EVAL_12(T,F) F(F(F(F(F(F(F(F(F(F(F(F(T))))))))))))
#define Tuple_EVAL_13(T,F) F(F(F(F(F(F(F(F(F(F(F(F(F(T)))))))))))))
#define Tuple_EVAL_14(T,F) F(F(F(F(F(F(F(F(F(F(F(F(F(F(T))))))))))))))
#define Tuple_EVAL_15(T,F) F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T)))))))))))))))
#define Tuple_EVAL_16(T,F) F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T))))))))))))))))
#define Tuple_EVAL_17(T,F) F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T)))))))))))))))))
#define Tuple_EVAL_18(T,F) \
                         F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T))))))))))))))))))
#define Tuple_EVAL_19(T,F) \
                      F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T)))))))))))))))))))
#define Tuple_EVAL_20(T,F) \
                   F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T))))))))))))))))))))
#define Tuple_EVAL_21(T,F) \
                F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(F(T)))))))))))))))))))))
#define Tuple_EVAL_22(T,F) Tuple_EVAL_1( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_23(T,F) Tuple_EVAL_2( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_24(T,F) Tuple_EVAL_3( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_25(T,F) Tuple_EVAL_4( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_26(T,F) Tuple_EVAL_5( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_27(T,F) Tuple_EVAL_6( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_28(T,F) Tuple_EVAL_7( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_29(T,F) Tuple_EVAL_8( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_30(T,F) Tuple_EVAL_9( Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_31(T,F) Tuple_EVAL_10(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_32(T,F) Tuple_EVAL_11(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_33(T,F) Tuple_EVAL_12(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_34(T,F) Tuple_EVAL_13(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_35(T,F) Tuple_EVAL_14(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_36(T,F) Tuple_EVAL_15(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_37(T,F) Tuple_EVAL_16(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_38(T,F) Tuple_EVAL_17(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_39(T,F) Tuple_EVAL_18(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_40(T,F) Tuple_EVAL_19(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_41(T,F) Tuple_EVAL_20(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_42(T,F) Tuple_EVAL_21(Tuple_EVAL_21(T,F),F)
#define Tuple_EVAL_43(T,F) Tuple_EVAL_1( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_44(T,F) Tuple_EVAL_2( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_45(T,F) Tuple_EVAL_3( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_46(T,F) Tuple_EVAL_4( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_47(T,F) Tuple_EVAL_5( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_48(T,F) Tuple_EVAL_6( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_49(T,F) Tuple_EVAL_7( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_50(T,F) Tuple_EVAL_8( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_51(T,F) Tuple_EVAL_9( Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_52(T,F) Tuple_EVAL_10(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_53(T,F) Tuple_EVAL_11(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_54(T,F) Tuple_EVAL_12(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_55(T,F) Tuple_EVAL_13(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_56(T,F) Tuple_EVAL_14(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_57(T,F) Tuple_EVAL_15(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_58(T,F) Tuple_EVAL_16(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_59(T,F) Tuple_EVAL_17(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_60(T,F) Tuple_EVAL_18(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_61(T,F) Tuple_EVAL_19(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_62(T,F) Tuple_EVAL_20(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)
#define Tuple_EVAL_63(T,F) Tuple_EVAL_21(Tuple_EVAL_21(Tuple_EVAL_21(T,F),F),F)

// ----- SPLIT -----

#define Tuple_SPLIT_0( ...)\
                      ((),(__VA_ARGS__))
#define Tuple_SPLIT_1( a,...)\
                      ((a),(__VA_ARGS__))
#define Tuple_SPLIT_2( a,b,...)\
                      ((a,b),(__VA_ARGS__))
#define Tuple_SPLIT_3( a,b,c,...)\
                      ((a,b,c),(__VA_ARGS__))
#define Tuple_SPLIT_4( a,b,c,d,...)\
                      ((a,b,c,d),(__VA_ARGS__))
#define Tuple_SPLIT_5( a,b,c,d,e,...)\
                      ((a,b,c,d,e),(__VA_ARGS__))
#define Tuple_SPLIT_6( a,b,c,d,e,f,...)\
                      ((a,b,c,d,e,f),(__VA_ARGS__))
#define Tuple_SPLIT_7( a,b,c,d,e,f,g,...)\
                      ((a,b,c,d,e,f,g),(__VA_ARGS__))
#define Tuple_SPLIT_8( a,b,c,d,e,f,g,h,...)\
                      ((a,b,c,d,e,f,g,h),(__VA_ARGS__))
#define Tuple_SPLIT_9( a,b,c,d,e,f,g,h,i,...)\
                      ((a,b,c,d,e,f,g,h,i),(__VA_ARGS__))
#define Tuple_SPLIT_10(a,b,c,d,e,f,g,h,i,j,...)\
                      ((a,b,c,d,e,f,g,h,i,j),(__VA_ARGS__))
#define Tuple_SPLIT_11(a,b,c,d,e,f,g,h,i,j,k,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k),(__VA_ARGS__))
#define Tuple_SPLIT_12(a,b,c,d,e,f,g,h,i,j,k,l,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l),(__VA_ARGS__))
#define Tuple_SPLIT_13(a,b,c,d,e,f,g,h,i,j,k,l,m,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m),(__VA_ARGS__))
#define Tuple_SPLIT_14(a,b,c,d,e,f,g,h,i,j,k,l,m,n,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n),(__VA_ARGS__))
#define Tuple_SPLIT_15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o),(__VA_ARGS__))
#define Tuple_SPLIT_16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p),(__VA_ARGS__))
#define Tuple_SPLIT_17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q),(__VA_ARGS__))
#define Tuple_SPLIT_18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r),(__VA_ARGS__))
#define Tuple_SPLIT_19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s),(__VA_ARGS__))
#define Tuple_SPLIT_20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t),(__VA_ARGS__))
#define Tuple_SPLIT_21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u),(__VA_ARGS__))
#define Tuple_SPLIT_22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v),(__VA_ARGS__))
#define Tuple_SPLIT_23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w),\
(__VA_ARGS__))
#define Tuple_SPLIT_24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x),\
(__VA_ARGS__))
#define Tuple_SPLIT_25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y),\
(__VA_ARGS__))
#define Tuple_SPLIT_26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,...)\
                      ((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z),\
(__VA_ARGS__))
#define Tuple_SPLIT_27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A),(__VA_ARGS__))
#define Tuple_SPLIT_28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B),(__VA_ARGS__))
#define Tuple_SPLIT_29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C),(__VA_ARGS__))
#define Tuple_SPLIT_30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D),(__VA_ARGS__))
#define Tuple_SPLIT_31(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E),(__VA_ARGS__))
#define Tuple_SPLIT_32(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F),(__VA_ARGS__))
#define Tuple_SPLIT_33(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G),(__VA_ARGS__))
#define Tuple_SPLIT_34(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H),(__VA_ARGS__))
#define Tuple_SPLIT_35(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I),(__VA_ARGS__))
#define Tuple_SPLIT_36(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J),(__VA_ARGS__))
#define Tuple_SPLIT_37(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K),(__VA_ARGS__))
#define Tuple_SPLIT_38(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L),(__VA_ARGS__))
#define Tuple_SPLIT_39(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M),(__VA_ARGS__))
#define Tuple_SPLIT_40(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N),(__VA_ARGS__))
#define Tuple_SPLIT_41(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O),(__VA_ARGS__))
#define Tuple_SPLIT_42(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P),(__VA_ARGS__))
#define Tuple_SPLIT_43(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q),(__VA_ARGS__))
#define Tuple_SPLIT_44(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R),(__VA_ARGS__))
#define Tuple_SPLIT_45(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S),(__VA_ARGS__))
#define Tuple_SPLIT_46(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T),(__VA_ARGS__))
#define Tuple_SPLIT_47(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U),(__VA_ARGS__))
#define Tuple_SPLIT_48(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V),(__VA_ARGS__))
#define Tuple_SPLIT_49(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W),(__VA_ARGS__))
#define Tuple_SPLIT_50(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X),(__VA_ARGS__))
#define Tuple_SPLIT_51(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y),(__VA_ARGS__))
#define Tuple_SPLIT_52(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z),(__VA_ARGS__))
#define Tuple_SPLIT_53(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1),(__VA_ARGS__))
#define Tuple_SPLIT_54(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1),(__VA_ARGS__))
#define Tuple_SPLIT_55(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1),(__VA_ARGS__))
#define Tuple_SPLIT_56(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1),(__VA_ARGS__))
#define Tuple_SPLIT_57(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1),(__VA_ARGS__))
#define Tuple_SPLIT_58(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,f1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1,f1),(__VA_ARGS__))
#define Tuple_SPLIT_59(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,f1,g1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1,f1,g1),(__VA_ARGS__))
#define Tuple_SPLIT_60(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,f1,g1,h1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1,f1,g1,h1),(__VA_ARGS__))
#define Tuple_SPLIT_61(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,f1,g1,h1,i1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1,f1,g1,h1,i1),(__VA_ARGS__))
#define Tuple_SPLIT_62(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1,f1,g1,h1,i1,j1),(__VA_ARGS__))
#define Tuple_SPLIT_63(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                        a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,...)\
((a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1),(__VA_ARGS__))

// ----- REV -----

#define Tuple_REV_1( a)\
                     (a)
#define Tuple_REV_2( a,b)\
                     (b,a)
#define Tuple_REV_3( a,b,c)\
                     (c,b,a)
#define Tuple_REV_4( a,b,c,d)\
                     (d,c,b,a)
#define Tuple_REV_5( a,b,c,d,e)\
                     (e,d,c,b,a)
#define Tuple_REV_6( a,b,c,d,e,f)\
                     (f,e,d,c,b,a)
#define Tuple_REV_7( a,b,c,d,e,f,g)\
                     (g,f,e,d,c,b,a)
#define Tuple_REV_8( a,b,c,d,e,f,g,h)\
                     (h,g,f,e,d,c,b,a)
#define Tuple_REV_9( a,b,c,d,e,f,g,h,i)\
                     (i,h,g,f,e,d,c,b,a)
#define Tuple_REV_10(a,b,c,d,e,f,g,h,i,j)\
                     (j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_11(a,b,c,d,e,f,g,h,i,j,k)\
                     (k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_12(a,b,c,d,e,f,g,h,i,j,k,l)\
                     (l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_13(a,b,c,d,e,f,g,h,i,j,k,l,m)\
                     (m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_14(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
                     (n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)\
                     (o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\
                     (p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)\
                     (q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r)\
                     (r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s)\
                     (s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t)\
                     (t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)\
                     (u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v)\
                     (v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w)\
                     (w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x)\
                     (x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y)\
                     (y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)\
                     (z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A)\
(A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B)\
(B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C)\
(C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D)\
(D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_31(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E)\
(E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_32(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F)\
(F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_33(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G)\
(G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_34(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H)\
(H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_35(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I)\
(I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_36(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J)\
(J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_37(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K)\
(K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_38(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L)\
(L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_39(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M)\
(M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_40(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N)\
(N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_41(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O)\
(O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_42(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P)\
(P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_43(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q)\
(Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a
#define Tuple_REV_44(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R)\
(R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_45(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S)\
(S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_46(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T)\
(T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_47(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U)\
(U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_48(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V)\
(V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_49(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W)\
(W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_50(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X)\
(X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_51(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y)\
(Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_52(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z)\
(Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_53(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1)\
(a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_54(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                   A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                   a1,b1)\
(b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_55(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1)\
(c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_56(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1)\
(d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_57(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1)\
(e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_58(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1,f1)\
(f1,e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_59(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1,f1,g1)\
(g1,f1,e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_60(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1,f1,g1,h1)\
(h1,g1,f1,e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_61(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1,f1,g1,h1,i1)\
(i1,h1,g1,f1,e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_62(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1,f1,g1,h1,i1,j1)\
(j1,i1,h1,g1,f1,e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define Tuple_REV_63(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,\
                      a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1)\
(k1,j1,i1,h1,g1,f1,e1,d1,c1,b1,a1,\
Z,Y,X,W,V,U,T,S,R,Q,P,O,N,M,L,K,J,I,H,G,F,E,D,C,B,A,\
z,y,x,w,v,u,t,s,r,q,p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)

#endif
