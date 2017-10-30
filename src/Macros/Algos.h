#ifndef ALGOS_H
#define ALGOS_H

/*
 * Preprocessor algos: adapted from 
 * http://cern.ch/laurent.deniau/cos.html
 */


#include "Arg.h"
#include "Utils.h"
#include "Tuple.h"
#include "Logic.h"


// fold left elements of tuple using function F(a,b)
#define Algos_FOLDL(T,a0,F) \
        Tuple_1ST(Tuple_EVAL(Tuple_LEN(T), \
          (a0,(Utils_ID T,),F),Algos_FOLDL_0))

// fold right elements of tuple using function F(a,b)
#define Algos_FOLDR(T,a0,F) \
        Tuple_1ST(Tuple_EVAL(Tuple_LEN(T), \
          (a0,(Tuple_SEQ(Tuple_REV(T)),),F),Algos_FOLDR_0))

// scan left elements of tuple using function F(a,b)
#define Algos_SCANL(T,a0,F) \
        Algos_RRES_(Tuple_EVAL(Tuple_LEN(T), \
          ((a0,),(Utils_ID T,),F),Algos_SCANL_0))

// scan right elements of tuple using function F(a,b)
#define Algos_SCANR(T,a0,F) \
        Tuple_REV(Algos_RRES_(Tuple_EVAL(Tuple_LEN(T), \
          ((a0,),(Tuple_SEQ(Tuple_REV(T)),),F),Algos_SCANR_0)))

// filter elements of tuple using the predicate function PF(A)
#define Algos_FILTER(T,PF) \
        Algos_RES1_(Tuple_EVAL(Tuple_LEN(T), \
          ((),(Utils_ID T,),PF),Algos_FILTER_0))

// map function F(a) on elements of tuple T
#define Algos_MAP(T,F) \
        Algos_RES_(Tuple_EVAL(Tuple_LEN(T), \
          ((),(Utils_ID T,),F),Algos_MAP_0))

// map function F(a1,a2) on elements of tuples T1,T2 up to length(T1)
#define Algos_MAP2(T1,T2,F) \
        Algos_RES_(Tuple_EVAL(Tuple_LEN(T1), \
          ((),(Utils_ID T1,),(Utils_ID T2,),F),Algos_MAP2_0))

// map function F(a1,a2,a3) on elements of tuples T1,T2,T3 up to length(T1)
#define Algos_MAP3(T1,T2,T3,F) \
        Algos_RES_(Tuple_EVAL(Tuple_LEN(T1), \
          ((),(Utils_ID T1,),(Utils_ID T2,),(Utils_ID T3,),F),Algos_MAP3_0))

// flatten tuple T to sequence using s as separator, i.e. (a,b,c) -> a b c
#define Algos_SEP(T) \
        Algos_FOLDL(T,,Utils_PAIR)

// flatten tuple T to sequence using s as separator, i.e. (a,b,c),; -> a ; b ; c
#define Algos_SEPWITH(T,s) \
        Logic_IF(Logic_2ARGS(Utils_ID T))( \
          Utils_ARG1 T Algos_SEP( \
          Algos_MAP2((Utils_REST T),(Arg_DUPSEQ_N(s)),Utils_SWAP)), \
          Utils_ID T)

/***********************************************************
 * Implementation
 */

#define Algos_RES1_(T) \
        Algos_RES1__ T
#define Algos_RES1__(R,...) \
        Logic_IF(Logic_2ARGS(Utils_ID R))((Utils_REST R),R)

#define Algos_RES_(T) \
        Algos_RES__ T
#define Algos_RES__(R,...) \
        (Utils_PAIR(Utils_REST,R))

#define Algos_RRES_(T) \
        Algos_RRES__ T
#define Algos_RRES__(R,...) \
        (Utils_PAIR(Utils_REST,Tuple_REV(R)))

#define Algos_FOLDL_0(T) \
        Algos_FOLDL_1 T
#define Algos_FOLDL_1(R,T,F) \
        (F(R,Utils_ARG1 T),(Utils_REST T),F)

#define Algos_FOLDR_0(T) \
        Algos_FOLDR_1 T
#define Algos_FOLDR_1(R,T,F) \
        (F(Utils_ARG1 T,R),(Utils_REST T),F)

#define Algos_SCANL_0(T) \
        Algos_SCANL_1 T
#define Algos_SCANL_1(R,T,F) \
        ((F(Utils_ARG1 R,Utils_ARG1 T),Utils_ID R),(Utils_REST T),F)

#define Algos_SCANR_0(T) \
        Algos_SCANR_1 T
#define Algos_SCANR_1(R,T,F) \
        ((F(Utils_ARG1 T,Utils_ARG1 R),Utils_ID R),(Utils_REST T),F)

#define Algos_FILTER_0(T) \
        Algos_FILTER_1 T
#define Algos_FILTER_1(R,T,PF) \
        (Logic_IF(PF(Utils_ARG1 T))((Utils_ID R,Utils_ARG1 T),R),(Utils_REST T),PF)

#define Algos_MAP_0(T) \
        Algos_MAP_1 T
#define Algos_MAP_1(R,T,F) \
        ((Utils_ID R,F(Utils_ARG1 T)),(Utils_REST T),F)

#define Algos_MAP2_0(T) \
        Algos_MAP2_1 T
#define Algos_MAP2_1(R,T1,T2,F) \
        ((Utils_ID R,F(Utils_ARG1 T1,Utils_ARG1 T2)),(Utils_REST T1),(Utils_REST T2),F)

#define Algos_MAP3_0(T) \
        Algos_MAP3_1 T
#define Algos_MAP3_1(R,T1,T2,T3,F) \
        ((Utils_ID R,F(Utils_ARG1 T1,Utils_ARG1 T2,Utils_ARG1 T3)), \
          (Utils_REST T1),(Utils_REST T2),(Utils_REST T3),F)

#endif
