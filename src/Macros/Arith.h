#ifndef ARITH_H
#define ARITH_H

/*
 * Preprocessor arithmetic macros: adapted from 
 * http://cern.ch/laurent.deniau/cos.html
 */

#include "Arg.h"
#include "Utils.h"
#include "Tuple.h"
#include "Logic.h"

// increment integer 0 <= n <= Arg_MAX_N, saturates at Arg_MAX_N
#define Arith_INCR(n) \
        Tuple_1ST(Tuple_2ND(Tuple_SPLIT(n,(Arg_NUMSEQ_N(),Arg_MAX_N))))

// decrement integer 0 <= n <= Arg_MAX_N, saturates at zero
#define Arith_DECR(n) \
        Tuple_1ST(Tuple_2ND(Tuple_SPLIT(n,(0,0,Arg_NUMSEQ_N()))))

// add two integers m >= 0, n >= 0, m+n <= Arg_MAX_N
#define Arith_ADD(m,n) \
        Tuple_1ST(Tuple_2ND(Tuple_SPLIT(n, \
          Tuple_2ND(Tuple_SPLIT(m,(0,Arg_NUMSEQ_N()))))))

// substract two integers m >= 0, n >= 0, , m <= n <= Arg_MAX_N
#define Arith_SUB(m,n) \
        Logic_IF(Logic_ISZERO(m))(0, \
          Tuple_1ST(Tuple_2ND(Tuple_SPLIT(n,Tuple_RCONS(Tuple_REV( \
            Tuple_1ST(Tuple_SPLIT(m,(Arg_NUMSEQ_N(),)))),0)))))

// multiply two integers m >= 0, n >= 0, (m,n,m*n)<=Arg_MAX_N
#define Arith_MUL(m,n) \
        Arg_IF(Logic_OR(Logic_ISZERO(m),Logic_ISZERO(n)))(0, \
          Arg_NARG(Utils_PAIR(Utils_REST,(Utils_DUP(m,Utils_DUP(n,,))))))

#endif
