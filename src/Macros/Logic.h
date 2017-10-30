#ifndef LOGIC_H
#define LOGIC_H

/*
 * Preprocessor logic/relational/predicate macros: adapted from 
 * http://cern.ch/laurent.deniau/cos.html
 */
 
 
#include "Arg.h"
#include "Utils.h"
#include "Tuple.h"


// the following macros expect b{1,2} to be 0 or 1
// see Logic_BOOL for conversion

// evaluate first element of the following tuple, otherwise evaluate the rest
#define Logic_IF(b) \
        Utils_CAT_(Logic_IF_,b)

#define Logic_IFDEF(m) \
        Logic_IF(Logic_NOT(Logic_ISBLANK(m)))

#define Logic_IFNDEF(m) \
        Logic_IF(Logic_ISBLANK(m))

#define Logic_NOT(b) \
        Utils_CAT_(Logic_NOT_,b)()

#define Logic_AND(b1,b2) \
        Utils_CAT3_(Logic_AND_,b1,b2)()

#define Logic_OR(b1,b2) \
        Utils_CAT3_(Logic_OR_,b1,b2)()

#define Logic_XOR(b1,b2) \
        Utils_CAT3_(Logic_XOR_,b1,b2)()

// the following macros expect 0 <= m,n <= 63

// m >= n
#define Logic_GE(m,n) \
        Logic_ISTUPLE(Tuple_1ST(Tuple_2ND( \
          Tuple_SPLIT(n,((),Tuple_DUPSEQ(m,()),Arg_DUPSEQ_N())) )))

// m <= n
#define Logic_LE(m,n) \
        Logic_GE(n,m)

// m > n
#define Logic_GT(m,n) \
        Logic_NOT(Logic_GE(n,m))

// m < n
#define Logic_LT(m,n) \
        Logic_NOT(Logic_GE(m,n))

// m == n
#define Logic_EQ(m,n) \
        Logic_AND(Logic_GE(m,n),Logic_GE(n,m))

// min(m,n)
#define Logic_MIN(m,n) \
        Logic_IF(Logic_GE(m,n))(n,m)

// max(m,n)
#define Logic_MAX(m,n) \
        Logic_IF(Logic_GE(m,n))(m,n)

// the following macros expect 'a' to be concatenable (not a symbol)

// return 1 if 'a' is non-zero, including numbers and tokens, 0 otherwise
#define Logic_BOOL(a) \
        Logic_NOT(Logic_OR(Logic_ISZERO(a),Logic_ISBLANK(a)))

// return 1 if 'a' is 0 or 1, 0 otherwise
#define Logic_ISBOOL(a) \
        Logic_ISTUPLE(Utils_CAT_(Logic_ISBOOL_,a))

// return 1 if 'a' is blank (empty), 0 otherwise
#define Logic_ISBLANK(a) \
        Logic_ISTUPLE(Utils_CAT_(Logic_ISBLANK_,a))

// return 1 if 'a' is 0, 0 otherwise
#define Logic_ISZERO(a) \
        Logic_ISTUPLE(Utils_CAT_(Logic_ISZERO_,a))

// return 1 if 'a' is 1, 0 otherwise
#define Logic_ISONE(a) \
        Logic_ISTUPLE(Utils_CAT_(Logic_ISONE_,a))

// return 1 if 'a' is 2, 0 otherwise
#define Logic_ISTWO(a) \
        Logic_ISTUPLE(Utils_CAT_(Logic_ISTWO_,a))

// return 1 if the token 'a' is a recognized token, 0 otherwise
// the '#define Logic_TOKEN_tok ()' define the set of recognized tokens
#define Logic_ISTOKEN(a) \
        Logic_ISTUPLE(Utils_CAT_(Logic_TOKEN_,a))

// return 1 if its *first* argument is/starts-with a tuple, 0 otherwise
#define Logic_ISTUPLE(...) \
        Utils_PAIR(Utils_ARG1, \
          (Utils_CAT(Logic_ISTUPLE_RET_,Logic_ISTUPLE_TST_ __VA_ARGS__)))

#define Logic_ISNTUPLE(...) \
        Logic_NOT(Logic_ISTUPLE(__VA_ARGS__))

// return 1 if it has no 'effective' argument, 0 otherwise
// WARNING: if the last argument is the name of a function-like macro,
//          the latter will be evaluated
#define Logic_NOARG(...) \
        Logic_AND(Logic_ISTUPLE(__VA_ARGS__ /* Warning macro-eval */ ()), \
                  Logic_NOT(Logic_ISTUPLE(__VA_ARGS__)))

// return 1 if at least one 'effective' argument is present, 0 otherwise
// WARNING: uses Logic_NOARG()
#define Logic_1ARG(...) \
        Logic_NOT(Logic_NOARG(__VA_ARGS__))

// return 1 if at least two arguments are present, 0 otherwise
#define Logic_2ARGS(...) \
        Logic_NOT(Logic_ISONE(Arg_NARG(__VA_ARGS__)))

/***********************************************************
 * Implementation
 */

#define Logic_IF_0(t, ...) __VA_ARGS__
#define Logic_IF_1(t, ...) t

#define Logic_NOT_0() 1
#define Logic_NOT_1() 0

#define Logic_AND_00() 0
#define Logic_AND_01() 0
#define Logic_AND_10() 0
#define Logic_AND_11() 1

#define Logic_OR_00() 0
#define Logic_OR_01() 1
#define Logic_OR_10() 1
#define Logic_OR_11() 1

#define Logic_XOR_00() 0
#define Logic_XOR_01() 1
#define Logic_XOR_10() 1
#define Logic_XOR_11() 0

#define Logic_ISBOOL_0 ()
#define Logic_ISBOOL_1 ()
#define Logic_ISBLANK_ ()
#define Logic_ISZERO_0 ()
#define Logic_ISONE_1  ()
#define Logic_ISTWO_2  ()

#define Logic_ISTUPLE_TST_(...) 1
#define Logic_ISTUPLE_RET_Logic_ISTUPLE_TST_ 0,
#define Logic_ISTUPLE_RET_1 1,

#endif
