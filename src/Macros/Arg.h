#ifndef ARG_H
#define ARG_H

/* 
 * Preprocessor arguments count: adapted from
 * Laurent Deniau http://cern.ch/laurent.deniau/cos.html
 */
 
#define Arg_MAX_N   63

#define Arg_NARG(...) \
        Arg_NARG_(__VA_ARGS__,Arg_REVSEQ_N(),)

#define Arg_DUPSEQ_N(...) \
        Arg_DUPSEQ_N_(__VA_ARGS__)

#define Arg_NUMSEQ_N(...) \
        Arg_NUMSEQ_N_(__VA_ARGS__)

#define Arg_REVSEQ_N(...) \
        Arg_REVSEQ_N_(__VA_ARGS__)



/**
 * Implementation
 */

#define Arg_NARG_(...) \
        Arg_NARG_N_(__VA_ARGS__)

#define Arg_NARG_N_( \
 _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
_11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
_21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
_31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
_41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
_51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
_61,_62,_63,N,...) N

#define Arg_DUPSEQ_N_(...) \
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,__VA_ARGS__,\
__VA_ARGS__,__VA_ARGS__,__VA_ARGS__

#define Arg_NUMSEQ_N_(...) \
__VA_ARGS__## 1,__VA_ARGS__## 2,__VA_ARGS__## 3,__VA_ARGS__## 4,__VA_ARGS__## 5,\
__VA_ARGS__## 6,__VA_ARGS__## 7,__VA_ARGS__## 8,__VA_ARGS__## 9,__VA_ARGS__##10,\
__VA_ARGS__##11,__VA_ARGS__##12,__VA_ARGS__##13,__VA_ARGS__##14,__VA_ARGS__##15,\
__VA_ARGS__##16,__VA_ARGS__##17,__VA_ARGS__##18,__VA_ARGS__##19,__VA_ARGS__##20,\
__VA_ARGS__##21,__VA_ARGS__##22,__VA_ARGS__##23,__VA_ARGS__##24,__VA_ARGS__##25,\
__VA_ARGS__##26,__VA_ARGS__##27,__VA_ARGS__##28,__VA_ARGS__##29,__VA_ARGS__##30,\
__VA_ARGS__##31,__VA_ARGS__##32,__VA_ARGS__##33,__VA_ARGS__##34,__VA_ARGS__##35,\
__VA_ARGS__##36,__VA_ARGS__##37,__VA_ARGS__##38,__VA_ARGS__##39,__VA_ARGS__##40,\
__VA_ARGS__##41,__VA_ARGS__##42,__VA_ARGS__##43,__VA_ARGS__##44,__VA_ARGS__##45,\
__VA_ARGS__##46,__VA_ARGS__##47,__VA_ARGS__##48,__VA_ARGS__##49,__VA_ARGS__##50,\
__VA_ARGS__##51,__VA_ARGS__##52,__VA_ARGS__##53,__VA_ARGS__##54,__VA_ARGS__##55,\
__VA_ARGS__##56,__VA_ARGS__##57,__VA_ARGS__##58,__VA_ARGS__##59,__VA_ARGS__##60,\
__VA_ARGS__##61,__VA_ARGS__##62,__VA_ARGS__##63

#define Arg_REVSEQ_N_(...) \
__VA_ARGS__##63,__VA_ARGS__##62,__VA_ARGS__##61,\
__VA_ARGS__##60,__VA_ARGS__##59,__VA_ARGS__##58,__VA_ARGS__##57,__VA_ARGS__##56,\
__VA_ARGS__##55,__VA_ARGS__##54,__VA_ARGS__##53,__VA_ARGS__##52,__VA_ARGS__##51,\
__VA_ARGS__##50,__VA_ARGS__##49,__VA_ARGS__##48,__VA_ARGS__##47,__VA_ARGS__##46,\
__VA_ARGS__##45,__VA_ARGS__##44,__VA_ARGS__##43,__VA_ARGS__##42,__VA_ARGS__##41,\
__VA_ARGS__##40,__VA_ARGS__##39,__VA_ARGS__##38,__VA_ARGS__##37,__VA_ARGS__##36,\
__VA_ARGS__##35,__VA_ARGS__##34,__VA_ARGS__##33,__VA_ARGS__##32,__VA_ARGS__##31,\
__VA_ARGS__##30,__VA_ARGS__##29,__VA_ARGS__##28,__VA_ARGS__##27,__VA_ARGS__##26,\
__VA_ARGS__##25,__VA_ARGS__##24,__VA_ARGS__##23,__VA_ARGS__##22,__VA_ARGS__##21,\
__VA_ARGS__##20,__VA_ARGS__##19,__VA_ARGS__##18,__VA_ARGS__##17,__VA_ARGS__##16,\
__VA_ARGS__##15,__VA_ARGS__##14,__VA_ARGS__##13,__VA_ARGS__##12,__VA_ARGS__##11,\
__VA_ARGS__##10,__VA_ARGS__## 9,__VA_ARGS__## 8,__VA_ARGS__## 7,__VA_ARGS__## 6,\
__VA_ARGS__## 5,__VA_ARGS__## 4,__VA_ARGS__## 3,__VA_ARGS__## 2,__VA_ARGS__## 1



/* Example
 * 
#define VFUNC(fct,...) \
        Utils_CAT_NARG(fct_,__VA_ARGS__)(__VA_ARGS__)



// definition for FOO
#define FOO(...) VFUNC(FOO, __VA_ARGS__)
// Implementation
#define FOO_2(x, y) ((x) + (y))
#define FOO_3(x, y, z) ((x) + (y) + (z))

// it also works with C functions:
int FOO_4(int a, int b, int c, int d) { return a + b + c + d; }

// Examples
FOO(42, 42) // will use macro function FOO_2
FOO(42, 42, 42) // will use macro function FOO_3
FOO(42, 42, 42, 42) // will call FOO_4 function
*/

/*
Limitations

    Only up to 63 arguments (but expandable)
    Function for no argument only in GCC possible
*/


#endif
