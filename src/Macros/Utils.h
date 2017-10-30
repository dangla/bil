#ifndef UTILS_H
#define UTILS_H

/* 
 * Preprocessor utility macros: adapted from
 * Laurent Deniau http://cern.ch/laurent.deniau/cos.html
 */
 
#include "Arg.h"


#define Utils_ID(...)       __VA_ARGS__
#define Utils_EAT(...)      /* nothing */
#define Utils_PART(...)     __VA_ARGS__,
#define Utils_PPART(...)   (__VA_ARGS__),
#define Utils_PAIR(a,...) a __VA_ARGS__
#define Utils_SWAP(a,...)   __VA_ARGS__ a
#define Utils_REST(a,...)   __VA_ARGS__
#define Utils_EMPTY()       /* empty */
#define Utils_COMMA()       ,
#define Utils_LPAR()        (
#define Utils_RPAR()        )

#define Utils_ARG1(a,                ...)  a
#define Utils_ARG2(a,b,              ...)  b
#define Utils_ARG3(a,b,c,            ...)  c
#define Utils_ARG4(a,b,c,d,          ...)  d
#define Utils_ARG5(a,b,c,d,e,        ...)  e
#define Utils_ARG6(a,b,c,d,e,f,      ...)  f
#define Utils_ARG7(a,b,c,d,e,f,g,    ...)  g
#define Utils_ARG8(a,b,c,d,e,f,g,h,  ...)  h
#define Utils_ARG9(a,b,c,d,e,f,g,h,i,...)  i

#define Utils_STR(...)        Utils_STR_ (      __VA_ARGS__)
#define Utils_CAT(a,...)      Utils_CAT_ (a,    __VA_ARGS__)
#define Utils_CAT3(a,b,...)   Utils_CAT3_(a,b,  __VA_ARGS__)
#define Utils_CAT4(a,b,c,...) Utils_CAT4_(a,b,c,__VA_ARGS__)
#define Utils_CAT_NARG(a,...) Utils_CAT  (a,Arg_NARG(__VA_ARGS__))

#define Utils_DUP(n,...)      Utils_CAT_(Utils_DUP_,n)(__VA_ARGS__)



/**
 * Implementation
 */

#define Utils_STR_(...)        #__VA_ARGS__
#define Utils_CAT_(a,...)      a##__VA_ARGS__
#define Utils_CAT3_(a,b,...)   a##b##__VA_ARGS__
#define Utils_CAT4_(a,b,c,...) a##b##c##__VA_ARGS__

#define Utils_DUP_0( ...)
#define Utils_DUP_1( ...) __VA_ARGS__
#define Utils_DUP_2( ...) Utils_DUP_1 (__VA_ARGS__)Utils_DUP_1 (__VA_ARGS__)
#define Utils_DUP_3( ...) Utils_DUP_1 (__VA_ARGS__)Utils_DUP_1 (__VA_ARGS__)Utils_DUP_1 (__VA_ARGS__)
#define Utils_DUP_4( ...) Utils_DUP_1 (__VA_ARGS__)Utils_DUP_1 (__VA_ARGS__)Utils_DUP_2 (__VA_ARGS__)
#define Utils_DUP_5( ...) Utils_DUP_1 (__VA_ARGS__)Utils_DUP_2 (__VA_ARGS__)Utils_DUP_2 (__VA_ARGS__)
#define Utils_DUP_6( ...) Utils_DUP_2 (__VA_ARGS__)Utils_DUP_2 (__VA_ARGS__)Utils_DUP_2 (__VA_ARGS__)
#define Utils_DUP_7( ...) Utils_DUP_2 (__VA_ARGS__)Utils_DUP_2 (__VA_ARGS__)Utils_DUP_3 (__VA_ARGS__)
#define Utils_DUP_8( ...) Utils_DUP_2 (__VA_ARGS__)Utils_DUP_3 (__VA_ARGS__)Utils_DUP_3 (__VA_ARGS__)
#define Utils_DUP_9( ...) Utils_DUP_3 (__VA_ARGS__)Utils_DUP_3 (__VA_ARGS__)Utils_DUP_3 (__VA_ARGS__)
#define Utils_DUP_10(...) Utils_DUP_3 (__VA_ARGS__)Utils_DUP_3 (__VA_ARGS__)Utils_DUP_4 (__VA_ARGS__)
#define Utils_DUP_11(...) Utils_DUP_3 (__VA_ARGS__)Utils_DUP_4 (__VA_ARGS__)Utils_DUP_4 (__VA_ARGS__)
#define Utils_DUP_12(...) Utils_DUP_4 (__VA_ARGS__)Utils_DUP_4 (__VA_ARGS__)Utils_DUP_4 (__VA_ARGS__)
#define Utils_DUP_13(...) Utils_DUP_4 (__VA_ARGS__)Utils_DUP_4 (__VA_ARGS__)Utils_DUP_5 (__VA_ARGS__)
#define Utils_DUP_14(...) Utils_DUP_4 (__VA_ARGS__)Utils_DUP_5 (__VA_ARGS__)Utils_DUP_5 (__VA_ARGS__)
#define Utils_DUP_15(...) Utils_DUP_5 (__VA_ARGS__)Utils_DUP_5 (__VA_ARGS__)Utils_DUP_5 (__VA_ARGS__)
#define Utils_DUP_16(...) Utils_DUP_5 (__VA_ARGS__)Utils_DUP_5 (__VA_ARGS__)Utils_DUP_6 (__VA_ARGS__)
#define Utils_DUP_17(...) Utils_DUP_5 (__VA_ARGS__)Utils_DUP_6 (__VA_ARGS__)Utils_DUP_6 (__VA_ARGS__)
#define Utils_DUP_18(...) Utils_DUP_6 (__VA_ARGS__)Utils_DUP_6 (__VA_ARGS__)Utils_DUP_6 (__VA_ARGS__)
#define Utils_DUP_19(...) Utils_DUP_6 (__VA_ARGS__)Utils_DUP_6 (__VA_ARGS__)Utils_DUP_7 (__VA_ARGS__)
#define Utils_DUP_20(...) Utils_DUP_6 (__VA_ARGS__)Utils_DUP_7 (__VA_ARGS__)Utils_DUP_7 (__VA_ARGS__)
#define Utils_DUP_21(...) Utils_DUP_7 (__VA_ARGS__)Utils_DUP_7 (__VA_ARGS__)Utils_DUP_7 (__VA_ARGS__)
#define Utils_DUP_22(...) Utils_DUP_7 (__VA_ARGS__)Utils_DUP_7 (__VA_ARGS__)Utils_DUP_8 (__VA_ARGS__)
#define Utils_DUP_23(...) Utils_DUP_7 (__VA_ARGS__)Utils_DUP_8 (__VA_ARGS__)Utils_DUP_8 (__VA_ARGS__)
#define Utils_DUP_24(...) Utils_DUP_8 (__VA_ARGS__)Utils_DUP_8 (__VA_ARGS__)Utils_DUP_8 (__VA_ARGS__)
#define Utils_DUP_25(...) Utils_DUP_8 (__VA_ARGS__)Utils_DUP_8 (__VA_ARGS__)Utils_DUP_9 (__VA_ARGS__)
#define Utils_DUP_26(...) Utils_DUP_8 (__VA_ARGS__)Utils_DUP_9 (__VA_ARGS__)Utils_DUP_9 (__VA_ARGS__)
#define Utils_DUP_27(...) Utils_DUP_9 (__VA_ARGS__)Utils_DUP_9 (__VA_ARGS__)Utils_DUP_9 (__VA_ARGS__)
#define Utils_DUP_28(...) Utils_DUP_9 (__VA_ARGS__)Utils_DUP_9 (__VA_ARGS__)Utils_DUP_10(__VA_ARGS__)
#define Utils_DUP_29(...) Utils_DUP_9 (__VA_ARGS__)Utils_DUP_10(__VA_ARGS__)Utils_DUP_10(__VA_ARGS__)
#define Utils_DUP_30(...) Utils_DUP_10(__VA_ARGS__)Utils_DUP_10(__VA_ARGS__)Utils_DUP_10(__VA_ARGS__)
#define Utils_DUP_31(...) Utils_DUP_10(__VA_ARGS__)Utils_DUP_10(__VA_ARGS__)Utils_DUP_11(__VA_ARGS__)
#define Utils_DUP_32(...) Utils_DUP_10(__VA_ARGS__)Utils_DUP_11(__VA_ARGS__)Utils_DUP_11(__VA_ARGS__)
#define Utils_DUP_33(...) Utils_DUP_11(__VA_ARGS__)Utils_DUP_11(__VA_ARGS__)Utils_DUP_11(__VA_ARGS__)
#define Utils_DUP_34(...) Utils_DUP_11(__VA_ARGS__)Utils_DUP_11(__VA_ARGS__)Utils_DUP_12(__VA_ARGS__)
#define Utils_DUP_35(...) Utils_DUP_11(__VA_ARGS__)Utils_DUP_12(__VA_ARGS__)Utils_DUP_12(__VA_ARGS__)
#define Utils_DUP_36(...) Utils_DUP_12(__VA_ARGS__)Utils_DUP_12(__VA_ARGS__)Utils_DUP_12(__VA_ARGS__)
#define Utils_DUP_37(...) Utils_DUP_12(__VA_ARGS__)Utils_DUP_12(__VA_ARGS__)Utils_DUP_13(__VA_ARGS__)
#define Utils_DUP_38(...) Utils_DUP_12(__VA_ARGS__)Utils_DUP_13(__VA_ARGS__)Utils_DUP_13(__VA_ARGS__)
#define Utils_DUP_39(...) Utils_DUP_13(__VA_ARGS__)Utils_DUP_13(__VA_ARGS__)Utils_DUP_13(__VA_ARGS__)
#define Utils_DUP_40(...) Utils_DUP_13(__VA_ARGS__)Utils_DUP_13(__VA_ARGS__)Utils_DUP_14(__VA_ARGS__)
#define Utils_DUP_41(...) Utils_DUP_13(__VA_ARGS__)Utils_DUP_14(__VA_ARGS__)Utils_DUP_14(__VA_ARGS__)
#define Utils_DUP_42(...) Utils_DUP_14(__VA_ARGS__)Utils_DUP_14(__VA_ARGS__)Utils_DUP_14(__VA_ARGS__)
#define Utils_DUP_43(...) Utils_DUP_14(__VA_ARGS__)Utils_DUP_14(__VA_ARGS__)Utils_DUP_15(__VA_ARGS__)
#define Utils_DUP_44(...) Utils_DUP_14(__VA_ARGS__)Utils_DUP_15(__VA_ARGS__)Utils_DUP_15(__VA_ARGS__)
#define Utils_DUP_45(...) Utils_DUP_15(__VA_ARGS__)Utils_DUP_15(__VA_ARGS__)Utils_DUP_15(__VA_ARGS__)
#define Utils_DUP_46(...) Utils_DUP_15(__VA_ARGS__)Utils_DUP_15(__VA_ARGS__)Utils_DUP_16(__VA_ARGS__)
#define Utils_DUP_47(...) Utils_DUP_15(__VA_ARGS__)Utils_DUP_16(__VA_ARGS__)Utils_DUP_16(__VA_ARGS__)
#define Utils_DUP_48(...) Utils_DUP_16(__VA_ARGS__)Utils_DUP_16(__VA_ARGS__)Utils_DUP_16(__VA_ARGS__)
#define Utils_DUP_49(...) Utils_DUP_16(__VA_ARGS__)Utils_DUP_16(__VA_ARGS__)Utils_DUP_17(__VA_ARGS__)
#define Utils_DUP_50(...) Utils_DUP_16(__VA_ARGS__)Utils_DUP_17(__VA_ARGS__)Utils_DUP_17(__VA_ARGS__)
#define Utils_DUP_51(...) Utils_DUP_17(__VA_ARGS__)Utils_DUP_17(__VA_ARGS__)Utils_DUP_17(__VA_ARGS__)
#define Utils_DUP_52(...) Utils_DUP_17(__VA_ARGS__)Utils_DUP_17(__VA_ARGS__)Utils_DUP_18(__VA_ARGS__)
#define Utils_DUP_53(...) Utils_DUP_17(__VA_ARGS__)Utils_DUP_18(__VA_ARGS__)Utils_DUP_18(__VA_ARGS__)
#define Utils_DUP_54(...) Utils_DUP_18(__VA_ARGS__)Utils_DUP_18(__VA_ARGS__)Utils_DUP_18(__VA_ARGS__)
#define Utils_DUP_55(...) Utils_DUP_18(__VA_ARGS__)Utils_DUP_18(__VA_ARGS__)Utils_DUP_19(__VA_ARGS__)
#define Utils_DUP_56(...) Utils_DUP_18(__VA_ARGS__)Utils_DUP_19(__VA_ARGS__)Utils_DUP_19(__VA_ARGS__)
#define Utils_DUP_57(...) Utils_DUP_19(__VA_ARGS__)Utils_DUP_19(__VA_ARGS__)Utils_DUP_19(__VA_ARGS__)
#define Utils_DUP_58(...) Utils_DUP_19(__VA_ARGS__)Utils_DUP_19(__VA_ARGS__)Utils_DUP_20(__VA_ARGS__)
#define Utils_DUP_59(...) Utils_DUP_19(__VA_ARGS__)Utils_DUP_20(__VA_ARGS__)Utils_DUP_20(__VA_ARGS__)
#define Utils_DUP_60(...) Utils_DUP_20(__VA_ARGS__)Utils_DUP_20(__VA_ARGS__)Utils_DUP_20(__VA_ARGS__)
#define Utils_DUP_61(...) Utils_DUP_20(__VA_ARGS__)Utils_DUP_20(__VA_ARGS__)Utils_DUP_21(__VA_ARGS__)
#define Utils_DUP_62(...) Utils_DUP_20(__VA_ARGS__)Utils_DUP_21(__VA_ARGS__)Utils_DUP_21(__VA_ARGS__)
#define Utils_DUP_63(...) Utils_DUP_21(__VA_ARGS__)Utils_DUP_21(__VA_ARGS__)Utils_DUP_21(__VA_ARGS__)
#define Utils_DUP_N( ...) Utils_DUP_63(__VA_ARGS__)

#endif
