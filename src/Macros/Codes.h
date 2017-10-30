#ifndef CODES_H
#define CODES_H

/*
 * Preprocessor coding/decoding macros: adapted from 
 * http://cern.ch/laurent.deniau/cos.html
 */

#include "Utils.h"
#include "Tuple.h"
#include "Logic.h"

// predicates
#define Codes_ISCODE(a) \
        Logic_ISTUPLE(Utils_CAT_(Codes_CODE_,a))

#define Codes_ISDIGIT(a) \
        Tuple_1ST(Utils_CAT_(Codes_CODE_,a))

#define Codes_ISXDIGIT(a) \
        Tuple_2ND(Utils_CAT_(Codes_CODE_,a))

#define Codes_ISUPPER(a) \
        Tuple_3RD(Utils_CAT_(Codes_CODE_,a))

#define Codes_ISLOWER(a) \
        Tuple_4TH(Utils_CAT_(Codes_CODE_,a))

#define Codes_ISALPHA(a) \
        Logic_OR(Codes_ISUPPER(a),Codes_ISLOWER(a))

#define Codes_ISALNUM(a) \
        Logic_OR(Codes_ISALPHA(a),Codes_ISDIGIT(a))

#define Codes_ISUNDERSCORE(a) \
        Logic_EQ(Codes_ORD(a),62)

#define Codes_ISSPACE(a) \
        Logic_EQ(Codes_ORD(a),63)

// conversions
#define Codes_ORD(a) \
        Tuple_5TH(Utils_CAT_(Codes_CODE_,a))

#define Codes_TOUPPER(a) \
        Tuple_6TH(Utils_CAT_(Codes_CODE_,a))

#define Codes_TOLOWER(a) \
        Tuple_7TH(Utils_CAT_(Codes_CODE_,a))

#define Codes_SPELL(a) \
        Tuple_8TH(Utils_CAT_(Codes_CODE_,a))

#define Codes_MORSE(a) \
        Tuple_9TH(Utils_CAT_(Codes_CODE_,a))

/***********************************************************
 * Implementation
 */

/* NOTE-DEV: data structure
   isdigit, isxdigit, isupper, islower, ord, upper, lower, NATO, Morse
 */
#define Codes_CODE_0 (1, 1, 0, 0,  0, 0, 0, Zero    , .____, )
#define Codes_CODE_1 (1, 1, 0, 0,  1, 1, 1, One     , ..___, )
#define Codes_CODE_2 (1, 1, 0, 0,  2, 2, 2, Two     , ...__, )
#define Codes_CODE_3 (1, 1, 0, 0,  3, 3, 3, Three   , ...._, )
#define Codes_CODE_4 (1, 1, 0, 0,  4, 4, 4, Four    , ....., )
#define Codes_CODE_5 (1, 1, 0, 0,  5, 5, 5, Five    , _...., )
#define Codes_CODE_6 (1, 1, 0, 0,  6, 6, 6, Six     , __..., )
#define Codes_CODE_7 (1, 1, 0, 0,  7, 7, 7, Seven   , ___.., )
#define Codes_CODE_8 (1, 1, 0, 0,  8, 8, 8, Eight   , ____., )
#define Codes_CODE_9 (1, 1, 0, 0,  9, 9, 9, Nine    , _____, )
#define Codes_CODE_A (0, 1, 1, 0, 10, A, a, Alfa    , ._   , )
#define Codes_CODE_B (0, 1, 1, 0, 11, B, b, Bravo   , _... , )
#define Codes_CODE_C (0, 1, 1, 0, 12, C, c, Charlie , _._. , )
#define Codes_CODE_D (0, 1, 1, 0, 13, D, d, Delta   , _..  , )
#define Codes_CODE_E (0, 1, 1, 0, 14, E, e, Echo    , .    , )
#define Codes_CODE_F (0, 1, 1, 0, 15, F, f, Foxtrot , .._. , )
#define Codes_CODE_G (0, 0, 1, 0, 16, G, g, Golf    , __.  , )
#define Codes_CODE_H (0, 0, 1, 0, 17, H, h, Hotel   , .... , )
#define Codes_CODE_I (0, 0, 1, 0, 18, I, i, India   , ..   , )
#define Codes_CODE_J (0, 0, 1, 0, 19, J, j, Juliett , .___ , )
#define Codes_CODE_K (0, 0, 1, 0, 20, K, k, Kilo    , _._  , )
#define Codes_CODE_L (0, 0, 1, 0, 21, L, l, Lima    , ._.. , )
#define Codes_CODE_M (0, 0, 1, 0, 22, M, m, Mike    , __   , )
#define Codes_CODE_N (0, 0, 1, 0, 23, N, n, November, _.   , )
#define Codes_CODE_O (0, 0, 1, 0, 24, O, o, Oscar   , ___  , )
#define Codes_CODE_P (0, 0, 1, 0, 25, P, p, Papa    , .__. , )
#define Codes_CODE_Q (0, 0, 1, 0, 26, Q, q, Quebec  , __._ , )
#define Codes_CODE_R (0, 0, 1, 0, 27, R, r, Romeo   , ._.  , )
#define Codes_CODE_S (0, 0, 1, 0, 28, S, s, Sierra  , ...  , )
#define Codes_CODE_T (0, 0, 1, 0, 29, T, t, Tango   , _    , )
#define Codes_CODE_U (0, 0, 1, 0, 30, U, u, Uniform , .._  , )
#define Codes_CODE_V (0, 0, 1, 0, 31, V, v, Victor  , ..._ , )
#define Codes_CODE_W (0, 0, 1, 0, 32, W, w, Whiskey , .__  , )
#define Codes_CODE_X (0, 0, 1, 0, 33, X, x, Xray    , _.._ , )
#define Codes_CODE_Y (0, 0, 1, 0, 34, Y, y, Yankee  , _.__ , )
#define Codes_CODE_Z (0, 0, 1, 0, 35, Z, z, Zulu    , __.. , )
#define Codes_CODE_a (0, 1, 0, 1, 36, A, a, Alfa    , ._   , )
#define Codes_CODE_b (0, 1, 0, 1, 37, B, b, Bravo   , _... , )
#define Codes_CODE_c (0, 1, 0, 1, 38, C, c, Charlie , _._. , )
#define Codes_CODE_d (0, 1, 0, 1, 39, D, d, Delta   , _..  , )
#define Codes_CODE_e (0, 1, 0, 1, 40, E, e, Echo    , .    , )
#define Codes_CODE_f (0, 1, 0, 1, 41, F, f, Foxtrot , .._. , )
#define Codes_CODE_g (0, 0, 0, 1, 42, G, g, Golf    , __.  , )
#define Codes_CODE_h (0, 0, 0, 1, 43, H, h, Hotel   , .... , )
#define Codes_CODE_i (0, 0, 0, 1, 44, I, i, India   , ..   , )
#define Codes_CODE_j (0, 0, 0, 1, 45, J, j, Juliett , .___ , )
#define Codes_CODE_k (0, 0, 0, 1, 46, K, k, Kilo    , _._  , )
#define Codes_CODE_l (0, 0, 0, 1, 47, L, l, Lima    , ._.. , )
#define Codes_CODE_m (0, 0, 0, 1, 48, M, m, Mike    , __   , )
#define Codes_CODE_n (0, 0, 0, 1, 49, N, n, November, _.   , )
#define Codes_CODE_o (0, 0, 0, 1, 50, O, o, Oscar   , ___  , )
#define Codes_CODE_p (0, 0, 0, 1, 51, P, p, Papa    , .__. , )
#define Codes_CODE_q (0, 0, 0, 1, 52, Q, q, Quebec  , __._ , )
#define Codes_CODE_r (0, 0, 0, 1, 53, R, r, Romeo   , ._.  , )
#define Codes_CODE_s (0, 0, 0, 1, 54, S, s, Sierra  , ...  , )
#define Codes_CODE_t (0, 0, 0, 1, 55, T, t, Tango   , _    , )
#define Codes_CODE_u (0, 0, 0, 1, 56, U, u, Uniform , .._  , )
#define Codes_CODE_v (0, 0, 0, 1, 57, V, v, Victor  , ..._ , )
#define Codes_CODE_w (0, 0, 0, 1, 58, W, w, Whiskey , .__  , )
#define Codes_CODE_x (0, 0, 0, 1, 59, X, x, Xray    , _.._ , )
#define Codes_CODE_y (0, 0, 0, 1, 60, Y, y, Yankee  , _.__ , )
#define Codes_CODE_z (0, 0, 0, 1, 61, Z, z, Zulu    , __.. , )
#define Codes_CODE__ (0, 0, 0, 0, 62, _, _, Stop    , /    , )
#define Codes_CODE_  (0, 0, 0, 0, 63,  ,  , Space   ,      , )

#endif
