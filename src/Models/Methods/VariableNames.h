#ifndef VARIABLENAMES_H
#define VARIABLENAMES_H

/* Inputs given as a tuple T where each element of T is 
 * a tuple composed of 2 elements, the name and its length: 
 * (N,L) where L = 1,3,9 so T = ((N1,L1),(N2,L2),...) */

#include "Tuple.h"
#include "Arg.h"

#define VariableNames_NbOfVariables(VN) \
        Tuple_LEN(VN)

#define VariableNames_TotalLength(VN) \


#define VariableName_Length(V) \
        Tuple_2ND(V)
        



#endif
