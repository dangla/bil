#ifndef LOG10EQUILIBRIUMCONSTANTOFHOMOGENEOUSREACTIONINWATER_H
#define LOG10EQUILIBRIUMCONSTANTOFHOMOGENEOUSREACTIONINWATER_H


extern void Log10EquilibriumConstantOfHomogeneousReactionInWater_Print(double) ;


#define Log10EquilibriumConstantOfHomogeneousReactionInWater(R,T) \
        (Log10EquilibriumConstantOfHomogeneousReactionInWater_##R(T))


#include "Log10EquilibriumConstantOfHomogeneousReactionInWater_DEFAULT.h"


#endif
