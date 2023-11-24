#ifndef LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_H
#define LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_H


extern void Log10DissociationConstantOfCementHydrationProduct_Print(double) ;



#define Log10DissociationConstantOfCementHydrationProduct(R,T) \
        (Log10DissociationConstantOfCementHydrationProduct_##R(T))


#include "Log10DissociationConstantOfCementHydrationProduct_DEFAULT.h.in"
 

#endif
