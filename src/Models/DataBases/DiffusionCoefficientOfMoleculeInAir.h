#ifndef DIFFUSIONCOEFFICIENTOFMOLECULEINAIR_H
#define DIFFUSIONCOEFFICIENTOFMOLECULEINAIR_H



#include "InternationalSystemOfUnits.h"

#define DiffusionCoefficientOfMoleculeInAir_Unit \
       (InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter/InternationalSystemOfUnits_OneSecond)



#include "AirViscosity.h"

/* Air viscosity at T = 291.15 (Pa.s) and P = 1 atm */
#define AirViscosity291      AirViscosity(291.15)     //(1.827e-5)

#define DiffusionCoefficientOfMoleculeInAir(I,T) \
       (DiffusionCoefficientOfMoleculeInAir_##I*DiffusionCoefficientOfMoleculeInAir_Unit*((T*AirViscosity(291.15))/(291.15*AirViscosity(T))))



/*
 * Molecular Diffusion Coefficient in Air At Temperature T = 291.15 K (m2/s)
 */
#define DiffusionCoefficientOfMoleculeInAir_CO2 \
       (1.6e-5)

#define DiffusionCoefficientOfMoleculeInAir_H2O  \
      ((2.82e-5)*((291.15*AirViscosity(298))/(298*AirViscosity(291.15))))


/* References
 * 
 * Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems (2nd ed.). New York: Cambridge University Press. ISBN 0-521-45078-0.
 * 
 */

#endif
