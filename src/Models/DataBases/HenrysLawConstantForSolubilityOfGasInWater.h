#ifndef HENRYSLAWCONSTANTFORSOLUBILITYOFGASINWATER_H
#define HENRYSLAWCONSTANTFORSOLUBILITYOFGASINWATER_H

#include "InternationalSystemOfUnits.h"

#define HenrysLawConstantForSolubilityOfGasInWater_Unit \
       (InternationalSystemOfUnits_OneMole/(InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OnePascal))

/*
 * Henry's law constant is defined usually by
 * 
 *   k_H = c_a / p_g
 * 
 *   c_a = concentration of the species in the aqueous phase
 *   p_g = partial pressure of that species in the gas phase
 *   k_H refers to standard conditions T0 = 298.15
 * 
 * Henry's law constant can also be expressed as the dimensionless ratio
 * between the aqueous phase concentration c_a of a species and its gas
 * phase concentration c_g:
 * 
 *   kc_H = c_a / c_g = k_H * RT
 * 
 * A simple way to describe Henry's law as a function of temperature is
 * 
 *   k_H = k_H0 * exp( A * (1/T - 1/T0) )
 * 
 * where A = - DH/R with DH the enthalpy of solution.
 * 
 * Below the k_H definition is used.
 * 
 * The official SI unit is mol/m3/Pa. The common used unit is mol/L/bar.
 * Conversions: 
 *              1  mol/(m3*Pa) = 101.325  mol/(L*bar)
 *              1  mol/(L*bar) = 9.87e-3  mol/(m3*Pa)
 */


#define T298   (298.)

#define HenrysLawConstantForSolubilityOfGasInWater(I,T) \
       (HenrysLawConstantForSolubilityOfGasInWater_##I(T)*HenrysLawConstantForSolubilityOfGasInWater_Unit)


/*
 * Henry's law constant for the solubility of gas in water ( mol/(m3*Pa) )
 */
/* Oxygen compounds */
#define HenrysLawConstantForSolubilityOfGasInWater_O2(T)   (1.28e-5*exp(1500*(1/T - 1/T298)))

/* Carbon compounds */
#define HenrysLawConstantForSolubilityOfGasInWater_CO2(T)  (3.45e-4*exp(2400*(1/T - 1/T298)))

/* Sulfur compounds */
#define HenrysLawConstantForSolubilityOfGasInWater_H2S(T)  (9.87e-4*exp(2100*(1/T - 1/T298)))


/* References:
 * Sander R., Compilation of Henry's law constants, version 3.99,
 * Atmospheric Chemistry and Physics Discussions, 
 * 14 (21), p. 29615-30521, 2014. doi: 10.5194/acpd-14-29615-2014.
 * (1) R. Sander, Compilation of Henry's Law Constants  for Inorganic and
 * Organic Species of Potential Importance in Environmental Chemistry
 * (http://www.mpch-mainz.mpg.de/~sander/res/henry.html).
 */
#endif
