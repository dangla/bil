#ifndef IONIZATIONCONSTANTOFWATER_H
#define IONIZATIONCONSTANTOFWATER_H

/* The International Association for the Properties of Water and Steam
 * Lucerne, Switzerland, August 2007
 * Release on the Ionization Constant of H2O
 * Analytical equation for pK_w over ranges of 
 * water density from 0 to 1.25 g/cm3
 * and temperature from 0 to 800 °C.
 * 
 * pK_w = -2*n*[log(1 + Q) - Q/(1 + Q)*rho*(b0 + b1/T +  b2*rho)] + pK_w^G
 * 
 * with Q = (rho/rho_0)*exp(a0 + a1/T + a2/T^2*rho^(2/3))
 * 
 * Units:  rho given in (g/cm3) and T  given in (K)
 * 
 * At rho = 1 g/cm3
 * 
 * pK_w = -2*n*[log(1 + Q) - Q/(1 + Q)*(b3 + b1/T)] + pK_w^G
 * 
 * with Q = exp(a0 + a1/T + a2/T^2)
*/

#define IonizationConstantOfWater(T)  (pow(10,-pIonizationConstantOfWater(T)))

#define pIonizationConstantOfWater_G(T) (-2.87530499 + (4.825133e4 + (-6.770793e4 + 1.010210e7/(T))/(T))/(T))

#define pIonizationConstantOfWater(T)   (pIonizationConstantOfWater_L1(T) + \
                                         pIonizationConstantOfWater_L2(T) + \
                                         pIonizationConstantOfWater_G(T))

#define pIonizationConstantOfWater_L0(T) (exp(-0.864671 + (8659.19 - 22786.2/(T))/(T)))

#define pIonizationConstantOfWater_L1(T) (-12*log10(1 + pIonizationConstantOfWater_L0(T)))

#define pIonizationConstantOfWater_L2(T) ((1 - 1/(1 + pIonizationConstantOfWater_L0(T)))*(3.19548 - 682.2408/(T)))

#endif
