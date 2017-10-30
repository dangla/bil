#ifndef AIRVISCOSITY_H
#define AIRVISCOSITY_H

#include "InternationalSystemOfUnits.h"

#define AirViscosity_Unit \
       (InternationalSystemOfUnits_OnePascal*InternationalSystemOfUnits_OneSecond)

/* Air viscosity as function of T (Pa.s) and P (Pa) */
#define AirViscosity_T291 \
       (1.827e-5*(AirViscosity_Unit))


/* Sutherland formula: mu(T) = lambda*T^1.5/(T + C)
 * with lambda = mu0*(T0 + C)/T0^1.5 = 1.512041288e-6 
 * for mu0 = 1.827e-5, T0 = 291.15, C = 120
 * Validity range: 0 < T < 555 K
 */

#define AirViscosity(T) \
       (1.512041288e-6*(T)*sqrt(T)/((T) + 120)*(AirViscosity_Unit))


#endif
