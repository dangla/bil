#ifndef WATERVISCOSITY_H
#define WATERVISCOSITY_H

#include "InternationalSystemOfUnits.h"

#define WaterViscosity_Unit \
       (InternationalSystemOfUnits_OnePascal*InternationalSystemOfUnits_OneSecond)

/* Water viscosity as function of T (Pa.s) */
/* www.engineeringtoolbox.com, equation from chemistry handbook */
#define WaterViscosity_T293  (1.002e-3*WaterViscosity_Unit)

#include <math.h>

#define WaterViscosity(T)    (WaterViscosity_T293*pow(10.,(-1.3272*((T)-293.) - 0.001053*((T)-293.)*((T)-293.))/((T)-168.)))

#endif
