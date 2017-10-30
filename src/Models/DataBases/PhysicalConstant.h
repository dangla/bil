#ifndef PHYSICALCONSTANT_H
#define PHYSICALCONSTANT_H

#include "InternationalSystemOfUnits.h"

#define PhysicalConstant_FaradayUnit \
       (InternationalSystemOfUnits_OneCoulomb/InternationalSystemOfUnits_OneMole)
       
#define PhysicalConstant_PerfectGasConstantUnit \
       (InternationalSystemOfUnits_OneJoule/(InternationalSystemOfUnits_OneMole*InternationalSystemOfUnits_OneKelvin))
       
#define PhysicalConstant_BoltzmannConstantUnit \
       (InternationalSystemOfUnits_OneJoule/(InternationalSystemOfUnits_OneKelvin))



#define PhysicalConstant(I)       PhysicalConstant_##I

/* Physical constants (Units SI: Length = m, Mass = kg, Time = s, Degree = K) */

#define PhysicalConstant_Faraday \
       (9.64846e4)*PhysicalConstant_FaradayUnit            /* (C/mol) */
       
#define PhysicalConstant_PerfectGasConstant \
       (8.3143)*PhysicalConstant_PerfectGasConstantUnit    /* (J/mol/K) */
       
#define PhysicalConstant_BoltzmannConstant \
       (1.38e-23)*PhysicalConstant_BoltzmannConstantUnit   /* (J/K) */

#endif
