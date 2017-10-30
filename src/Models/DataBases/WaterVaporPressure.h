#ifndef WATERVAPORPRESSURE_H
#define WATERVAPORPRESSURE_H

#include "InternationalSystemOfUnits.h"

#define WaterVaporPressure_Unit \
       (InternationalSystemOfUnits_OnePascal)


#define WaterVaporPressure(T)  WaterVaporPressure_Antoine(T)


/* Validity range: -80:50 °C */
#define WaterVaporPressure_Buck(T) \
        ((611.21*exp((19.842819 - (T)/234.5)*(((T)-273.15)/(-16.01+(T)))))*InternationalSystemOfUnits_OnePascal)

/* Validity range: -50:102 °C */
#define WaterVaporPressure_GoffGratch(T) \
        (((101324.6)*pow(10,-7.90298*(373.16/(T) - 1)) * pow(373.16/(T),5.02808) * pow(10,-1.3816e-7*(pow(10,11.344*(1 - (T)/373.16)) - 1) + 8.1328e-3*(pow(10,-3.49149*(373.16/(T)  - 1)) - 1)))*InternationalSystemOfUnits_OnePascal)

/* Validity range: 1:100 °C */
#define WaterVaporPressure_Antoine(T) \
        ((pow(10,10.196213 - 1730.63/(-39.724 + (T))))*InternationalSystemOfUnits_OnePascal)


/* Validity range: 0:360 °C */
#define WaterVaporPressure_Wagner(T) \
        ((2.2064e7*exp((-7.85951783*(1 - (T)/647.096) + 1.84408259*pow((1 - (T)/647.096),1.5) - 11.7866497*pow((1 - (T)/647.096),3) + 22.6807411*pow((1 - (T)/647.096),3.5) - 15.9618719*pow((1 - (T)/647.096),4) + 1.80122502*pow((1 - (T)/647.096),7.5))*647.096/(T)))*InternationalSystemOfUnits_OnePascal)


/* 
 * References 
 * ----------
 * 
 * Buck, A. L. (1981), "New equations for computing vapor pressure and enhancement factor", J. Appl. Meteorol. 20: 1527-1532

 * Antoine, C. (1888), "Tensions des vapeurs; nouvelle relation entre les tensions et les températures" [Vapor Pressure: a new relationship between pressure and temperature], Comptes Rendus des Séances de l'Académie des Sciences (in French) 107: 681–684, 778–780, 836–837

 * W. Wagner and A. Pruss (1993) J. Phys. Chem. Reference Data, 22, 783-787.
 * 
 */

#endif
