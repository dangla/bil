#ifndef PARTIALMOLARVOLUMEOFMOLECULEINWATER_H
#define PARTIALMOLARVOLUMEOFMOLECULEINWATER_H

#include "InternationalSystemOfUnits.h"

#define PartialMolarVolumeOfMoleculeInWater_Unit \
       (InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter/InternationalSystemOfUnits_OneMole)



#define PartialMolarVolumeOfMoleculeInWater(I) \
      (PartialMolarVolumeOfMolecule_##I*PartialMolarVolumeOfMoleculeInWater_Unit)



/* Partial Molar Volume of Molecules in Water (m3/mole) at 25°C
 * ref [1] */
/* After Lothenbach and Kulik (2001) */

/** Compounds of type I */

/* Water */
#define PartialMolarVolumeOfMolecule_H           (-5.5e-6)
#define PartialMolarVolumeOfMolecule_OH          (-4.71e-6) 
#define PartialMolarVolumeOfMolecule_H2O         (18.e-6)

/* Calcium compounds */
#define PartialMolarVolumeOfMolecule_Ca          (-18.44e-6)
#define PartialMolarVolumeOfMolecule_CaOH        (5.76e-6)
#define PartialMolarVolumeOfMolecule_CaO2H2      (26.20e-6)

/* Silicon compounds */
#define PartialMolarVolumeOfMolecule_H4SiO4      (16.06e-6)
#define PartialMolarVolumeOfMolecule_H3SiO4      (4.53e-6) 
#define PartialMolarVolumeOfMolecule_H2SiO4      (34.13e-6)

/* Sodium compounds */
#define PartialMolarVolumeOfMolecule_Na          (-1.21e-6)
#define PartialMolarVolumeOfMolecule_NaOH        (3.511e-6)

/* Potassium compounds */
#define PartialMolarVolumeOfMolecule_K           (9.06e-6)
#define PartialMolarVolumeOfMolecule_KOH         (3.511e-6)

/* Carbon compounds */
#define PartialMolarVolumeOfMolecule_CO2         (32.81e-6)
#define PartialMolarVolumeOfMolecule_H2CO3       (50.e-6)
#define PartialMolarVolumeOfMolecule_HCO3        (24.21e-6) 
#define PartialMolarVolumeOfMolecule_CO3         (-6.06e-6)

/* Sulfur compounds */
#define PartialMolarVolumeOfMolecule_H2SO4       (50.e-6)
#define PartialMolarVolumeOfMolecule_HSO4        (24.21e-6) 
#define PartialMolarVolumeOfMolecule_SO4         (22.9e-6)

/* Chlorine compounds */
#define PartialMolarVolumeOfMolecule_Cl          (17.34e-6)

/* Aluminium compounds */
#define PartialMolarVolumeOfMolecule_Al          (20.e-6)    /* Guess */
#define PartialMolarVolumeOfMolecule_AlO4H4      (42.3e-6)   /* [2] */


/** Compounds of type II */

/* Calcium-Carbon compounds */
#define PartialMolarVolumeOfMolecule_CaCO3       (-15.65e-6)
#define PartialMolarVolumeOfMolecule_CaHCO3      (13.33e-6)

/* Calcium-Silicon compounds */
#define PartialMolarVolumeOfMolecule_CaH2SiO4    (15.69e-6)
#define PartialMolarVolumeOfMolecule_CaH3SiO4    (-6.74e-6)

/* Sodium-Carbon compounds */
#define PartialMolarVolumeOfMolecule_NaCO3       (-0.42e-6)
#define PartialMolarVolumeOfMolecule_NaHCO3      (32.32e-6)

/* Calcium-Sulfur compounds */
#define PartialMolarVolumeOfMolecule_CaSO4       (26.20e-6)
#define PartialMolarVolumeOfMolecule_CaHSO4      (27.e-6)


/*
 * [1] Y. Marcus. The Standard Partial Molar Volumes of Ions in Solution. Part 4. Ionic Volumes in Water at
0-100 °C, J. Phys. Chem. B 2009, 113, 10285-10291.
 * [2] B. Sanjuan, G. Michard. Determination of the partial molar volume of aluminate ion (Al(OH)4-) at infinite dilution in water at 25.degree.C, J. Chem. Eng. Data, 1988, 33 (2), pp 78-80, DOI: 10.1021/je00052a003.
 * 
 */


#endif
