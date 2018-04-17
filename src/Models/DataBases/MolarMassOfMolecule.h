#ifndef MOLARMASSOFMOLECULE_H
#define MOLARMASSOFMOLECULE_H

#include "InternationalSystemOfUnits.h"

#define MolarMassOfMolecule_Unit \
       (InternationalSystemOfUnits_OneKilogram/InternationalSystemOfUnits_OneMole)



#define MolarMassOfMolecule(I) \
       (MolarMassOfMolecule_##I*MolarMassOfMolecule_Unit)



/* Molar Mass of Chemicals (kg/mol) */
/* Water */
#define MolarMassOfMolecule_H           (1.e-3)
#define MolarMassOfMolecule_O           (16.e-3)
#define MolarMassOfMolecule_OH          (17.e-3)
#define MolarMassOfMolecule_H2O         (18.e-3)

/* Calcium compounds */
#define MolarMassOfMolecule_Ca          (40.e-3)
#define MolarMassOfMolecule_CaO2H2      (74.e-3)
#define MolarMassOfMolecule_CaO         (56.e-3)
#define MolarMassOfMolecule_CaOH        (57.e-3)

/* Silicium compounds */
#define MolarMassOfMolecule_Si          (28.e-3)
#define MolarMassOfMolecule_SiO2        (60.e-3)
#define MolarMassOfMolecule_H2SiO4      (94.e-3)
#define MolarMassOfMolecule_H3SiO4      (95.e-3)
#define MolarMassOfMolecule_H4SiO4      (96.e-3)

/* Sodium compounds */
#define MolarMassOfMolecule_Na          (23.e-3)
#define MolarMassOfMolecule_NaOH        (40.e-3)

/* Potassium compounds */
#define MolarMassOfMolecule_K           (39.e-3)
#define MolarMassOfMolecule_KOH         (56.e-3)

/* Carbon compounds */
#define MolarMassOfMolecule_C           (12.e-3)
#define MolarMassOfMolecule_CO2         (44.e-3)
#define MolarMassOfMolecule_H2CO3       (62.e-3)
#define MolarMassOfMolecule_HCO3        (61.e-3)
#define MolarMassOfMolecule_CO3         (60.e-3)

/* Aluminium compounds */
#define MolarMassOfMolecule_Al          (26.98e-3)
#define MolarMassOfMolecule_AlO4H4      (95.01e-3)
#define MolarMassOfMolecule_Al2O3       (101.96e-3)

/* Chlorine compounds */
#define MolarMassOfMolecule_Cl          (35.45e-3)

/* Sulfur compounds */
#define MolarMassOfMolecule_S           (32.e-3)
#define MolarMassOfMolecule_HS          (33.e-3)
#define MolarMassOfMolecule_H2S         (34.e-3)
#define MolarMassOfMolecule_H2SO4       (98.e-3)
#define MolarMassOfMolecule_HSO4        (97.e-3)
#define MolarMassOfMolecule_SO4         (96.e-3)
#define MolarMassOfMolecule_SO3         (80.e-3)

/* Iron compounds */
#define MolarMassOfMolecule_Fe          (55.85e-3)
#define MolarMassOfMolecule_Fe2O3       (159.7e-3)

/* Calcium-Silicium compounds */
#define MolarMassOfMolecule_CaH2SiO4    (134.e-3)
#define MolarMassOfMolecule_CaH3SiO4    (135.e-3)

/* Calcium-Carbon compounds */
#define MolarMassOfMolecule_CaCO3       (100.e-3)
#define MolarMassOfMolecule_CaHCO3      (101.e-3)

/* Sodium-Carbon compounds */
#define MolarMassOfMolecule_NaHCO3      (84.e-3)
#define MolarMassOfMolecule_NaCO3       (83.e-3)

/* Calcium-Sulfur compounds */
#define MolarMassOfMolecule_CaSO4       (136.e-3)
#define MolarMassOfMolecule_CaHSO4      (137.e-3)
#endif
