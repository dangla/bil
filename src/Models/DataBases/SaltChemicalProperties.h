#ifndef SALTCHEMICALPROPERTIES_H
#define SALTCHEMICALPROPERTIES_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct SaltChemicalProperties_s     ; 
typedef struct SaltChemicalProperties_s     SaltChemicalProperties_t ;


extern SaltChemicalProperties_t* SaltChemicalProperties(char*,double) ;


#define SaltChemicalProperties_GetValenceNumberOfCation(scp)  ((scp)->z_cation)
#define SaltChemicalProperties_GetValenceNumberOfAnion(scp)   ((scp)->z_anion)
#define SaltChemicalProperties_GetStoichiometricCoefficientOfCation(scp)  ((scp)->nu_cation)
#define SaltChemicalProperties_GetStoichiometricCoefficientOfAnion(scp)   ((scp)->nu_anion)
#define SaltChemicalProperties_GetStoichiometricCoefficientOfWater(scp)   ((scp)->nu_h2o)
#define SaltChemicalProperties_GetMolarMass(scp)                    ((scp)->molarmass)
#define SaltChemicalProperties_GetMolarVolume(scp)                  ((scp)->molarvolume)
#define SaltChemicalProperties_GetSolubilityProduct(scp)            ((scp)->solubilityproduct)
#define SaltChemicalProperties_GetPartialMolarVolumeOfCation(scp)   ((scp)->v_cation)
#define SaltChemicalProperties_GetPartialMolarVolumeOfAnion(scp)    ((scp)->v_anion)


struct SaltChemicalProperties_s {
  double z_cation ;
  double z_anion ;
  double nu_cation ;
  double nu_anion ;
  double nu_h2o ;
  double molarmass;
  double molarvolume ;
  double solubilityproduct ;
  double v_cation ;
  double v_anion ; 
} ;


#endif
