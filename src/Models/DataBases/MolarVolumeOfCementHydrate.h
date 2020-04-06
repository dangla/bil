#ifndef MOLARVOLUMEOFCEMENTHYDRATE_H
#define MOLARVOLUMEOFCEMENTHYDRATE_H

#include "InternationalSystemOfUnits.h"

#define MolarVolumeOfCementHydrate_Unit \
        (InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter*InternationalSystemOfUnits_OneMeter/InternationalSystemOfUnits_OneMole)



#define MolarVolumeOfCementHydrate(I) \
        (MolarVolumeOfCementHydrate_##I*MolarVolumeOfCementHydrate_Unit)



/* Molar Volume of Cement Hydrates (m3/mole)
 * -----------------------------------------
 * 
 * Reference:
 * 
 * Lothenbach, Kulik, Matschei, Balonis, Baquerizo, Dilnesa, Miron, Myers
 * Cemdata18: A chemical thermodynamic database for hydrated Portland cements and 
 * alkali-activated materials, 
 * CCR 115 (2019) 472-506. 
 */
 
 
 


/*
 * Cement chemistry notation:
 * H  = H2O  ; K  = K2O  ; N  = Na2O  ;
 * C  = CaO  ; S  = SiO2 ; A  = Al2O3 ; 
 * c  = CO2  ; s  = SO3  ; F  = Fe2O3 ;
 */
 
 
 /* Convention
  * The symbol 'd' between two integers 'N' and 'M' like 'NdM' means 'N.M'.
  * The symbol 'h' after an integer 'N' like 'Nh' means 'N+0.5'.
  */
 

 
 
/** AFt-phases */

/* Al-Ettringite */
/* C3A.3(Cs).H32 = Ca6.Al2.3(SO4).12(OH)·26(H2O) */
#define MolarVolumeOfCementHydrate_AlEttringite       (707.e-6)

/* Fe-Ettringite */
/* C3F.3(Cs).H32 = Ca6.Fe2.3(SO4).12(OH)·26(H2O) */
#define MolarVolumeOfCementHydrate_FeEttringite       (717.e-6)

/* Thaumasite */
/* CS.Cs.Cc.H15 = Ca3.(SiO3).(SO4).(CO3)·15(H2O) */
#define MolarVolumeOfCementHydrate_Thaumasite         (330.e-6)



/** AFm-phases */

/* Friedel's salt */
/* C3AH10.CaCl2 = C4ACl2H10 = Ca4.Al2.Cl2.12(OH)·4(H2O) */
#define MolarVolumeOfCementHydrate_FriedelSalt        (272.e-6)

/* Kuzel's salt */
/* C3AH12.0.5(Cs).0.5(CaCl2) = Ca4.Al2.Cl(SO4)0.5.12(OH)·6(H2O) */
#define MolarVolumeOfCementHydrate_KuzelSalt          (289.e-6)

#define MolarVolumeOfCementHydrate_C4AH19             (369.e-6)

#define MolarVolumeOfCementHydrate_C4AH13             (274.e-6)

#define MolarVolumeOfCementHydrate_C4AH11             (257.e-6)

#define MolarVolumeOfCementHydrate_C2AH7h             (180.e-6)

#define MolarVolumeOfCementHydrate_CAH10              (193.e-6)

#define MolarVolumeOfCementHydrate_C4Ac0hH12          (285.e-6)

#define MolarVolumeOfCementHydrate_C4Ac0hH10.5        (261.e-6)

#define MolarVolumeOfCementHydrate_C4Ac0hH9           (249.e-6)

#define MolarVolumeOfCementHydrate_C4AcH11            (262.e-6)

#define MolarVolumeOfCementHydrate_C4AcH9             (234.e-6)

#define MolarVolumeOfCementHydrate_C4AsH16            (351.e-6)

#define MolarVolumeOfCementHydrate_C4AsH14            (332.e-6)

#define MolarVolumeOfCementHydrate_C4AsH12            (310.e-6)

#define MolarVolumeOfCementHydrate_C4AsH10h           (282.e-6)

#define MolarVolumeOfCementHydrate_C4AsH9             (275.e-6)

#define MolarVolumeOfCementHydrate_C2ASH8             (216.e-6)

#define MolarVolumeOfCementHydrate_C2ASH7             (215.e-6)

#define MolarVolumeOfCementHydrate_C2ASH5h            (213.e-6)

#define MolarVolumeOfCementHydrate_C4As0hClH12        (289.e-6)

#define MolarVolumeOfCementHydrate_C4ACl2H10          (272.e-6)

#define MolarVolumeOfCementHydrate_C4A(NO3)2H10       (296.e-6)

#define MolarVolumeOfCementHydrate_C4A(NO2)2H10       (275.e-6)

#define MolarVolumeOfCementHydrate_C4FH13             (286.e-6)

#define MolarVolumeOfCementHydrate_C4Fc0hH10          (273.e-6)

#define MolarVolumeOfCementHydrate_C4FcH12            (292.e-6)

#define MolarVolumeOfCementHydrate_C4FsH12            (321.e-6)

#define MolarVolumeOfCementHydrate_C2FSH8             (??)

#define MolarVolumeOfCementHydrate_C4FCl2H10          (278.e-6)



/** Hydrogarnets */

#define MolarVolumeOfCementHydrate_C3AH6              (150.e-6)

#define MolarVolumeOfCementHydrate_C3AS0d41H5d18      (146.e-6)

#define MolarVolumeOfCementHydrate_C3AS0d84H4d32      (142.e-6)

#define MolarVolumeOfCementHydrate_C3FH6              (155.e-6)




/** Sulfates */

/* Anhydrite */
#define MolarVolumeOfCementHydrate_Cs                 (46.e-6)

/* Gypsum */
#define MolarVolumeOfCementHydrate_CsH2               (75.e-6)
#define MolarVolumeOfCementHydrate_Gypsum             (75.e-6)

/* Hemihydrate */
#define MolarVolumeOfCementHydrate_CsH0h              (62.e-6)  

/* Syngenite K2Ca(SO4)2·H2O */
#define MolarVolumeOfCementHydrate_KsCsH2             (128.e-6)
#define MolarVolumeOfCementHydrate_Syngenite          (128.e-6)



/** Hydroxides */

/* Gibbsite (Al(OH)3)*/
#define MolarVolumeOfCementHydrate_Gibbsite           (32.e-6)
#define MolarVolumeOfCementHydrate_AH3                (64.e-6)

/* Amorphous Gibbsite */
#define MolarVolumeOfCementHydrate_AH3am              (64.e-6)

/* Microcrystalline Gibbsite */
#define MolarVolumeOfCementHydrate_AH3mic             (64.e-6)

/* Amorphous Fe(OH)3 */
#define MolarVolumeOfCementHydrate_FH3am              (2x)

/* Microcrystalline Fe(OH)3 */
#define MolarVolumeOfCementHydrate_FH3mic             (2x)

/* Goethite (FeOOH)  */
#define MolarVolumeOfCementHydrate_Geothite           (21.e-6)
#define MolarVolumeOfCementHydrate_FH                 (42.e-6)

/* Microcrystalline Geothite */
#define MolarVolumeOfCementHydrate_FHmic              (42.e-6)

/* Portlandite */
#define MolarVolumeOfCementHydrate_CH                 (33.e-6)
#define MolarVolumeOfCementHydrate_Portlandite        (33.e-6)

/* Amorphous Silica */
#define MolarVolumeOfCementHydrate_Sam                (29.e-6)

/* Quartz */
#define MolarVolumeOfCementHydrate_S                  (29.e-6)
#define MolarVolumeOfCementHydrate_Quartz             (29.e-6)



/** Clinkers */

#define MolarVolumeOfCementHydrate_C3S                (73.e-6)

#define MolarVolumeOfCementHydrate_C2S                (52.e-6)

#define MolarVolumeOfCementHydrate_C3A                (89.e-6)

#define MolarVolumeOfCementHydrate_C12A7              (518.e-6) 

#define MolarVolumeOfCementHydrate_CA                 (54.e-6)

#define MolarVolumeOfCementHydrate_CA2                (89.e-6)

#define MolarVolumeOfCementHydrate_C4AF               (130.e-6)

/* Lime (CaO) */
#define MolarVolumeOfCementHydrate_C                  (17.e-6)
#define MolarVolumeOfCementHydrate_Lime               (17.e-6)

/* Arcanite (K2SO4) */
#define MolarVolumeOfCementHydrate_Ks                 (66.e-6)
#define MolarVolumeOfCementHydrate_Arcanite           (66.e-6)

/* K2O */
#define MolarVolumeOfCementHydrate_K                  (40.e-6)

/* Thenardite (Na2SO4)  */
#define MolarVolumeOfCementHydrate_Ns                 (53.e-6)
#define MolarVolumeOfCementHydrate_Thenardite         (53.e-6)

/* Na2O */
#define MolarVolumeOfCementHydrate_N                  (25.e-6)





#endif
