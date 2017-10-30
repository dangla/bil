#include "Algos.h"


#define CementSolutionChemistry_GetIndexOf(CPD) \
        Utils_CAT(CementSolutionChemistry_B_,CPD)


#define CementSolutionChemistry_ENUM \
        Tuple_SEQ(Algos_MAP(CementSolutionChemistry_ListOfCompounds,CementSolutionChemistry_GetIndexOf))

//#define CementSolutionChemistry_ListOfCompounds CementSolutionChemistry_ListOfCompounds_

#define CementSolutionChemistry_ListOfCompounds (\
        H2O,H,OH \
       ,Ca,CaOH,CaO2H2 \
       ,H2SiO4,H3SiO4,H4SiO4 \
       ,CaH2SiO4,CaH3SiO4 \
       ,Na,NaOH \
       ,K,KOH \
       ,H2CO3,HCO3,CO3,CO2 \
       ,CaHCO3,CaCO3 \
       ,NaHCO3,NaCO3 \
       ,H2SO4,HSO4,SO4 \
       ,SO3,S2O3 \
       ,H2S,HS,S,S0 \
       ,CaHSO4,CaSO4 \
       ,Cl \
       ,Al,AlO4H4 \
       )




enum CementSolutionChemistry_e {
  CementSolutionChemistry_ENUM
} ;
