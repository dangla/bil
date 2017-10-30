/* General features of the model:
 * To be completed!!!!!
 * This model is from M70con
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "FVM.h"

#define TITLE   "Sulfuric acid attack of concrete (Jan. 2015)" 
#define AUTHORS "Yuan-Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ     (6)
#define NVE     (31)
#define NVI     (28)
#define NV0     (2)

/* Equation Indexes */
#define E_S     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)
#define E_K     (4)
#define E_Cl    (5)

/* Primary Unknown Indexes */
#define U_C_H2SO4  (0)
#define U_PSI      (1)
#define U_ZN_Ca_S  (2)
#define U_ZN_Si_S  (3)
#define U_C_K      (4)
#define U_C_Cl     (5)


#define NOLOG_U     1
#define LOG_U       2
#define Ln10        (2.302585093)
#define U_H2SO4     LOG_U
#define EXPLICIT  1
#define IMPLICIT  2
#define U_PHI     IMPLICIT


/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWN_n(n,i)   (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])



#if (U_H2SO4 == LOG_U)
  #define LogC_H2SO4(n)   (UNKNOWN(n,U_C_H2SO4))
  #define LogC_H2SO4n(n)  (UNKNOWN_n(n,U_C_H2SO4))
  #define C_H2SO4(n)      (exp(Ln10*UNKNOWN(n,U_C_H2SO4)))
  #define C_H2SO4n(n)     (exp(Ln10*UNKNOWN_n(n,U_C_H2SO4)))
#else
  #define C_H2SO4(n)      (UNKNOWN(n,U_C_H2SO4))
  #define C_H2SO4n(n)     (UNKNOWN_n(n,U_C_H2SO4))
  #define LogC_H2SO4(n)   (log10(UNKNOWN(n,U_C_H2SO4)))
  #define LogC_H2SO4n(n)  (log10(UNKNOWN_n(n,U_C_H2SO4)))
#endif
#define ZN_Ca_S(n)   (UNKNOWN(n,U_ZN_Ca_S))
#define ZN_Si_S(n)   (UNKNOWN(n,U_ZN_Si_S))
#define PSI(n)       (UNKNOWN(n,U_PSI))
#define C_K(n)       (UNKNOWN(n,U_C_K))
#define C_Cl(n)      (UNKNOWN(n,U_C_Cl))

#define ZN_Ca_Sn(n)  (UNKNOWN_n(n,U_ZN_Ca_S))
#define ZN_Si_Sn(n)  (UNKNOWN_n(n,U_ZN_Si_S))
#define PSIn(n)      (UNKNOWN_n(n,U_PSI))
#define C_Kn(n)      (UNKNOWN_n(n,U_C_K))
#define C_Cln(n)     (UNKNOWN_n(n,U_C_Cl))

#define N_S(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define N_Si(n)    (f[(6+n)])
#define N_K(n)     (f[(8+n)])
#define N_Cl(n)    (f[(10+n)])
#define W_S        (f[12])
#define W_q        (f[13])
#define W_Ca       (f[14])
#define W_Si       (f[15])
#define W_K        (f[16])
#define W_Cl       (f[17])
#define N_CH(n)    (f[(18+n)])
#define N_CSH2(n)  (f[(20+n)])
#define PHI(n)     (f[(22+n)])
#define PHI_C(n)   (f[(24+n)])
#define D_CON(n)   (f[(26+n)])

#define N_Sn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_Cln(n)   (f_n[(10+n)])
#define N_CHn(n)   (f_n[(18+n)])
#define N_CSH2n(n) (f_n[(20+n)])
#define PHIn(n)    (f_n[(22+n)])
#define PHI_Cn(n)  (f_n[(24+n)])
#define D_CONn(n)  (f_n[(26+n)])

#define KF_OH       (va[(0)])
#define KF_H        (va[(1)])
#define KF_CO2      (va[(2)])
#define KF_H2SO4    (va[(3)])
#define KF_HSO4     (va[(4)])
#define KF_SO4      (va[(5)])
#define KF_Ca       (va[(6)])
#define KF_CaHSO4   (va[(7)])
#define KF_CaH3SiO4 (va[(8)])
#define KF_H3SiO4   (va[(9)])
#define KF_H4SiO4   (va[(10)])
#define KF_H2SiO4   (va[(11)])
#define KF_CaH2SiO4 (va[(12)])
#define KF_CaSO4aq  (va[(13)])
#define KF_CaOH     (va[(14)])
#define KF_K        (va[(15)])
#define KF_Cl       (va[(16)])

#define Kpsi_OH       (va[(17)])
#define Kpsi_H        (va[(18)])
#define Kpsi_HSO4     (va[(19)])
#define Kpsi_SO4      (va[(20)])
#define Kpsi_Ca       (va[(21)])
#define Kpsi_CaHSO4   (va[(22)])
#define Kpsi_CaH3SiO4 (va[(23)])
#define Kpsi_H3SiO4   (va[(24)])
#define Kpsi_q        (va[(25)])
#define Kpsi_H2SiO4   (va[(26)])
#define Kpsi_CaOH     (va[(27)])
#define Kpsi_K        (va[(28)])
#define Kpsi_Cl       (va[(29)])
#define KD_CSH2       (va[(30)])

#define V_Cem0(n)     (v0[(0+n)])


/*
  Aqueous solution
*/

/* Charge of ions */
#include "ElectricChargeOfIonInWater.h"
#define z_ca          ElectricChargeOfIonInWater(Ca)
#define z_h           ElectricChargeOfIonInWater(H)
#define z_oh          ElectricChargeOfIonInWater(OH)
#define z_h3sio4      ElectricChargeOfIonInWater(H3SiO4)
#define z_cah3sio4    ElectricChargeOfIonInWater(CaH3SiO4)
#define z_h2sio4      ElectricChargeOfIonInWater(H2SiO4)
#define z_caoh        ElectricChargeOfIonInWater(CaOH)
#define z_k           ElectricChargeOfIonInWater(K)
#define z_cl          ElectricChargeOfIonInWater(Cl)
#define z_hso4        ElectricChargeOfIonInWater(HSO4)
#define z_so4         ElectricChargeOfIonInWater(SO4)
#define z_cahso4      ElectricChargeOfIonInWater(CaHSO4)


/* volumes molaires partiels des ions (dm3/mole) from [Millero F J,Partial molar volum of ions in seawater]*/
#ifdef NOTDEFINED
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_h2o      (18.e-3)
#define v_co2      (32.81e-3)
#define v_h2so4    (50.e-3)
#define v_hso4     (24.21e-3)   /* d'apres Lothenbach */
#define v_so4      (22.9e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)
#define v_caso4aq  (26.20e-3)    /* a modifier */
#define v_caoh     (26.20e-3)    /* a modifier */
#define v_k        (43.93e-3)    /* d'apres Antoine */
#define v_cl       (43.93e-3)    /*a modifier*/
#define v_koh      (27.44e-3)    /* d'apres Antoine */
#endif


/* Equilibrium constant of homogeneous aqueous reaction */
#define K_h2o      (1.e-14)          /* H2O = OH[-] + H[+] */

#define K_hso4     (1.e6)            /* H2SO4   = HSO4[-] + H[+] */
#define K_so4      (1.e-2)           /* HSO4[-] = SO4[2-] + H[+] */

#define K_h2sio4   (4.68)            /* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define K_h4sio4   (6.45e9)          /* H3SiO4[-] + H[+]  = H4SiO4 */
#define K_h3sio4   (1.55e-10)        /* H4SiO4    = H3SiO4[-] + H[+] */

#define K_cahso4   (1.276e+1)        /* Ca[2+] + HSO4[-]    = CaHSO4[+] */
#define K_caso4aq  (1.4e+3)          /* Ca[2+] + SO4[2-]    = CaSO4[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-]      = CaOH[+] */


/*
  Solids
  CH   = Calcium Hydroxide (Portlandite)
  CSH2 = Calcium Sulfate Dihydrate (Gypsum)
  CSH  = Calcium Silicates Hydrate
  SH   = Amorphous Silica Gel
*/

/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(q) (Curve_ComputeValue(Element_GetCurve(el),q))
#define WaterSiliconRatioInCSH(q)   (Curve_ComputeValue(Element_GetCurve(el) + 1,q))
#define MolarVolumeOfCSH(q)    (Curve_ComputeValue(Element_GetCurve(el) + 2,q))



/* S-H Properties */
/* Equilibrium constant */
#define K_SH       (1.93642e-3)      /* SHt = H4SiO4 + (t-2)H2O */
/* Saturation Degree of Dissolved S-H */
#define S_SHeq(q)     (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define SaturationDegreeOfSH(s_ch,zn_si_s)       (NEGEXP(zn_si_s)*S_SHeq(s_ch))
/* Ion Activity Product of Dissolved S-H */
#define IonActivityProductOfSH(s_sh)             (K_SH*(s_sh))


/* CH Properties */
/* Equilibrium constant */
#define K_CH       (6.456e-6)        /* CH  = Ca[2+] + 2OH[-] */
/* Molar volume of CH solid (dm3/mole) */
#define V_CH       (33.e-3)      /* (33.e-3) */
/* Saturation Degree of Dissolved CH */
#define SaturationDegreeOfCH(zc_h2so4,zn_ca_s)  (NEGEXP(zn_ca_s)/MAX(zc_h2so4,1.))
/* Ion Activity Product of Dissolved CH */
#define IonActivityProductOfCH(s_ch)            (K_CH*(s_ch))


/* CSH2 Properties */
/* Equilibrium constant */
#define K_CSH2     (2.5e-5)         /* CSH2 = Ca[2+] + SO4[2-] + 2H2O */
/* Threshold of H2SO4 concentration */
#define C_H2SO4_eq    (K_h2o*K_h2o*K_CSH2/(K_hso4*K_so4*K_CH))
/* Molar volume of CSH2 crystal (dm3/mole) */
#define V_CSH2     (75.e-3)      /* (75.e-3) */
/* Saturation Degree of Dissolved CSH2 */
//#define SaturationDegreeOfCSH2(zc_h2so4,zn_ca_s)   (NEGEXP(zn_ca_s)*MIN(zc_h2so4,1.))
#define SaturationDegreeOfCSH2(zc_h2so4,s_ch)      ((s_ch)*(zc_h2so4))
/* Ion Activity Product of Dissolved CSH2 */
#define IonActivityProductOfCSH2(s_csh2)           (K_CSH2*(s_csh2))


/* Element contents in solid phases  */
#define CalciumContentInCHAndCSH2(zn_ca_s) (n_ca_ref*MAX(zn_ca_s,0.))
#define SiliconContentInCSH(zn_si_s)       (n_si_ref*MAX(zn_si_s,0.))


/* Gypsum-based porous material properties */
/* Porosity of gypsum-based porous material (-) */
#define PHI_Gyp    (0.85)
/* Molar volume of gypsum-based porous material */
#define V_Gyp      (V_CSH2/(1 - PHI_Gyp))
/* Compressibility of gypsum-based porous material (MPa) */
#define CC_Gyp     (5.e+06)      /* (??) */


/* Concrete properties */
#define CC_Con     (5.e+9)       /* (??) */
#define CS_Con     (3.5e+6)      /* Attention the pressure should be lower than CC_Gyp * PHI_Gyp = 4.25e6 here */


/* Concentration of OH computed from electroneutrality */
#define ConcentrationOfOHInLiquid(A,B,C,D,E)  concentration_oh(A,B,C,D,E)


/* Function used above */
#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])



//#include "PCM.h"
/* Fonctions */
static int    pm(const char *s) ;


static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
//static PCM_ComputeComponents_t ComputeComponents ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
//static PCM_ComputeSecondaryComponents_t ComputeSecondaryComponents ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;
/* #define ComputeComponentDerivatives(a,b,c,d,e) (PCM_ComputeComponentDerivatives(a,ComputeSecondaryComponents,b,c,d,e)) */


static void    ComputeTransferCoefficients(Element_t*) ;
static double* ComputeComponentFluxes(Element_t*,int,int) ;
//#define ComputeComponentFluxes(a,i,j) (PCM_ComputeComponentFluxes(a,ComputeFluxes,i,j))

static double* ComputeFluxes(Element_t*,double*,int,int) ;
//static PCM_ComputeFluxes_t ComputeFluxes ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;

static double  concentration_oh(double,double,double,double,double) ;
static double  poly4(double,double,double,double,double) ;



/* Parameters */
static double phi0,c_h2so4_eq,t_ch,t_csh2 ;
static double n_ca_ref,n_si_ref ;

static double d_h,d_oh ;
static double d_ca,d_caoh ;
static double d_h2so4,d_hso4,d_so4 ;
static double d_cahso4 ;
static double d_h4sio4,d_h3sio4,d_h2sio4 ;
static double d_cah2sio4,d_cah3sio4 ;
static double d_caso4aq ;
static double d_k,d_koh ;
static double d_cl ;

static double FsRT ;






#define dm2      (1.e2)   /* This is 1 m2 (SI surface unit) in dm2 */
#include "DiffusionCoefficientOfMoleculeInWater.h"
#include "PhysicalConstant.h"

void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (dm2/s) */
  d_oh         = DiffusionCoefficientOfMoleculeInWater(OH,TK)*dm2 ;
  d_h          = DiffusionCoefficientOfMoleculeInWater(H,TK)*dm2 ;
  
  d_ca         = DiffusionCoefficientOfMoleculeInWater(Ca,TK)*dm2 ;
  d_caoh       = DiffusionCoefficientOfMoleculeInWater(CaOH,TK)*dm2 ;
  
  d_h4sio4     = DiffusionCoefficientOfMoleculeInWater(H4SiO4,TK)*dm2 ;
  d_h3sio4     = DiffusionCoefficientOfMoleculeInWater(H3SiO4,TK)*dm2 ;
  d_h2sio4     = DiffusionCoefficientOfMoleculeInWater(H2SiO4,TK)*dm2 ;
  
  d_h2so4      = (1.5e-7) ;
  d_hso4       = (1.385e-7) ;
  d_so4        = (5.32e-8) ;
  
  d_cah2sio4   = DiffusionCoefficientOfMoleculeInWater(CaH2SiO4,TK)*dm2 ;
  d_cah3sio4   = DiffusionCoefficientOfMoleculeInWater(CaH3SiO4,TK)*dm2 ;

  d_cahso4     = (1.07e-7) ;
  d_caso4aq    = (1.43e-7) ;
  
  d_k          = DiffusionCoefficientOfMoleculeInWater(K,TK)*dm2 ;
  d_koh        = DiffusionCoefficientOfMoleculeInWater(KOH,TK)*dm2 ;
  
  d_cl         = DiffusionCoefficientOfMoleculeInWater(Cl,TK)*dm2 ;
  
  /* Physical constants */
  {
    double RT      = PhysicalConstant(PerfectGasConstant)*TK*1.e3 ;
    double Faraday = PhysicalConstant(Faraday)*1.e3 ;
    
    FsRT     = Faraday/RT ;
  }
}


#define NN                (2)
#define NbOfComponents    (46)
static double Components[NN][NbOfComponents] ;
static double dComponents[NbOfComponents] ;

#define I_C_OH         (6)
#define I_C_H          (7)

#define I_C_H2SO4      (8)
#define I_C_HSO4       (9)
#define I_C_SO4        (10)

#define I_C_Ca         (11)
#define I_C_CaOH       (12)

#define I_C_CaHSO4     (13)
#define I_C_CaSO4aq    (14)

#define I_C_H2SiO4     (15)
#define I_C_H3SiO4     (16)
#define I_C_H4SiO4     (17)

#define I_C_CaH2SiO4   (18)
#define I_C_CaH3SiO4   (19)

#define I_C_K          (20)

#define I_C_Cl         (21)

#define I_ZN_Ca_S      (22)
#define I_ZN_Si_S      (23)

#define I_S_CH         (24)

#define I_N_Q          (25)
#define I_N_S          (26)
#define I_N_Ca         (27)
#define I_N_Si         (28)
#define I_N_K          (29)
#define I_N_Cl         (30)

#define I_PSI          (31)

#define I_N_CH         (32)
#define I_N_CSH2       (33)
#define I_D_CON        (34)

#define I_N_CHn        (35)
#define I_N_CSH2n      (36)
#define I_D_CONn       (37)

#define I_PHI          (38)
#define I_PHI_C        (39)

#define I_PHIn         (40)
#define I_PHI_Cn       (41)

#define I_V_Cem        (42)

#define I_V_Cem0       (43)

#define I_N_Si_S       (44)

#define I_P_CSH2       (45)
  
  

#define NbOfComponentFluxes    (6)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_S           (0)
#define I_W_Ca          (1)
#define I_W_Si          (2)
#define I_W_K           (3)
#define I_W_Cl          (4)
#define I_W_q           (5)


int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0)        return (0) ;
  else if(strcmp(s,"N_CH") == 0)       return (1) ;
  else if(strcmp(s,"N_Si") == 0)       return (2) ;
  else if(strcmp(s,"C_H2SO4_eq") == 0) return (3) ;
  else if(strcmp(s,"T_CH") == 0)       return (4) ;
  else if(strcmp(s,"T_CSH2") == 0)     return (5) ;
  else if(strcmp(s,"surface") == 0)    return (6) ;
  else return(-1) ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  Model_GetNbOfVariables(model) = NbOfComponents ;
  Model_GetNbOfVariableFluxes(model) = NbOfComponentFluxes ;
  
  Model_CopyNameOfEquation(model,E_S, "sulfur") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_Si,"silicium") ;
  Model_CopyNameOfEquation(model,E_q, "charge") ;
  Model_CopyNameOfEquation(model,E_K, "potassium") ;
  Model_CopyNameOfEquation(model,E_Cl,"chlorine") ;

#if (U_H2SO4 == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_H2SO4,"logc_h2so4") ;
#else
  Model_CopyNameOfUnknown(model,U_C_H2SO4,"c_h2so4") ;
#endif
  Model_CopyNameOfUnknown(model,U_ZN_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,U_PSI,    "psi") ;
  Model_CopyNameOfUnknown(model,U_ZN_Si_S,"z_si") ;
  Model_CopyNameOfUnknown(model,U_C_K,    "c_k") ;
  Model_CopyNameOfUnknown(model,U_C_Cl,   "c_cl") ;

  {
    double temperature = 293 ;
    
    ComputePhysicoChemicalProperties(temperature) ;
  }
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 7 ;

  {
    /* Self-initialization */
    Material_GetProperty(mat)[pm("N_CH")]   = 1 ;
    Material_GetProperty(mat)[pm("N_Si")]   = 1 ;
    Material_GetProperty(mat)[pm("T_CH")]   = 600 ;
    Material_GetProperty(mat)[pm("T_CSH2")] = 600 ;

    Material_ScanProperties(mat,datafile,pm) ;

    t_ch      = Material_GetProperty(mat)[pm("T_CH")] ;
    t_csh2    = Material_GetProperty(mat)[pm("T_CSH2")] ;

    if(t_csh2  == 0.) Material_GetProperty(mat)[pm("T_CSH2")]  = t_ch ;

    Material_GetProperty(mat)[pm("C_H2SO4_eq")] = C_H2SO4_eq ;
  }
  
  return(NbOfProp) ;
}



int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The 6 equations are:\n") ;
  printf("\t- Mass balance of S      (sulfur)\n") ;
  printf("\t- Charge balance         (charge)\n") ;
  printf("\t- Mass balance of Ca     (calcium)\n") ;
  printf("\t- Mass balance of Si     (silicium)\n") ;
  printf("\t- Mass balance of K      (potassium)\n") ;
  printf("\t- Mass balance of Cl     (chlorine)\n") ;
  
  printf("\n") ;
  printf("The 6 primary unknowns are:\n") ;
  printf("\t- Sulfuric acid concentration     (c_h2so4 or logc_h2so4)\n") ;
  printf("\t- Electric potential              (psi)\n") ;
  printf("\t- Zeta unknown for calcium        (z_ca)\n") ;
  printf("\t- Zeta unknown for silicium       (z_si)\n") ;
  printf("\t- Potassium concentration         (c_k)\n") ;
  printf("\t- Chlorine concentration          (c_cl)\n") ;

  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length : dm !\n") ;
  printf("\t time   : s !\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1       # Contenu en CH (moles/L)\n") ;
  fprintf(ficd,"N_Si  = 2.4       # contenu en Si solide (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4       # contenu en K (moles/L)\n") ;
  fprintf(ficd,"N_Cl   = 0.4      # contenu en Cl (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5      # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"T_CSH2  = 1.e5    # Cinetique de dissolution de CSH2 (s)\n") ;
  fprintf(ficd,"Curves = solid    # File name: S_CH  C/S  H/S  V  S_SHt\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  return(0) ;
}


int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double *r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/* Initialise les variables du systeme (f,va) */ 
{
  double *f = Element_GetImplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  //PCM_t* pcm = PCM_GetInstance(el,u,u,f) ;
  int i ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2so4_eq  = GetProperty("C_H2SO4_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_csh2    = GetProperty("T_CSH2") ;
  
  
  /* Pre-initialization */
  for(i = 0 ; i < nn ; i++) {
    double c_h2so4    = (C_H2SO4(i) > 0.) ? C_H2SO4(i) : c_h2so4_eq ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zn_si_s    = ZN_Si_S(i) ;
    
    /* Liquid components */
    double zc_h2so4   = c_h2so4/c_h2so4_eq ;
    double s_ch       = SaturationDegreeOfCH(zc_h2so4,zn_ca_s) ;
    
    /* Solid contents */
    /* ... as components: CH, CSH2, CSH */
    double n_ch_csh2  = CalciumContentInCHAndCSH2(zn_ca_s) ;
    double n_ch       = (zc_h2so4 <= 1) ? n_ch_csh2  : 0 ;
    double n_csh2     = (zc_h2so4 >  1) ? n_ch_csh2  : 0 ;
    /* ... as elements: S, Ca, Si */
    double n_si_s     = SiliconContentInCSH(zn_si_s) ;
    /* ... as volume */
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double v_cem      = V_CH*n_ch + v_csh*n_si_s ;

    /* Porosity */
    double phi_c = phi0 ;
    double phi   = phi_c - V_CSH2*n_csh2 ;


    /* Back up what is needed to compute components */
    N_CH(i)    = n_ch ;
    N_CSH2(i)  = n_csh2 ;
    PHI(i)     = phi ;
    PHI_C(i)   = phi_c ;
    D_CON(i)   = 0 ;

    V_Cem0(i)  = v_cem ;


#if (U_H2SO4 == LOG_U)
    LogC_H2SO4(i) = log10(c_h2so4) ;
#else
    C_H2SO4(i)   = c_h2so4 ;
#endif
  }
  
  
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    //double *x = ComputeComponents(pcm,0,i) ;
    double *x = ComputeComponents(el,u,f,0,i) ;

    
    /* Back up */
    N_S(i)  = x[I_N_S] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ;
    N_Cl(i) = x[I_N_Cl] ;
    /* charge density */
    N_q(i)  = x[I_N_Q] ;

    
    N_CH(i)    = x[I_N_CH] ;
    N_CSH2(i)  = x[I_N_CSH2] ;
    PHI(i)     = x[I_PHI] ;
    PHI_C(i)   = x[I_PHI_C] ;
    D_CON(i)   = x[I_D_CON] ;
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  ComputeTransferCoefficients(el) ;

  /* Flux */
  {
    double* w = ComputeComponentFluxes(el,0,1) ;

    W_S     = w[I_W_S  ] ;
    W_Ca    = w[I_W_Ca ] ;
    W_Si    = w[I_W_Si ] ;
    W_q     = w[I_W_q  ] ;
    W_K     = w[I_W_K  ] ;
    W_Cl    = w[I_W_Cl ] ;
  }
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Thermes explicites (va)  */
{
  double* f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  //PCM_t* pcm = PCM_GetInstance(el,u,u,f) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2so4_eq  = GetProperty("C_H2SO4_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_csh2    = GetProperty("T_CSH2") ;
  
  /*
    Coefficients de transfert
  */
  {
    int i ;
    
    for(i = 0 ; i < 2 ; i++) {
      //ComputeComponents(pcm,0,i) ;
      ComputeComponents(el,u,f,0,i) ;
    }
  }
  ComputeTransferCoefficients(el) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  //PCM_t* pcm = PCM_GetInstance(el,u,u_n,f_n) ;
  int i ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2so4_eq  = GetProperty("C_H2SO4_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_csh2    = GetProperty("T_CSH2") ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
    //double *x = ComputeComponents(pcm,dt,i) ;
    double *x = ComputeComponents(el,u,f_n,dt,i) ;

    /* Back up */

    /* Molar contents */
    N_S(i)  = x[I_N_S] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ;
    N_Cl(i) = x[I_N_Cl] ;
    
    /* Charge density */
    N_q(i)  = x[I_N_Q] ;

    
    N_CH(i)    = x[I_N_CH] ;
    N_CSH2(i)  = x[I_N_CSH2] ;
    PHI(i)     = x[I_PHI] ;
    PHI_C(i)   = x[I_PHI_C] ;
    D_CON(i)   = x[I_D_CON] ;


    {
      double c_h2so4    = x[I_C_H2SO4] ;
      double zn_si_s    = x[I_ZN_Si_S] ;
      double n_si_s     = SiliconContentInCSH(zn_si_s) ;
      double n_ch       = x[I_N_CH] ;
      double n_csh2     = x[I_N_CSH2] ;
      double s_ch       = x[I_S_CH] ;
      double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
      double n_ca_s     = n_ch + n_csh2 + x_csh*n_si_s ;
      double phi        = x[I_PHI] ;
      
      if(c_h2so4 < 0. || n_ca_s < 0. || n_si_s < 0.) {
        double xx = Element_GetNodeCoordinate(el,i)[0] ;
        double c_h3sio4   = x[I_C_H3SiO4] ;
        double c_oh       = x[I_C_OH] ;
        double zn_ca_s    = x[I_ZN_Ca_S] ;
        printf("x         = %e\n",xx) ;
        printf("c_h2so4   = %e\n",c_h2so4) ;
        printf("n_csh2    = %e\n",n_csh2) ;
        printf("n_ca_s    = %e\n",n_ca_s) ;
        printf("n_si_s    = %e\n",n_si_s) ;
        printf("zn_si_s   = %e\n",zn_si_s) ;
        printf("zn_ca_s   = %e\n",zn_ca_s) ;
        printf("c_h3sio4  = %e\n",c_h3sio4) ;
        printf("c_oh      = %e\n",c_oh) ;
        return(-1) ;
      }
 
      if(phi < 0.) {
        double xx = Element_GetNodeCoordinate(el,i)[0] ;
        double phi_c      = x[I_PHI_C] ;
        double p_csh2     = x[I_P_CSH2] ;
        printf("phi       = %e\n",phi) ;
        printf("phi_c     = %e\n",phi_c) ;
        printf("CH        = %e\n",n_ch) ;
        printf("CSH2      = %e\n",n_csh2) ;
        printf("Si        = %e\n",n_si_s) ;
        printf("c_h2so4   = %e\n",c_h2so4) ;
        printf("x         = %e\n",xx) ;
        printf("p_csh2    = %e\n",p_csh2) ;
        return(-1) ;
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Flux */
  {
    double* w = ComputeComponentFluxes(el,0,1) ;

    W_S     = w[I_W_S  ] ;
    W_Ca    = w[I_W_Ca ] ;
    W_Si    = w[I_W_Si ] ;
    W_q     = w[I_W_q  ] ;
    W_K     = w[I_W_K  ] ;
    W_Cl    = w[I_W_Cl ] ;
  }

  return(0) ;
}



int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double c[4*NEQ*NEQ] ;
  int    i ;
  
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2so4_eq  = GetProperty("C_H2SO4_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_csh2    = GetProperty("T_CSH2") ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }
  


#if (U_H2SO4 == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_H2SO4)     *= Ln10*C_H2SO4(0) ;
      K(i,U_C_H2SO4+NEQ) *= Ln10*C_H2SO4(1) ;
    }
  }
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Boundary Surface Area */
  {
    double *area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }
  
  /*
    Mass balance of elements S, Ca, Si, K, Cl
  */
  R(0,E_S)  -= volume[0]*(N_S(0)  - N_Sn(0))  + dt*surf*W_S ;
  R(1,E_S)  -= volume[1]*(N_S(1)  - N_Sn(1))  - dt*surf*W_S ;
  
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  
  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;
  
  R(0,E_K)  -= volume[0]*(N_K(0)  - N_Kn(0))  + dt*surf*W_K ;
  R(1,E_K)  -= volume[1]*(N_K(1)  - N_Kn(1))  - dt*surf*W_K ; 
  
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*W_Cl ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*W_Cl ; 
  /*
    Conservation charge
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  //PCM_t* pcm = PCM_GetInstance(el,u,u,f) ;
  int    nso = 34 ;
  int    i ;

  /* if(Element_IsSubmanifold(el)) return(0) ; */
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2so4_eq  = GetProperty("C_H2SO4_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_csh2    = GetProperty("T_CSH2") ;


  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }


  /* output quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    //double *x = ComputeComponents(pcm,0,j) ;
    double *x = ComputeComponents(el,u,f,0,j) ;
    
    double c_h2so4    = x[I_C_H2SO4] ;
    double zn_si_s    = x[I_ZN_Si_S] ;
    double zn_ca_s    = x[I_ZN_Ca_S] ;
    double c_k        = x[I_C_K] ;
    double c_cl       = x[I_C_Cl] ;
    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hso4     = x[I_C_HSO4] ;
    double c_so4      = x[I_C_SO4] ;
    double c_ca       = x[I_C_Ca] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_h4sio4   = x[I_C_H4SiO4] ;
    double c_cah2sio4 = x[I_C_CaH2SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_cahso4   = x[I_C_CaHSO4] ;
    double c_caso4aq  = x[I_C_CaSO4aq] ;
    double c_caoh     = x[I_C_CaOH] ;
    double s_ch       = x[I_S_CH] ;
    double zc_h2so4   = c_h2so4/c_h2so4_eq ;
    
    /* charge density */
    double c_q        = x[I_N_Q] ;
    /* solid contents */
    double n_ch       = x[I_N_CH] ;
    double n_csh2     = x[I_N_CSH2] ;
    double n_si_s     = x[I_N_Si_S]  ;
    double w_s        = W_S ;
    double w_ca       = W_Ca ;
    double p_csh2     = x[I_P_CSH2] ;
    
    /* porosity */
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double phi        = x[I_PHI] ;

    double psi        = x[I_PSI] ;
    double ph         = 14 + log(c_oh)/log(10.) ;
    double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
    
    double Q_CH       = IonActivityProductOfCH(s_ch) ;
    double test       = Q_CH/(c_oh*c_oh) ;
    double damage     = x[I_D_CON] * 1234567;

    i = 0 ;
    Result_Store(r + i++,&ph,"ph",1) ;
    Result_Store(r + i++,&c_h2so4,"c_h2so4",1) ;
    Result_Store(r + i++,&c_hso4,"c_hso4",1) ;
    Result_Store(r + i++,&c_so4,"c_so4",1) ;
    Result_Store(r + i++,&c_ca,"c_ca",1) ;
    Result_Store(r + i++,&c_caoh,"c_caoh",1) ;
    Result_Store(r + i++,&c_h2sio4,"c_h2sio4",1) ;
    Result_Store(r + i++,&c_h3sio4,"c_h3sio4",1) ;
    Result_Store(r + i++,&c_h4sio4,"c_h4sio4",1) ;
    Result_Store(r + i++,&c_cah2sio4,"c_cah2sio4",1) ;
    Result_Store(r + i++,&c_cah3sio4,"c_cah3sio4",1) ;
    Result_Store(r + i++,&c_caso4aq,"c_caso4aq",1) ;
    Result_Store(r + i++,&c_cahso4,"c_cahso4",1) ;
    Result_Store(r + i++,&c_k,"c_k",1) ;
    Result_Store(r + i++,&c_cl,"c_cl",1) ;
    Result_Store(r + i++,&c_oh,"c_oh",1) ;
    Result_Store(r + i++,&c_h,"c_h",1) ;
    Result_Store(r + i++,&n_ch,"n_ch",1) ;
    Result_Store(r + i++,&n_csh2,"n_csh2",1) ;
    Result_Store(r + i++,&n_si_s,"n_si_s",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&psi,"potentiel_electrique",1) ;
    Result_Store(r + i++,&c_q,"charge",1) ;
    Result_Store(r + i++,&zn_ca_s,"zn_ca_s",1) ;
    Result_Store(r + i++,&zc_h2so4,"zc_h2so4",1) ;
    Result_Store(r + i++,&s_ch,"s_ch",1) ;
    Result_Store(r + i++,&zn_si_s,"zn_si_s",1) ;
    Result_Store(r + i++,&v_csh,"V_CSH",1) ;
    Result_Store(r + i++,&x_csh,"C/S",1) ;
    Result_Store(r + i++,&w_s,"W_S",1) ;
    Result_Store(r + i++,&w_ca,"W_Ca",1) ;
    Result_Store(r + i++,&p_csh2,"P_CSH2",1) ;
    Result_Store(r + i++,&damage,"Damage",1) ;
    Result_Store(r + i++,&test,"TEST",1) ;
  }
  
  if(i != nso) arret("ComputeOutputs (Yuan2)") ;

  return(nso) ;
}


//void ComputeTransferCoefficients(PCM_t* pcm)
void ComputeTransferCoefficients(Element_t* el)
/* Termes explicites (va)  */
{
  //Element_t* el = PCM_GetElement(pcm) ;
  double *va = Element_GetExplicitTerm(el) ;
  int i ;
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  /*
    Coefficients de transfert
  */
  for(i = 0 ; i < 2 ; i++) {
    /* molarities */
    /* double *x = ComputeComponents(pcm,0,i) ; */
    //double* x = PCM_GetPointerToComponents(pcm)[i] ;
    double* x = Components[i] ;

    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hso4     = x[I_C_HSO4] ;
    double c_so4      = x[I_C_SO4] ;
    double c_ca       = x[I_C_Ca] ;
    double c_cahso4   = x[I_C_CaHSO4] ;
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_caoh     = x[I_C_CaOH] ;
    double c_k        = x[I_C_K] ;
    double c_cl       = x[I_C_Cl] ;

    /* porosity */
    double phi        = x[I_PHI] ;
    
    /* tortuosite liquide */
    double iff    = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_H2SO4      += d_h2so4*iff;
    KF_HSO4       += d_hso4*iff;
    KF_SO4        += d_so4*iff;

    KF_Ca         += d_ca*iff ;
    KF_CaHSO4     += d_cahso4*iff;
    KF_CaH3SiO4   += d_cah3sio4*iff;

    KF_H3SiO4     += d_h3sio4*iff ;
    KF_H4SiO4     += d_h4sio4*iff ;
    KF_H2SiO4     += d_h2sio4*iff ;
    KF_CaH2SiO4   += d_cah2sio4*iff ;
    KF_CaSO4aq    += d_caso4aq*iff;
    KF_CaOH       += d_caoh*iff ;
    
    KF_K          += d_k*iff;
    KF_Cl         += d_cl*iff;


    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HSO4     += FsRT*KF_HSO4*z_hso4*c_hso4 ;
    Kpsi_SO4      += FsRT*KF_SO4*z_so4*c_so4 ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaHSO4   += FsRT*KF_CaHSO4*z_cahso4*c_cahso4 ;
    Kpsi_CaH3SiO4 += FsRT*KF_CaH3SiO4*z_cah3sio4*c_cah3sio4 ;
    Kpsi_H3SiO4   += FsRT*KF_H3SiO4*z_h3sio4*c_h3sio4 ;
    Kpsi_H2SiO4   += FsRT*KF_H2SiO4*z_h2sio4*c_h2sio4 ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    
    Kpsi_K        += FsRT*KF_K*z_k*c_k ;
    Kpsi_Cl       += FsRT*KF_Cl*z_cl*c_cl ;
    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hso4*Kpsi_HSO4 + z_so4*Kpsi_SO4 + z_ca*Kpsi_Ca + z_cahso4*Kpsi_CaHSO4 + z_h3sio4*Kpsi_H3SiO4 \
                   + z_cah3sio4*Kpsi_CaH3SiO4 + z_caoh*Kpsi_CaOH + z_h2sio4*Kpsi_H2SiO4 + z_k*Kpsi_K + z_cl*Kpsi_Cl ;
    /*KD_CSH2       += n_csh2*k_g;*/
  }
  
  
  {
    FVM_t *fvm   = FVM_GetInstance(el) ;
    double* dist = FVM_ComputeIntercellDistances(fvm) ;
    double dij   = dist[1] ;
    double lij = 0.5/dij ;
  
    /* Averaging */
    for(i = 0 ; i < NVE ; i++) va[i] *= lij ;
  }
}


//double* ComputeFluxes(PCM_t* pcm,double* grd,int i,int j)
double* ComputeFluxes(Element_t* el,double* grd,int i,int j)
{
  //Element_t* el = PCM_GetElement(pcm) ;
  double *va = Element_GetExplicitTerm(el) ;
  double* w  = ComponentFluxes ;
  /* double* w  = PCM_GetPointerToComponentFluxes(pcm)[i] ; */

  /* Gradients */
  double grd_h        = grd[I_C_H] ;
  double grd_oh       = grd[I_C_OH] ;
  
  double grd_h2so4    = grd[I_C_H2SO4] ;
  double grd_hso4     = grd[I_C_HSO4] ;
  double grd_so4      = grd[I_C_SO4] ;
  
  double grd_ca       = grd[I_C_Ca] ;
  double grd_caoh     = grd[I_C_CaOH] ;
  double grd_cahso4   = grd[I_C_CaHSO4] ;
  double grd_caso4aq  = grd[I_C_CaSO4aq] ;
  double grd_cah3sio4 = grd[I_C_CaH3SiO4] ;
  double grd_cah2sio4 = grd[I_C_CaH2SiO4] ;
  
  double grd_h2sio4   = grd[I_C_H2SiO4] ;
  double grd_h3sio4   = grd[I_C_H3SiO4] ;
  double grd_h4sio4   = grd[I_C_H4SiO4] ;
  
  double grd_k        = grd[I_C_K] ;
  
  double grd_cl       = grd[I_C_Cl] ;
  
  double grd_psi      = grd[I_PSI] ;
    
    
  /* Flux */
  double w_h2so4    = - KF_H2SO4*grd_h2so4  ;
  double w_hso4     = - KF_HSO4*grd_hso4          - Kpsi_HSO4*grd_psi ;
  double w_so4      = - KF_SO4*grd_so4            - Kpsi_SO4*grd_psi  ;
  double w_cahso4   = - KF_CaHSO4*grd_cahso4      - Kpsi_CaHSO4*grd_psi ;
  double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
  double w_cah3sio4 = - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
  double w_h3sio4   = - KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi ;
  double w_h2sio4   = - KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi ;
  double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi ;
  double w_h4sio4   = - KF_H4SiO4*grd_h4sio4 ;
  double w_cah2sio4 = - KF_CaH2SiO4*grd_cah2sio4 ;
  double w_caso4aq  = - KF_CaSO4aq*grd_caso4aq ;    
  double w_k        = - KF_K*grd_k                - Kpsi_K*grd_psi ;
  double w_cl       = - KF_Cl*grd_cl              - Kpsi_Cl*grd_psi ; 
    
  double w_q        = - z_h*KF_H*grd_h		      \
                      - z_oh*KF_OH*grd_oh		      \
                      - z_hso4*KF_HSO4*grd_hso4             \
                      - z_so4*KF_SO4*grd_so4		      \
                      - z_ca*KF_Ca*grd_ca		      \
                      - z_cahso4*KF_CaHSO4*grd_cahso4	      \
                      - z_h3sio4*KF_H3SiO4*grd_h3sio4	      \
                      - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                      - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                      - z_caoh*KF_CaOH*grd_caoh \
                      - z_k*KF_K*grd_k \
                      - z_cl*KF_Cl*grd_cl \
                      - Kpsi_q*grd_psi ;
   
  w[I_W_S ]  = w_h2so4 + w_hso4 + w_so4 + w_cahso4 + w_caso4aq ;
  w[I_W_Ca]  = w_ca + w_cahso4 + w_cah3sio4 + w_caso4aq + w_caoh + w_cah2sio4 ;
  w[I_W_Si]  = w_h3sio4 + w_h4sio4 + w_cah3sio4 + w_cah2sio4 + w_h2sio4 ;
  w[I_W_q ]  = w_q ;
  w[I_W_K ]  = w_k ;
  w[I_W_Cl]  = w_cl;
   
  return(w) ;
}



int TangentCoefficients(Element_t* el,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  //PCM_t* pcm = PCM_GetInstance(el,u,u_n,f_n) ;
  Model_t* model = Element_GetModel(el) ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < nn*nn*NEQ*NEQ ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < 2 ; i++) {
    int j = 1 - i ;
    /* Content terms at node i */
    double *cii = c + (i*nn + i)*NEQ*NEQ ;
    /* Transfer terms from node i to node j */
    double *cij = c + (i*nn + j)*NEQ*NEQ ;
    /* Components */
    //double *x         = ComputeComponents(pcm,dt,i) ;
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    double dxi[NEQ] ;
    int k ;
    
    dxi[U_C_H2SO4  ] = 1.e-6*x[U_C_H2SO4] ; 
    dxi[U_ZN_Si_S  ] = 1.e-8 ; 
    dxi[U_ZN_Ca_S  ] = 1.e-10 ;
    dxi[U_C_K      ] = 1.e-6 ;
    dxi[U_C_Cl     ] = 1.e-6 ;
    dxi[U_PSI      ] = 1. ;

    /*
    for(k = 0 ; k < NEQ ; k++) {
      dxi[k] =  1.e-2*ObVal_GetValue(obval + k) ;
    }
    * 
    #if (U_H2SO4 == LOG_U)
    dxi[U_C_H2SO4] *=  x[U_C_H2SO4] ;
    #endif
    * 
    */
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      //double *dx    = ComputeComponentDerivatives(pcm,dt,dxk,k,i) ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      //double *dw    = ComputeFluxes(pcm,dx,0,1) ;
      double *dw    = ComputeFluxes(el,dx,0,1) ;
    
      cii[E_S*NEQ    + k] = dx[I_N_S] ;
      cii[E_Ca*NEQ   + k] = dx[I_N_Ca] ;
      cii[E_Si*NEQ   + k] = dx[I_N_Si] ;
      cii[E_K*NEQ    + k] = dx[I_N_K] ;
      cii[E_Cl*NEQ   + k] = dx[I_N_Cl] ;
      
      cij[E_S*NEQ    + k] = - dt*dw[I_W_S] ;
      cij[E_Ca*NEQ   + k] = - dt*dw[I_W_Ca] ;
      cij[E_Si*NEQ   + k] = - dt*dw[I_W_Si] ;
      cij[E_K*NEQ    + k] = - dt*dw[I_W_K] ;
      cij[E_Cl*NEQ   + k] = - dt*dw[I_W_Cl] ;
      cij[E_q*NEQ    + k] = - dt*dw[I_W_q] ;
    }
  }

  return(dec) ;
}


double concentration_oh(double c_h2so4,double s_ch,double s_sh,double c_k,double c_cl)
/* Solve electroneutrality: SUM(z_i c_i) = 0
   for c_h as root of ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* The primary variables are considered as constant */
  double zc_h2so4     = c_h2so4/c_h2so4_eq ;
  /* Ion activity products are constant as well */
  double s_csh2     = SaturationDegreeOfCSH2(zc_h2so4,s_ch) ;
  double Q_CSH2     = IonActivityProductOfCSH2(s_csh2) ;
  double Q_SH       = IonActivityProductOfSH(s_sh) ;
  /* Other constant concentrations */
  double c_h4sio4   = Q_SH ;
  /* Electroneutrality is written as  Sum z_i*x_i = 0, with
   * ------------------------------------------------------------------
   * x_i        =  A_i*(x_h)**n                              : n    z_i
   * ------------------------------------------------------------------
   * 
   * x_h        = x_h                                        : +1   +1
   * x_oh       = k_e/x_h                                    : -1   -1
   * 
   * c_h2so4    = c_h2so4                                    :  0    0
   * c_hso4     = K_hso4*c_h2so4/c_h                         : -1   -1
   * c_so4      = K_so4*c_hso4/c_h                           : -2   -2
   * 
   * c_ca       = Q_CSH2/c_so4                               : +2   +2
   * c_caoh     = K_caoh*c_ca*c_oh ;                         : +1   +1
   * x_caoh2aq  = k_caoh2*Q_CH                               :  0    0
   * 
   * c_caso4aq  = K_caso4aq*c_ca*c_so4                       :  0    0
   * c_cahso4   = K_cahso4*c_ca*c_hso4                       : +1   +1
   * 
   * c_h4sio4   = Q_SH                                       :  0    0
   * c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)                    : -1   -1
   * c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;                   : -2   -2
   * 
   * c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4                   : +1   +1
   * c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;                 :  0    0
   * 
   * x_k        = x_k                                        :  0   +1
   * x_koh      = k_koh*x_k/x_h                              : -1    0
   * 
   * c_cl       = c_cl                                       :  0   -1
   *
   * We compute the A_i for i having non zero z_i
   */
  double A_hso4     = K_hso4*c_h2so4 ;
  double A_so4      = K_so4*A_hso4 ;
  double A_ca       = Q_CSH2/A_so4 ;
  double A_cahso4   = K_cahso4*A_ca*A_hso4 ;
  double A_h3sio4   = c_h4sio4/K_h4sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;
  double A_h2sio4   = K_h2sio4*A_h3sio4*K_h2o ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahso4*A_cahso4 + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_k*c_k + z_cl*c_cl;
  double d = z_oh*K_h2o + z_hso4*A_hso4 + z_h3sio4*A_h3sio4 ;
  double e = z_so4*A_so4 + z_h2sio4*A_h2sio4 ;

  /* a > 0 ; b > 0 ; c > 0 ; d < 0 ; e < 0 */
  double c_h = poly4(a,b,c,d,e) ;
 
  return(K_h2o/c_h) ;
}


double poly4(double a,double b,double c,double d,double e)
/* on resout ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double tol = 1e-4 ;
  double y[5] ;
  double x ;
  
  y[0] = a ;
  y[1] = b ;
  y[2] = c ;
  y[3] = d ;
  y[4] = e ;
  
  {
    int n = Math_ComputePolynomialEquationRoots(y,4) ;
    int i ;
    for(i = 0 ; i < n ; i++) {
      if((x = y[i]) < 1) break ;
    }
  }
  
  y[0] = a ;
  y[1] = b ;
  y[2] = c ;
  y[3] = d ;
  y[4] = e ;
  
  x = Math_PolishPolynomialEquationRoot(y,4,x,tol*x,20) ;
  
  return(x) ;
}


//double* ComputeComponents(PCM_t* pcm,double dt,int n)
double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  //Element_t* el = PCM_GetElement(pcm) ;
  //double** u   = PCM_GetPointerToNodalUnknowns(pcm) ;
  //double* f_n  = PCM_GetPreviousImplicitTerms(pcm) ;
  double *v0   = Element_GetConstantTerm(el) ;
  /* double *x = PCM_GetPointerToComponents(pcm)[n] ; */
  double *x = Components[n] ;
  
  /* Primary Variables */
  x[U_C_H2SO4 ] = C_H2SO4(n) ;
  x[U_ZN_Ca_S ] = ZN_Ca_S(n) ;
  x[U_ZN_Si_S ] = ZN_Si_S(n) ;
  x[U_C_K     ] = C_K(n) ;
  x[U_C_Cl    ] = C_Cl(n) ;  
  x[U_PSI     ] = PSI(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn  ]  = N_CHn(n) ;
  x[I_N_CSH2n]  = N_CSH2n(n) ;
  x[I_D_CONn ]  = D_CONn(n) ;
  x[I_PHIn   ]  = PHIn(n) ;
  x[I_PHI_Cn ]  = PHI_Cn(n) ;

  x[I_V_Cem0  ] = V_Cem0(n) ;
    
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


//void  ComputeSecondaryComponents(PCM_t* pcm,double dt,double* x)
void  ComputeSecondaryComponents(Element_t *el,double dt,double *x)
{
  //Element_t* el = PCM_GetElement(pcm) ;
  double c_h2so4    = x[U_C_H2SO4] ;
  double zn_si_s    = x[U_ZN_Si_S] ;
  double zn_ca_s    = x[U_ZN_Ca_S] ;
  double c_k        = x[U_C_K] ;
  double c_cl       = x[U_C_Cl] ;      
    
    
  /* Liquid components */
  double zc_h2so4   = c_h2so4/c_h2so4_eq ;

  double s_ch       = SaturationDegreeOfCH(zc_h2so4,zn_ca_s) ;
  double s_sh       = SaturationDegreeOfSH(s_ch,zn_si_s) ;
  double Q_SH       = IonActivityProductOfSH(s_sh) ;
  double s_csh2     = SaturationDegreeOfCSH2(zc_h2so4,s_ch) ;
  double Q_CSH2     = IonActivityProductOfCSH2(s_csh2) ;
  
  double c_h4sio4   = Q_SH ;
  double c_oh       = ConcentrationOfOHInLiquid(c_h2so4,s_ch,s_sh,c_k,c_cl) ;
  double c_h        = K_h2o/c_oh ;
  double c_hso4     = K_hso4*c_h2so4/c_h ;
  double c_so4      = K_so4*c_hso4/c_h ;
  double c_ca       = Q_CSH2/c_so4 ;
  double c_cahso4   = K_cahso4*c_ca*c_hso4 ;
  double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
  double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
  double c_caso4aq  = K_caso4aq*c_ca*c_so4 ;
  double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
  double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
  double c_caoh     = K_caoh*c_ca*c_oh ;
  
  double c_q     = z_h*c_h + z_oh*c_oh \
                 + z_hso4*c_hso4 + z_so4*c_so4 \
                 + z_ca*c_ca + z_caoh*c_caoh \
                 + z_h3sio4*c_h3sio4 + z_h2sio4*c_h2sio4 \
                 + z_cah3sio4*c_cah3sio4 + z_cahso4*c_cahso4 \
                 + z_k*c_k + z_cl*c_cl;
  
  
  /* Solid contents */
  /* ... as components: CH, CSH2, CSH */
  double n_chn      = x[I_N_CHn] ;
  double n_csh2n    = x[I_N_CSH2n] ;
  double n_ch_ci    = n_chn*pow(zc_h2so4,-dt/t_ch) ;
  double n_csh2_ci  = n_csh2n*pow(zc_h2so4,dt/t_csh2) ;
  double n_ch_csh2  = CalciumContentInCHAndCSH2(zn_ca_s) ;
  double n_ch       = (zc_h2so4 <= 1) ? n_ch_csh2 - n_csh2_ci : n_ch_ci ;
  double n_csh2     = (zc_h2so4 >  1) ? n_ch_csh2 - n_ch_ci : n_csh2_ci ;
  /* ... as elements: S, Ca, Si */
  double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double n_si_s     = SiliconContentInCSH(zn_si_s) ;
  double n_ca_s     = n_ch_csh2 + x_csh*n_si_s ;
  double n_s_s      = n_csh2 ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(s_ch) ;
  double v_cem      = V_CH*n_ch + v_csh*n_si_s ;
  double v_gyp      = V_Gyp*n_csh2 ;
  double v_csh2     = V_CSH2*n_csh2 ;


  /* Porosities in unconfined conditions (no pressure) */
  double v_cem0     = x[I_V_Cem0] ;
  /* ... of concrete */
  double phi_con    = phi0 + v_cem0 - v_cem ;
  /* ... of gypsum */
  double phi_gyp    = PHI_Gyp ;


  /* Damage modeling */
  /* 
   * Confining pressure of gypsum calculated from the equality of:
   * volume of gypsum:  v_gyp(p)   = v_gyp   * ( 1 - p/CC_Gyp ) 
   * and
   * concrete porosity: phi_con(p) = phi_con * ( 1 + p/CC_Con ) 
   */
  double p_csh20    = (v_gyp - phi_con)/(phi_con/CC_Con + v_gyp/CC_Gyp) ;
  double p_csh2     = (v_gyp > phi_con) ? p_csh20 : 0 ;
  
  /* Saturation degree of gypsum */
  double s_gyp     = (v_gyp > phi_con) ? 1 : v_gyp/phi_con ;
    
  /* Damage */
  double damagen    = x[I_D_CONn] ;
  double damage     = (p_csh2 > CS_Con || damagen > 0.5) ? 1 : 0 ;
  
  if(damage > 0.5) {
    p_csh2    = 0 ;
    /* phi_con   = phi_gyp ; */
  }
    
  if(p_csh2 > CC_Gyp*phi_gyp) {
    arret("ComputeSecondaryComponents: pressure too high") ;
  }

  /* Porosities in confined conditions */
  /* ... of concrete */
  double phi_conp      = phi_con * (1 + p_csh2/CC_Con) ;
  /* ... of gypsum */
  /* ERROR FOUND in original M70con !!! */
  /* double phi_gypp    = 1 - (1 - phi_gyp) * (1 - p_csh2/CC_Gyp) ; */
  double phi_gypp    = (phi_gyp - p_csh2/CC_Gyp) / (1 - p_csh2/CC_Gyp) ;
  
  /* Volume of gypsum */
  double v_gypp        = v_gyp * (1 - p_csh2/CC_Gyp) ;
    
  /* Total porosity */
  /* ... in case of filled pores */
  /* double phi_fil = phi_conp * phi_gypp ; */
  double phi_fil = phi_con * phi_gyp ;
  /* double phi_fil = v_gypp - v_csh2 ; */
  /* ... in case of unfilled pores */
  /* double phi_ufil = phi_conp - v_csh2 ; */
  double phi_ufil = phi_con - v_csh2 ;
  /* ... in case of undamaged sample */
  double phi_udam = phi_conp - v_csh2 ;
  /* ... in case of damaged sample */
  double phi_dam = phi_con * phi_gyp ;

  double phi_t      = (v_gyp < phi_con) ? phi_ufil : phi_fil ;
  /* double phi_t      = (damage > 0.5) ? phi_dam : phi_udam ; */
  double phi_c      = phi_conp ;
  

#if (U_PHI == IMPLICIT)
  double phi_l        = phi_t ;
#else
  double phi_l        = x[I_PHIn] ;
#endif
    
    
  /* Liquid contents */
  /* ... as elements: S, Ca, Si, K, Cl */
  double c_s_l  = c_h2so4 + c_hso4 + c_so4 + c_cahso4 + c_caso4aq ;
  double c_ca_l = c_ca + c_cahso4 + c_cah3sio4 + c_cah2sio4 + c_caso4aq + c_caoh ;
  double c_si_l = c_h3sio4 + c_h4sio4 + c_cah3sio4 + c_h2sio4 + c_cah2sio4 ;
  double c_k_l  = c_k ;
  double c_cl_l = c_cl ;
  double n_s_l  = phi_l*c_s_l ;
  double n_ca_l = phi_l*c_ca_l ;
  double n_si_l = phi_l*c_si_l ;
  double n_k_l  = phi_l*c_k_l ;
  double n_cl_l = phi_l*c_cl_l ;

  

  /* Back up */
  
  /* Liquid components */
  x[I_C_OH      ] = c_oh ;
  x[I_C_H       ] = c_h ;
  
  x[I_C_H2SO4   ] = c_h2so4 ;
  x[I_C_HSO4    ] = c_hso4 ;
  x[I_C_SO4     ] = c_so4 ;
  
  x[I_C_Ca      ] = c_ca ;
  x[I_C_CaOH    ] = c_caoh ;
  
  x[I_C_H4SiO4  ] = c_h4sio4 ;
  x[I_C_H3SiO4  ] = c_h3sio4 ;
  x[I_C_H2SiO4  ] = c_h2sio4 ;
  
  x[I_C_CaH3SiO4] = c_cah3sio4 ;
  x[I_C_CaH2SiO4] = c_cah2sio4 ;
  
  x[I_C_CaHSO4  ] = c_cahso4 ;
  x[I_C_CaSO4aq ] = c_caso4aq ;
  
  x[I_C_K       ] = c_k ;
  
  x[I_C_Cl      ] = c_cl ;    
  
  x[I_S_CH      ] = s_ch ;


  /* Solid components */
  x[I_N_CH      ] = n_ch ;
  x[I_N_CSH2    ] = n_csh2 ;
  x[I_N_Si_S    ] = n_si_s ;
  x[I_P_CSH2    ] = p_csh2 ;
  x[I_D_CON     ] = damage ;
  
  x[I_ZN_Si_S   ] = zn_si_s ;
  x[I_ZN_Ca_S   ] = zn_ca_s ;  
  
  
  /* Porosities */
  x[I_PHI       ] = phi_t ;
  x[I_PHI_C     ] = phi_c ;
  
  
  /* Element contents */
  x[I_N_S       ] = n_s_l  + n_s_s ;
  x[I_N_Ca      ] = n_ca_l + n_ca_s ;
  x[I_N_Si      ] = n_si_l + n_si_s ;
  x[I_N_K       ] = n_k_l  ;
  x[I_N_Cl      ] = n_cl_l  ;

  
  /* Charge density */
  x[I_N_Q       ]  = c_q ;
  
  /* Electric potential */
  x[I_PSI       ]  = x[U_PSI] ;
}


//double* ComputeComponentDerivatives1(PCM_t* pcm,double dt,double dxi,int i,int n)
double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  //double* x  = PCM_GetPointerToComponents(pcm)[n] ;
  /* double* dx = PCM_GetPointerToComponentDerivatives(pcm)[n] ; */
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  dx[U_C_H2SO4 ] = x[U_C_H2SO4 ] ;
  dx[U_ZN_Ca_S ] = x[U_ZN_Ca_S ] ;
  dx[U_ZN_Si_S ] = x[U_ZN_Si_S ] ;
  dx[U_C_K     ] = x[U_C_K     ] ;
  dx[U_C_Cl    ] = x[U_C_Cl    ] ;
  dx[U_PSI     ] = x[U_PSI     ] ;
  
  /* Needed variables to compute secondary components */
  dx[I_N_CHn   ] = x[I_N_CHn  ] ;
  dx[I_N_CSH2n ] = x[I_N_CSH2n] ;
  dx[I_D_CONn  ] = x[I_D_CONn ] ;
  dx[I_PHIn    ] = x[I_PHIn   ] ;
  dx[I_PHI_Cn  ] = x[I_PHI_Cn ] ;

  dx[I_V_Cem0  ] = x[I_V_Cem0  ] ;
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryComponents(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfComponents ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}




//double* ComputeComponentFluxes1(PCM_t* pcm,int i,int j)
double* ComputeComponentFluxes(Element_t* el,int i,int j)
{
  //Element_t* el = PCM_GetElement(pcm) ;
  //Model_t* model = Element_GetModel(el) ;
  /* double* grdij = PCM_GetPointerToComponentDerivatives(pcm)[i] ; */
  double *grdij = dComponents ;


  {
    //int NbOfComponents = Model_GetNbOfVariables(model) ;
  
    /*
    if(NbOfComponents > PCM_MaxNbOfComponents) {
      arret("PCM_ComputeComponentFluxes") ;
    }
    */
    
    {
      //double* xi  = PCM_GetPointerToComponents(pcm)[i] ;
      //double* xj  = PCM_GetPointerToComponents(pcm)[j] ;
      double* xi  = Components[i] ;
      double* xj  = Components[j] ;
      int k ;
      
      for(k = 0 ; k < NbOfComponents ; k++)  {
        grdij[k] = xj[k] - xi[k] ;
      }
    }
  }
  
  /* Fluxes */
  {
    //int NbOfComponentFluxes = Model_GetNbOfVariableFluxes(model) ;
  
    /*
    if(NbOfComponentFluxes > PCM_MaxNbOfComponentFluxes) {
      arret("PCM_ComputeComponentFluxes") ;
    }
    */
    
    {
      double* w = ComputeFluxes(el,grdij,i,j) ;
    
      return(w) ;
    }
  }
}
