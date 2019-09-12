/* General features of the model:
 * To be completed!!!!!
 * This model is from M70congas
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the numerical method */
#include "FVM.h"

#define TITLE   "Hydrogen Sulfide attack of concrete (Jan. 2015)" 
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
#define I_C_H2S     (0)
#define I_PSI       (1)
#define I_ZN_Ca_sol (2)
#define I_ZN_Si_sol (3)
#define I_C_K       (4)
#define I_C_Cl      (5)


#define NOLOG_U     1
#define LOG_U       2
#define Ln10        (2.302585093)
#define U_H2S       LOG_U


/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWN_n(n,i)   (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])



#if (U_H2S == LOG_U)
  #define LogC_H2S(n)   (UNKNOWN(n,I_C_H2S))
  #define LogC_H2Sn(n)  (UNKNOWN_n(n,I_C_H2S))
  #define C_H2S(n)      (pow(10,UNKNOWN(n,I_C_H2S)))
  #define C_H2Sn(n)     (pow(10,UNKNOWN_n(n,I_C_H2S)))
#else
  #define C_H2S(n)      (UNKNOWN(n,I_C_H2S))
  #define C_H2Sn(n)     (UNKNOWN_n(n,I_C_H2S))
  #define LogC_H2S(n)   (log10(UNKNOWN(n,I_C_H2S)))
  #define LogC_H2Sn(n)  (log10(UNKNOWN_n(n,I_C_H2S)))
#endif
#define ZN_Ca_sol(n) (UNKNOWN(n,I_ZN_Ca_sol))
#define ZN_Si_sol(n) (UNKNOWN(n,I_ZN_Si_sol))
#define PSI(n)       (UNKNOWN(n,I_PSI))
#define C_K(n)       (UNKNOWN(n,I_C_K))
#define C_Cl(n)      (UNKNOWN(n,I_C_Cl))

#define ZN_Ca_soln(n) (UNKNOWN_n(n,I_ZN_Ca_sol))
#define ZN_Si_soln(n) (UNKNOWN_n(n,I_ZN_Si_sol))
#define PSIn(n)       (UNKNOWN_n(n,I_PSI))
#define C_Kn(n)       (UNKNOWN_n(n,I_C_K))
#define C_Cln(n)      (UNKNOWN_n(n,I_C_Cl))


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
#define N_CaS(n)   (f[(20+n)])

#define N_Sn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_Cln(n)   (f_n[(10+n)])
#define N_CHn(n)   (f_n[(18+n)])
#define N_CaSn(n)  (f_n[(20+n)])


#define KF_OH         (va[(0)])
#define KF_H          (va[(1)])
#define KF_CO2        (va[(2)])
#define KF_H2S        (va[(3)])
#define KF_HS         (va[(4)])
#define KF_S          (va[(5)])
#define KF_Ca         (va[(6)])
#define KF_CaHS       (va[(7)])
#define KF_CaH3SiO4   (va[(8)])
#define KF_H3SiO4     (va[(9)])
#define KF_H4SiO4     (va[(10)])
#define KF_H2SiO4     (va[(11)])
#define KF_CaH2SiO4   (va[(12)])
#define KF_CaSaq      (va[(13)])
#define KF_CaOH       (va[(14)])
#define KF_K          (va[(15)])
#define KF_Cl         (va[(16)])

#define Kpsi_OH       (va[(17)])
#define Kpsi_H        (va[(18)])
#define Kpsi_HS       (va[(19)])
#define Kpsi_S        (va[(20)])
#define Kpsi_Ca       (va[(21)])
#define Kpsi_CaHS     (va[(22)])
#define Kpsi_CaH3SiO4 (va[(23)])
#define Kpsi_H3SiO4   (va[(24)])
#define Kpsi_q        (va[(25)])
#define Kpsi_H2SiO4   (va[(26)])
#define Kpsi_CaOH     (va[(27)])
#define Kpsi_K        (va[(28)])
#define Kpsi_Cl       (va[(29)])
#define KD_CaS        (va[(30)])


#define V_Mat0(n)     (v0[(0+n)])

/*
  Aqueous solution
*/

/* Charge of ions */
#include "ElectricChargeOfIonInWater.h"
#define z_ca          ElectricChargeOfIonInWater(Ca)
#define z_h           ElectricChargeOfIonInWater(H)
#define z_oh          ElectricChargeOfIonInWater(OH)
#define z_hco3        ElectricChargeOfIonInWater(HCO3)
#define z_co3         ElectricChargeOfIonInWater(CO3)
#define z_h3sio4      ElectricChargeOfIonInWater(H3SiO4)
#define z_cahco3      ElectricChargeOfIonInWater(CaHCO3)
#define z_cah3sio4    ElectricChargeOfIonInWater(CaH3SiO4)
#define z_h2sio4      ElectricChargeOfIonInWater(H2SiO4)
#define z_caoh        ElectricChargeOfIonInWater(CaOH)
#define z_k           ElectricChargeOfIonInWater(K)
#define z_cl          ElectricChargeOfIonInWater(Cl)
#define z_hs          ElectricChargeOfIonInWater(HS)
#define z_s           ElectricChargeOfIonInWater(S)
#define z_cahs        ElectricChargeOfIonInWater(CaHS)

/* volumes molaires partiels des ions (dm3/mole) from [Millero F J,Partial molar volum of ions in seawater]*/
#ifdef NOTDEFINED
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_h2o      (18.e-3)
#define v_h2s      (34.9e-3)
#define v_hs       (27.3e-3)   /* d'apres Lothenbach */
#define v_s        (-0.2e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)
#define v_casaq    (26.20e-3)    /* a modifier */
#define v_caoh     (26.20e-3)    /* a modifier */
#define v_k        (43.93e-3)    /* d'apres Antoine */
#define v_cl       (43.93e-3)    /*a modifier*/
#define v_koh      (27.44e-3)    /* d'apres Antoine */
#endif



/* Equilibrium constant of homogeneous aqueous reaction */
#define K_h2o      (1.e-14)          /* H2O = OH[-] + H[+] */

#define K_hs       (8.9e-8)          /* H2S(aq) = HS[-] + H[+] */
#define K_s        (1.2e-13)         /* HS[-] = S[2-] + H[+] */

#define K_h2sio4   (4.68)            /* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define K_h4sio4   (6.45e9)          /* H3SiO4[-] + H[+]  = H4SiO4 */
#define K_h3sio4   (1.55e-10)        /* H4SiO4  = H3SiO4[-] + H[+] */

#define K_cahs     (1.276e+1)        /* Ca[2+] + HS[-]      = CaHS[+] */
#define K_casaq    (3.5e+3)          /* Ca[2+] + S[2-]      = CaS[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-]      = CaOH[+] */

/* Dissolution of H2S gas */
#define KHenry_H2S (0.987e-6)        /* Henry's law constant (mol/(L*Pa) */
/* c_h2s = KH*p_h2s where p_h2s is the partial pressure of H2S
 * ie p_h2s = ppm*0.101325, hence c_h2s = 1.e-7*ppm (mol/L) */
#define KF_A_H2S   (2e-9)            /* mol/dm2 s */



/*
  Solids
  CH  = Calcium Hydroxide (Portlandite)
  CaS = Calcium Sulfide
  CSH = Calcium Silicates Hydrate
  SH  = Amorphous Silica Gel
*/

/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(q) (Curve_ComputeValue(Element_GetCurve(el),q))
#define WaterSiliconRatioInCSH(q)   (Curve_ComputeValue(Element_GetCurve(el) + 1,q))
#define MolarVolumeOfCSH(q)         (Curve_ComputeValue(Element_GetCurve(el) + 2,q))


/* S-H Properties */
/* Equilibrium constant */
#define K_SH       (1.93642e-3)      /* SHt = H4SiO4 + (t-2)H2O */
/* Saturation Degree of Dissolved S-H */
#define S_SHeq(q)     (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define SaturationDegreeOfSH(s_ch,zn_si_sol)       (NEGEXP(zn_si_sol)*S_SHeq(s_ch))
/* Ion Activity Product of Dissolved S-H */
#define IonActivityProductOfSH(s_ch,zn_si_sol)     (K_SH*SaturationDegreeOfSH(s_ch,zn_si_sol))


/* CH Properties */
/* Equilibrium constant */
#define K_CH       (6.456e-6)        /* CH  = Ca[2+] + 2OH[-] */
/* Molar volume of CH solid (dm3/mole) */
#define V_CH       (33.e-3)      /* (33.e-3) */
/* Saturation Degree of Dissolved CH */
#define SaturationDegreeOfCH(zc_h2s,zn_ca_sol)    (NEGEXP(zn_ca_sol)/MAX(zc_h2s,1.))
/* Ion Activity Product of Dissolved CH */
#define IonActivityProductOfCH(zc_h2s,zn_ca_sol)  (K_CH*SaturationDegreeOfCH(zc_h2s,zn_ca_sol))


/* CaS Properties */
/* Equilibrium constant */
#define K_CaS     (7.9e-7)         /* CaS = Ca[2+] + S[2-] */
/*  Threshlod */
#define C_H2S_eq (K_h2o*K_h2o*K_CaS/(K_hs*K_s*K_CH))
/* Molar volume of CaS solid (dm3/mole) */
#define V_CaS     (25.e-3)
/* Saturation Degree of Dissolved CaS */
#define SaturationDegreeOfCaS(zc_h2s,zn_ca_sol)   (NEGEXP(zn_ca_sol)*MIN(zc_h2s,1.))
/* Ion Activity Product of Dissolved CaS */
#define IonActivityProductOfCaS(zc_h2s,zn_ca_sol) (K_CaS*SaturationDegreeOfCaS(zc_h2s,zn_ca_sol))


/* Element contents in solid phases  */
#define CalciumContentInCHAndCaS(zn_ca_sol)  (n_ca_ref*MAX(zn_ca_sol,0.))
#define SiliciumContentInCSH(zn_si_sol)      (n_si_ref*MAX(zn_si_sol,0.))


/* Concentration of OH computed from electroneutrality */
#define ConcentrationOfOHInLiquid(A,B,C,D,E)  concentration_oh(A,el,B,C,D,E)


/* Function used above */
#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])



/* Internal Functions */
static int    pm(const char *s) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void    ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;

static double concentration_oh(double,elem_t*,double,double,double,double) ;
static double poly4(double,double,double,double,double) ;


/* Parameters */
static double phi0,c_h2s_eq,t_ch,t_cas ;
static double n_ca_ref,n_si_ref ;

static double d_h,d_oh ;
static double d_ca,d_caoh,d_caoh2aq ;
static double d_h2s,d_hs,d_s ;
static double d_cahs,d_casaq ;
static double d_h4sio4,d_h3sio4,d_h2sio4 ;
static double d_cah2sio4,d_cah3sio4 ;
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
  
  d_h2s        = (1.36e-7) ;  /* (radius = 1.5e-10 m) */
  d_hs         = (1.73e-7) ;  /* (radius = 1.91e-10 m) */
  d_s          = (1.13e-7) ;  /*? (radius = 1.89e-10 m) */
  
  d_ca         = DiffusionCoefficientOfMoleculeInWater(Ca,TK)*dm2 ;
  d_caoh       = DiffusionCoefficientOfMoleculeInWater(CaOH,TK)*dm2 ;
  d_caoh2aq  	 = 7.92e-8 ;
  
  d_cahs       = (1.07e-7) ;   /* (radius = 2e-10 m) */
  d_casaq      = (1.43e-7) ;   /* (radius = 1.5e-10 m) */
  
  d_h4sio4     = DiffusionCoefficientOfMoleculeInWater(H4SiO4,TK)*dm2 ;
  d_h3sio4     = DiffusionCoefficientOfMoleculeInWater(H3SiO4,TK)*dm2 ;
  d_h2sio4     = DiffusionCoefficientOfMoleculeInWater(H2SiO4,TK)*dm2 ;
  
  d_cah2sio4   = DiffusionCoefficientOfMoleculeInWater(CaH2SiO4,TK)*dm2 ;
  d_cah3sio4   = DiffusionCoefficientOfMoleculeInWater(CaH3SiO4,TK)*dm2 ;
  
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


#define NbOfComponents    (34)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;

#define I_C_OH         (6)
#define I_C_H          (7)

#define I_C_HS         (8)
#define I_C_S          (9)

#define I_C_Ca         (10)
#define I_C_CaOH       (11)
#define I_C_CaHS       (12)
#define I_C_CaSaq      (13)

#define I_C_H2SiO4     (14)
#define I_C_H3SiO4     (15)
#define I_C_H4SiO4     (16)

#define I_C_CaH2SiO4   (17)
#define I_C_CaH3SiO4   (18)

#define I_S_CH         (19)

#define I_N_Q          (20)

#define I_N_S          (21)
#define I_N_Ca         (22)
#define I_N_Si         (23)
#define I_N_K          (24)
#define I_N_Cl         (25)

#define I_N_CH         (26)
#define I_N_CaS        (27)
#define I_V_Mat        (28)

#define I_N_CHn        (29)
#define I_N_CaSn       (30)
#define I_V_Mat0       (31)

#define I_Phi          (32)

#define I_N_Si_sol     (33)


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
  if(strcmp(s,"porosite") == 0)       return (0) ;
  else if(strcmp(s,"N_CH") == 0)      return (1) ;
  else if(strcmp(s,"N_Si") == 0)      return (2) ;
  else if(strcmp(s,"C_H2S_eq") == 0) return (3) ;
  else if(strcmp(s,"T_CH") == 0)      return (4) ;
  else if(strcmp(s,"T_CaS") == 0)     return (5) ;
  else return(-1) ;
}



int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_S, "sulfur"   ) ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium"  ) ;
  Model_CopyNameOfEquation(model,E_Si,"silicium" ) ;
  Model_CopyNameOfEquation(model,E_q, "charge"   ) ;
  Model_CopyNameOfEquation(model,E_K, "potassium") ;
  Model_CopyNameOfEquation(model,E_Cl,"chlorine" ) ;

  /** Names of the main (nodal) unknowns */
#if (U_H2S == LOG_U)
  Model_CopyNameOfUnknown(model,I_C_H2S,"logc_h2s") ;
#else
  Model_CopyNameOfUnknown(model,I_C_H2S,"c_h2s"   ) ;
#endif
  Model_CopyNameOfUnknown(model,I_ZN_Ca_sol,"z_ca"    ) ;
  Model_CopyNameOfUnknown(model,I_PSI,    "psi"     ) ;
  Model_CopyNameOfUnknown(model,I_ZN_Si_sol,"z_si"    ) ;
  Model_CopyNameOfUnknown(model,I_C_K,    "c_k"     ) ;
  Model_CopyNameOfUnknown(model,I_C_Cl,   "c_cl"    ) ;

  {
    double temperature = 293 ;
    
    ComputePhysicoChemicalProperties(temperature) ;
  }
  
  return(0) ;
}



int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 6 ;

  {
    /* Self-initialization */
    Material_GetProperty(mat)[pm("N_CH")]  = 1 ;
    Material_GetProperty(mat)[pm("N_Si")]  = 1 ;
    Material_GetProperty(mat)[pm("T_CH")]  = 600 ;
    Material_GetProperty(mat)[pm("T_CaS")] = 600 ;
  
    Material_ScanProperties(mat,datafile,pm) ;

    t_ch     = Material_GetProperty(mat)[pm("T_CH")] ;
    t_cas    = Material_GetProperty(mat)[pm("T_CaS")] ;

    if(t_cas  == 0.) Material_GetProperty(mat)[pm("T_CaS")]  = t_ch ;
    
    Material_GetProperty(mat)[pm("C_H2S_eq")] = C_H2S_eq ;
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
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length : dm !\n") ;
  printf("\t time   : s !\n") ;
  
  printf("\n") ;
  printf("The 6 primary unknowns are:\n") ;
  printf("\t- Hydrogen sulfide concentration  (c_h2s or logc_h2s)\n") ;
  printf("\t- Electric potential              (psi)\n") ;
  printf("\t- Zeta unknown for calcium        (z_ca)\n") ;
  printf("\t- Zeta unknown for silicium       (z_si)\n") ;
  printf("\t- Potassium concentration         (c_k)\n") ;
  printf("\t- Chlorine concentration          (c_cl)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH   = 6.1      # Contenu en CH (moles/L)\n") ;
  fprintf(ficd,"N_Si   = 2.4      # contenu en Si solide (moles/L)\n") ;
  fprintf(ficd,"N_K    = 0.4      # contenu en K (moles/L)\n") ;
  fprintf(ficd,"N_Cl   = 0.4      # contenu en Cl (moles/L)\n") ;
  fprintf(ficd,"T_CH   = 1.e5     # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"T_CaS  = 1.e5     # Cinetique de dissolution de CaS (s)\n") ;
  fprintf(ficd,"Curves = solid    # Nom du fichier: q_ch X Y V f_S\n") ;

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
  char* nameofeq = Load_GetNameOfEquation(cg) ;
  int ieq = Element_FindEquationPositionIndex(el,nameofeq) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
    
    
  if(ieq == E_S && Load_TypeIs(cg,"H2Sflowrate") && nn == 1) {
    int dim = Element_GetDimensionOfSpace(el) ;
    Node_t **no = Element_GetPointerToNode(el) ;
    double *x = Node_GetCoordinate(no[0]) ;
    Field_t *field = Load_GetField(cg) ;
    double p_h2s = Field_ComputeValueAtPoint(field,x,dim) ;
    double n = 0.55 ;
    double a_h2s = KF_A_H2S*pow(fabs(p_h2s),n) ;
    double w_h2s = (p_h2s > 0) ? a_h2s : - a_h2s ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = 0 ;
    
    {
      Function_t *function = Load_GetFunction(cg) ;
      double ft  = (function) ? Function_ComputeValue(function,t) : dt ;
      double ftn = (function) ? Function_ComputeValue(function,t - dt) : 0 ;
      double dft = ft - ftn ;
      
      r[ieq] = - dft*w_h2s ;
    }
    
  } else {
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
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int i ;
  
  
  /*
    Donnees
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2s_eq  = GetProperty("C_H2S_eq") ;
  
  
  /* Pre-initialization */
  for(i = 0 ; i < nn ; i++) {
    double c_h2s       = (C_H2S(i) > 0.) ? C_H2S(i) : c_h2s_eq ;
    double zn_ca_sol   = ZN_Ca_sol(i) ;
    double zn_si_sol   = ZN_Si_sol(i) ;
    
    /* Liquid components */
    double zc_h2s     = c_h2s/c_h2s_eq ;

    double s_ch       = SaturationDegreeOfCH(zc_h2s,zn_ca_sol) ;
    
    /* Solid contents */
    /* ... as components: CH, CaS, CSH */
    double n_ch_cas   = CalciumContentInCHAndCaS(zn_ca_sol) ;
    double n_ch       = (zc_h2s <= 1) ? n_ch_cas  : 0 ;
    double n_cas      = (zc_h2s >  1) ? n_ch_cas  : 0 ;
    /* ... as elements: Ca, Si, S */
    double n_si_sol   = SiliciumContentInCSH(zn_si_sol) ;
    /* ... as volume */
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double v_mat      = V_CH*n_ch + V_CaS*n_cas + v_csh*n_si_sol ;
      
    /* Back up what is needed to compute components */
    N_CH(i)    = n_ch ;
    N_CaS(i)   = n_cas ;
      
    V_Mat0(i)  = v_mat ;

#if (U_H2S == LOG_U)
    LogC_H2S(i) = log10(c_h2s) ;
#else
    C_H2S(i)   = c_h2s ;
#endif
  }
  
  
  for(i = 0 ; i < nn ; i++) {
    /* Components */
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
    N_CaS(i)   = x[I_N_CaS] ;
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  ComputeTransferCoefficients(el,u,f) ;

  /* Flux */
  ComputeFluxes(el,u) ;
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Thermes explicites (va)  */
{
  double *f = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int i ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Coefficients de transfert
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int i ;
  
  /*
    Donnees
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_h2s_eq  = GetProperty("C_H2S_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cas     = GetProperty("T_CaS") ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
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
    N_CaS(i)   = x[I_N_CaS] ;


    {
      double c_h2s     = x[I_C_H2S] ;
      double zn_si_sol = x[I_ZN_Si_sol] ;
      double zn_ca_sol = x[I_ZN_Ca_sol] ;
      double s_ch      = x[I_S_CH] ;
      double c_h3sio4  = x[I_C_H3SiO4] ;
      double c_oh      = x[I_C_OH] ;
      double n_ch      = x[I_N_CH] ;
      double n_cas     = x[I_N_CaS] ;
      double x_csh     = CalciumSiliconRatioInCSH(s_ch) ;
      double n_si_sol  = x[I_N_Si_sol] ;
      double n_ca_sol  = n_ch + n_cas + x_csh*n_si_sol ;
      double phi       = x[I_Phi] ;
    
      if(c_h2s < 0. || n_ca_sol < 0. || n_si_sol < 0.) {
        double xx = Element_GetNodeCoordinate(el,i)[0] ;
        printf("x         = %e\n",xx) ;
        printf("c_h2s    = %e\n",c_h2s) ;
        printf("n_cas    = %e\n",n_cas) ;
        printf("n_ca_sol    = %e\n",n_ca_sol) ;
        printf("n_si_sol    = %e\n",n_si_sol) ;
        printf("zn_si_sol   = %e\n",zn_si_sol) ;
        printf("zn_ca_sol   = %e\n",zn_ca_sol) ;
        printf("c_h3sio4  = %e\n",c_h3sio4) ;
        printf("c_oh      = %e\n",c_oh) ;
        return(-1) ;
      }
 
      if(phi < 0.) {
        double xx = Element_GetNodeCoordinate(el,i)[0] ;
        printf("phi = %e\n",phi) ;
        printf("CH = %e\n",n_ch) ;
        printf("CaS = %e\n",n_cas) ;
        printf("Si = %e\n",n_si_sol) ;
        printf("c_h2s    = %e\n",c_h2s) ;
        printf("x         = %e\n",xx) ;
        return(-1) ;
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Flux */
  ComputeFluxes(el,u) ;

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
  c_h2s_eq = GetProperty("C_H2S_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cas     = GetProperty("T_CaS") ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }
  

#if (U_H2S == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,I_C_H2S)     *= Ln10*C_H2S(0) ;
      K(i,I_C_H2S+NEQ) *= Ln10*C_H2S(1) ;
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
  Symmetry_t sym = Element_GetSymmetry(el) ;
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i=0;i<NEQ*nn;i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    CALCUL DE volume ET DE surf
  */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    dx = x1 - x0 ;
    xm = (x1 + x0)/deux ;
  }
  for(i=0;i<nn;i++) {
    double x = Element_GetNodeCoordinate(el,i)[0] ;
    volume[i] = fabs(dx)/deux ; 
    if(sym == AXIS) volume[i] *= M_PI*(x + xm) ; 
  }
  if(sym == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  /*
    Conservation de S (sulfur)
  */
  R(0,E_S) -= volume[0]*(N_S(0) - N_Sn(0)) + dt*surf*W_S ;
  R(1,E_S) -= volume[1]*(N_S(1) - N_Sn(1)) - dt*surf*W_S ;
  /*
    Conservation de Ca (calcium)
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Conservation de Si (silicium)
  */
  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;
  /*
      Conservation de K (potassium)
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ; 
  
    /*
      Conservation de Cl (chlorine)
  */
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*W_Cl ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*W_Cl ; 
  /*
    Conservation de la charge
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
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    nso = 29 ;
  int    i ;
  

  /* if(Element_IsSubmanifold(el)) return(0) ; */
  
  /*
    Input data
  */
  c_h2s_eq  = GetProperty("C_H2S_eq") ;

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }


  /* Quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    double *x = ComputeComponents(el,u,f,0,j) ;
    
    double c_h2s      = x[I_C_H2S] ;
    double zn_si_sol  = x[I_ZN_Si_sol] ;
    double zn_ca_sol  = x[I_ZN_Ca_sol] ;
    double c_k        = x[I_C_K] ;
    double c_cl       = x[I_C_Cl] ;
    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hs       = x[I_C_HS] ;
    double c_s        = x[I_C_S] ;
    double c_ca       = x[I_C_Ca] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_h4sio4   = x[I_C_H4SiO4] ;
    double c_cah2sio4 = x[I_C_CaH2SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_cahs     = x[I_C_CaHS] ;
    double c_casaq    = x[I_C_CaSaq] ;
    double c_caoh     = x[I_C_CaOH] ;
    double s_ch       = x[I_S_CH] ;
    double zc_h2s     = c_h2s/c_h2s_eq ;
    
    /* charge density */
    double c_q        = x[I_N_Q] ;
    /* solid contents */
    double n_ch       = x[I_N_CH] ;
    double n_cas      = x[I_N_CaS] ;
    double n_si_sol   = x[I_N_Si_sol]  ;

    
    /* porosity */
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double phi        = x[I_Phi] ;

    double psi        = x[I_PSI] ;
    double ph         = 14 + log(c_oh)/log(10.) ;
    double pk_ch      = s_ch ;
    double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;


    i = 0 ;
    Result_Store(r + i++,&ph,"ph",1) ;
    Result_Store(r + i++,&c_h2s,"c_h2s",1) ;
    Result_Store(r + i++,&c_hs,"c_hs",1) ;
    Result_Store(r + i++,&c_s,"c_s",1) ;
    Result_Store(r + i++,&c_ca,"c_ca",1) ;
    Result_Store(r + i++,&c_caoh,"c_caoh",1) ;
    Result_Store(r + i++,&c_h2sio4,"c_h2sio4",1) ;
    Result_Store(r + i++,&c_h3sio4,"c_h3sio4",1) ;
    Result_Store(r + i++,&c_h4sio4,"c_h4sio4",1) ;
    Result_Store(r + i++,&c_cah2sio4,"c_cah2sio4",1) ;
    Result_Store(r + i++,&c_cah3sio4,"c_cah3sio4",1) ;
    Result_Store(r + i++,&c_casaq,"c_casaq",1) ;
    Result_Store(r + i++,&c_cahs,"c_cahs",1) ;
    Result_Store(r + i++,&c_k,"c_k",1) ;
    Result_Store(r + i++,&c_cl,"c_cl",1) ;
    Result_Store(r + i++,&c_oh,"c_oh",1) ;
    Result_Store(r + i++,&c_h,"c_h",1) ;
    Result_Store(r + i++,&n_ch,"n_ch",1) ;
    Result_Store(r + i++,&n_cas,"n_cas",1) ;
    Result_Store(r + i++,&n_si_sol,"n_csh",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&psi,"potentiel_electrique",1) ;
    Result_Store(r + i++,&c_q,"charge",1) ;
    Result_Store(r + i++,&zn_ca_sol,"zn_ca_sol",1) ;
    Result_Store(r + i++,&zc_h2s,"zc_h2s",1) ;
    Result_Store(r + i++,&pk_ch,"pk_ch",1) ;
    Result_Store(r + i++,&zn_si_sol,"zn_si_sol",1) ;
    Result_Store(r + i++,&v_csh,"V_CSH",1) ;
    Result_Store(r + i++,&x_csh,"C/S",1) ;
  }
  
  if(i != nso) arret("ComputeOutputs") ;

  return(nso) ;
}


void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int i ;
  
  /*
    Donnees
  */
  phi0      = GetProperty("porosite") ;
  c_h2s_eq = GetProperty("C_H2S_eq") ;
  /*k_int     = GetProperty("k_int") ;*/
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i = 0 ; i < 2 ; i++) {
    /* molarities */
    double *x = ComputeComponents(el,u,f,0,i) ;

    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hs       = x[I_C_HS] ;
    double c_s        = x[I_C_S] ;
    double c_ca       = x[I_C_Ca] ;
    double c_cahs     = x[I_C_CaHS] ;
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_caoh     = x[I_C_CaOH] ;
    double c_k        = x[I_C_K] ;
    double c_cl       = x[I_C_Cl] ;

    /* porosity */
    double phi        = x[I_Phi] ;
    
    /* Gypsum contents */
    /* double n_cas     = x[I_N_CaS] ; */
    
    /* tortuosite liquide */
    double iff    = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    /*double iff    = 0 ;*/
    /* permeabilite */
    /* double k_g    = (k_int/mu_cas)*pow(phi_c/phi0,3.)*pow(((1-phi0)/(1-phi_c)),2.) ; */
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_H2S        += d_h2s*iff;
    KF_HS         += d_hs*iff;
    KF_S          += d_s*iff;

    KF_Ca         += d_ca*iff ;
    KF_CaHS       += d_cahs*iff;
    KF_CaH3SiO4   += d_cah3sio4*iff;

    KF_H3SiO4     += d_h3sio4*iff ;
    KF_H4SiO4     += d_h4sio4*iff ;
    KF_H2SiO4     += d_h2sio4*iff ;
    KF_CaH2SiO4   += d_cah2sio4*iff ;
    KF_CaSaq      += d_casaq*iff;
    KF_CaOH       += d_caoh*iff ;
    
    KF_K          += d_k*iff;
    KF_Cl         += d_cl*iff;


    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HS       += FsRT*KF_HS*z_hs*c_hs ;
    Kpsi_S        += FsRT*KF_S*z_s*c_s ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaHS     += FsRT*KF_CaHS*z_cahs*c_cahs ;
    Kpsi_CaH3SiO4 += FsRT*KF_CaH3SiO4*z_cah3sio4*c_cah3sio4 ;
    Kpsi_H3SiO4   += FsRT*KF_H3SiO4*z_h3sio4*c_h3sio4 ;
    Kpsi_H2SiO4   += FsRT*KF_H2SiO4*z_h2sio4*c_h2sio4 ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    
    Kpsi_K        += FsRT*KF_K*z_k*c_k ;
    Kpsi_Cl       += FsRT*KF_Cl*z_cl*c_cl ;
    
    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH \
                   + z_hs*Kpsi_HS + z_s*Kpsi_S \
                   + z_ca*Kpsi_Ca + z_caoh*Kpsi_CaOH \
                   + z_h2sio4*Kpsi_H2SiO4 + z_h3sio4*Kpsi_H3SiO4 \
                   + z_cah3sio4*Kpsi_CaH3SiO4 + z_cahs*Kpsi_CaHS \
                   + z_k*Kpsi_K + z_cl*Kpsi_Cl ;
    /*KD_CaS       += n_cas*k_g;*/
  }
  
  /* Averaging */
  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;
}




void ComputeFluxes(Element_t *el,double **u)
/* Les flux (f) */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *grd = dComponents ;

  /* Gradients (electric potential included) */
  {
    double *x1 = ComputeComponents(el,u,f,0.,1) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] = x1[i] ;
  }
  {
    double *x0 = ComputeComponents(el,u,f,0.,0) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] -= x0[i] ;
  }
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] /= dx ;
  }
  
  /* Fluxes */
  {
    double *w = Fluxes(el,grd) ;

    W_S     = w[I_W_S] ;
    W_Ca    = w[I_W_Ca] ;
    W_Si    = w[I_W_Si] ;
    W_q     = w[I_W_q] ;
    W_K     = w[I_W_K] ;
    W_Cl    = w[I_W_Cl] ;
  }
    
}


double* Fluxes(Element_t *el,double *grd)
{
  double *va = Element_GetExplicitTerm(el) ;
  double *w  = ComponentFluxes ;

  /* Gradients */
  double grd_h        = grd[I_C_H] ;
  double grd_oh       = grd[I_C_OH] ;
  
  double grd_h2s      = grd[I_C_H2S] ;
  double grd_hs       = grd[I_C_HS] ;
  double grd_s        = grd[I_C_S] ;
  
  double grd_ca       = grd[I_C_Ca] ;
  double grd_caoh     = grd[I_C_CaOH] ;
  double grd_cahs     = grd[I_C_CaHS] ;
  double grd_casaq    = grd[I_C_CaSaq] ;
  double grd_cah3sio4 = grd[I_C_CaH3SiO4] ;
  double grd_cah2sio4 = grd[I_C_CaH2SiO4] ;
  
  double grd_h2sio4   = grd[I_C_H2SiO4] ;
  double grd_h3sio4   = grd[I_C_H3SiO4] ;
  double grd_h4sio4   = grd[I_C_H4SiO4] ;
  
  double grd_k        = grd[I_C_K] ;
  
  double grd_cl       = grd[I_C_Cl] ;
  
  double grd_psi      = grd[I_PSI] ;
    
    
  /* Flux */
  double w_h2s      = - KF_H2S*grd_h2s ;
  double w_hs       = - KF_HS*grd_hs             - Kpsi_HS*grd_psi ;
  double w_s        = - KF_S*grd_s               - Kpsi_S*grd_psi  ;
    
  double w_cahs     = - KF_CaHS*grd_cahs          - Kpsi_CaHS*grd_psi ;
  double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
  double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi ;
  double w_casaq    = - KF_CaSaq*grd_casaq ;
  double w_cah2sio4 = - KF_CaH2SiO4*grd_cah2sio4 ;
  double w_cah3sio4 = - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
    
  double w_h3sio4   = - KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi ;
  double w_h2sio4   = - KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi ;
  double w_h4sio4   = - KF_H4SiO4*grd_h4sio4 ;
    
  double w_k        = - KF_K*grd_k                - Kpsi_K*grd_psi ;
  double w_cl       = - KF_Cl*grd_cl              - Kpsi_Cl*grd_psi ; 
    
  double w_q        = - z_h*KF_H*grd_h \
                      - z_oh*KF_OH*grd_oh	\
                      - z_hs*KF_HS*grd_hs \
                      - z_s*KF_S*grd_s	\
                      - z_ca*KF_Ca*grd_ca	\
                      - z_cahs*KF_CaHS*grd_cahs \
                      - z_h3sio4*KF_H3SiO4*grd_h3sio4	\
                      - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                      - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                      - z_caoh*KF_CaOH*grd_caoh \
                      - z_k*KF_K*grd_k \
                      - z_cl*KF_Cl*grd_cl \
                      - Kpsi_q*grd_psi ;

  w[I_W_S  ]  = w_h2s + w_hs + w_s + w_cahs + w_casaq ;
  w[I_W_Ca ]  = w_ca + w_cahs + w_cah3sio4 + w_casaq + w_caoh + w_cah2sio4 ; 
  w[I_W_Si ]  = w_h3sio4 + w_h4sio4 + w_cah3sio4 + w_cah2sio4 + w_h2sio4 ;
  w[I_W_q  ]  = w_q ;
  w[I_W_K  ]  = w_k ;
  w[I_W_Cl ]  = w_cl ;
   
  return(w) ;
}




int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double x1    = Element_GetNodeCoordinate(el,1)[0] ;
  double x0    = Element_GetNodeCoordinate(el,0)[0] ;
  double dij   = fabs(x1 - x0) ;
  double dtdij = dt/dij ;
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
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    double dxi[NEQ] ;
    int k ;
    
    dxi[I_C_H2S    ] = 1.e-8*x[I_C_H2S] ;
    dxi[I_ZN_Si_sol] = 1.e-8 ;
    dxi[I_ZN_Ca_sol] = 1.e-10 ;
    dxi[I_C_K      ] = 1.e-6 ;
    dxi[I_C_Cl     ] = 1.e-6 ;
    dxi[I_PSI      ] = 1. ;

    /*
    dxi[I_C_H2S    ] =  ObVal_GetValue(obval + I_C_H2S) ;
    dxi[I_ZN_Si_sol] =  ObVal_GetValue(obval + I_ZN_Si_sol) ;
    dxi[I_ZN_Ca_sol] =  ObVal_GetValue(obval + I_ZN_Ca_sol) ;
    dxi[I_C_K      ] =  ObVal_GetValue(obval + I_C_K) ;
    dxi[I_C_Cl     ] =  ObVal_GetValue(obval + I_C_Cl) ;
    dxi[I_PSI      ] =  ObVal_GetValue(obval + I_PSI) ;
    * 
    for(k = 0 ; k < NEQ ; k++) {
      dxi[k] =  1.e-2*ObVal_GetValue(obval + k) ;
    }
    */
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      double *dw    = Fluxes(el,dx) ;
    
      cii[E_S*NEQ    + k] = dx[I_N_S] ;
      cii[E_Ca*NEQ   + k] = dx[I_N_Ca] ;
      cii[E_Si*NEQ   + k] = dx[I_N_Si] ;
      cii[E_K*NEQ    + k] = dx[I_N_K] ;
      cii[E_Cl*NEQ   + k] = dx[I_N_Cl] ;
      
      cij[E_S*NEQ    + k] = - dtdij*dw[I_W_S] ;
      cij[E_Ca*NEQ   + k] = - dtdij*dw[I_W_Ca] ;
      cij[E_Si*NEQ   + k] = - dtdij*dw[I_W_Si] ;
      cij[E_K*NEQ    + k] = - dtdij*dw[I_W_K] ;
      cij[E_Cl*NEQ   + k] = - dtdij*dw[I_W_Cl] ;
      cij[E_q*NEQ    + k] = - dtdij*dw[I_W_q] ;
    }
  }

  return(dec) ;
}



double concentration_oh(double c_h2s,elem_t *el,double zn_ca_sol,double c_k,double c_cl,double zn_si_sol)
/* Solve electroneutrality: SUM(z_i c_i) = 0
   for c_oh as root of ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_h2s, zn_si_sol and zn_ca_sol are considered as constant */
  double zc_h2s    = c_h2s/c_h2s_eq ;
  /* Ion acitivity products are constant as well */
  double Q_CaS     = IonActivityProductOfCaS(zc_h2s,zn_ca_sol) ;
  double s_ch      = SaturationDegreeOfCH(zc_h2s,zn_ca_sol) ;
  double Q_SH      = IonActivityProductOfSH(s_ch,zn_si_sol) ;
  
  double c_h4sio4   = Q_SH ;
  /*
  rappel des expressions c_i = A_i*(c_h)**n   : n
     c_h        = K_h2o/c_oh                  : +1
     c_hs       = K_hs*c_h2s/c_h              : -1
     c_s        = K_s*c_hs/c_h                : -2
     c_ca       = Q_CaS/c_s                   : +2
     c_cahs     = K_cahs*c_ca*c_hs            : +1
     c_h4sio4   = Q_SH                        :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)     : -1
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4    : +1
     c_casaq    = K_casaq*c_ca*c_s            :  0      
     c_h2sio4   = K_h2sio4*c_h3sio4*c_oh      : -2       
     c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca    :  0      
     c_caoh     = K_caoh*c_ca*c_oh            : +1       
  */
  double A_hs       = K_hs*c_h2s ;
  double A_s        = K_s*A_hs ;
  double A_ca       = Q_CaS/A_s ;
  double A_cahs     = K_cahs*A_ca*A_hs ;
  double A_h3sio4   = c_h4sio4/K_h4sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;
  double A_h2sio4   = K_h2sio4*A_h3sio4*K_h2o ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahs*A_cahs + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_k*c_k + z_cl*c_cl;
  double d = z_oh*K_h2o + z_hs*A_hs + z_h3sio4*A_h3sio4 ;
  double e = z_s*A_s + z_h2sio4*A_h2sio4 ;

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
  
  Math_PolishPolynomialEquationRoot(y,4,&x,tol*x,20) ;
  
  return(x) ;
}



double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *v0 = Element_GetConstantTerm(el) ;
  double *x = Components ;
  
  /* Primary Variables */
  x[I_C_H2S     ] = C_H2S(n) ;
  x[I_ZN_Ca_sol ] = ZN_Ca_sol(n) ;
  x[I_ZN_Si_sol ] = ZN_Si_sol(n) ;
  x[I_C_K       ] = C_K(n) ;
  x[I_C_Cl      ] = C_Cl(n) ;  
  x[I_PSI       ] = PSI(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn  ]  = N_CHn(n) ;
  x[I_N_CaSn ]  = N_CaSn(n) ;
  x[I_V_Mat0 ]  = V_Mat0(n) ;
  
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  dx[I_C_H2S     ] = x[I_C_H2S   ] ;
  dx[I_ZN_Ca_sol ] = x[I_ZN_Ca_sol ] ;
  dx[I_ZN_Si_sol ] = x[I_ZN_Si_sol ] ;
  dx[I_C_K       ] = x[I_C_K     ] ;
  dx[I_C_Cl      ] = x[I_C_Cl    ] ;
  dx[I_PSI       ] = x[I_PSI     ] ;

  /* Needed variables to compute secondary components */
  dx[I_N_CHn   ] = x[I_N_CHn ] ;
  dx[I_N_CaSn  ] = x[I_N_CaSn] ;
  dx[I_V_Mat0  ] = x[I_V_Mat0  ] ;
  
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



void  ComputeSecondaryComponents(Element_t *el,double dt,double *x)
{
  double c_h2s      = x[I_C_H2S] ;
  double zn_si_sol  = x[I_ZN_Si_sol] ;
  double zn_ca_sol  = x[I_ZN_Ca_sol] ;
  double c_k        = x[I_C_K] ;
  double c_cl       = x[I_C_Cl] ;      
    
  /* Liquid components */
  double zc_h2s     = c_h2s/c_h2s_eq ;

  double Q_CaS      = IonActivityProductOfCaS(zc_h2s,zn_ca_sol) ;
  double s_ch       = SaturationDegreeOfCH(zc_h2s,zn_ca_sol) ;
  double Q_SH       = IonActivityProductOfSH(s_ch,zn_si_sol) ;

  double c_h4sio4   = Q_SH ;

  double c_oh       = ConcentrationOfOHInLiquid(c_h2s,zn_ca_sol,c_k,c_cl,zn_si_sol) ;
  double c_h        = K_h2o/c_oh ;
  double c_hs       = K_hs*c_h2s/c_h ;
  double c_s        = K_s*c_hs/c_h ;
  double c_ca       = Q_CaS/c_s ;
  double c_cahs     = K_cahs*c_ca*c_hs ;
  double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
  double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
  double c_casaq    = K_casaq*c_ca*c_s ;
  double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
  double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
  double c_caoh     = K_caoh*c_ca*c_oh ;
  
  double c_q        = z_h*c_h + z_oh*c_oh \
                    + z_hs*c_hs + z_s*c_s \
                    + z_ca*c_ca + z_caoh*c_caoh \
                    + z_h3sio4*c_h3sio4 + z_h2sio4*c_h2sio4 \
                    + z_cah3sio4*c_cah3sio4 + z_cahs*c_cahs \
                    + z_k*c_k + z_cl*c_cl ;
  
  /* Solid contents */
  /* ... as components: CH, CaS, CSH */
  double n_chn      = x[I_N_CHn] ;
  double n_casn     = x[I_N_CaSn] ;
  double n_ch_ci    = n_chn*pow(zc_h2s,-dt/t_ch) ;
  double n_cas_ci   = n_casn*pow(zc_h2s,dt/t_cas) ;
  double n_ch_cas   = CalciumContentInCHAndCaS(zn_ca_sol) ;
  double n_ch       = (zc_h2s <= 1) ? n_ch_cas - n_cas_ci : n_ch_ci ;
  double n_cas      = (zc_h2s >  1) ? n_ch_cas - n_ch_ci : n_cas_ci ;
  /* ... as elements: S, Ca, Si */
  double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double n_si_sol   = SiliciumContentInCSH(zn_si_sol) ;
  double n_ca_sol   = n_ch_cas + x_csh*n_si_sol ;
  double n_s_sol    = n_cas ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(s_ch) ;
  double v_mat      = V_CH*n_ch + V_CaS*n_cas + v_csh*n_si_sol ;
  
  /* Porosity */
  double v_mat0     = x[I_V_Mat0] ;
  double phi1       = phi0 + v_mat0 - v_mat ;
  double phi2       = (phi1 < 0 ) ? 0.00001 : phi1 ;
  double phi        = phi2 ;
  
  /* Liquid contents */
  /* ... as elements: S, Ca, Si, K, Cl */
  double c_s_liq  = c_h2s + c_hs + c_s + c_cahs + c_casaq ;
  double c_ca_liq = c_ca + c_cahs + c_cah3sio4 + c_cah2sio4 + c_casaq + c_caoh ;
  double c_si_liq = c_h3sio4 + c_h4sio4 + c_cah3sio4 + c_h2sio4 + c_cah2sio4 ;
  double c_k_liq  = c_k ;
  double c_cl_liq = c_cl ;
  double n_s_liq  = phi*c_s_liq ;
  double n_ca_liq = phi*c_ca_liq ;
  double n_si_liq = phi*c_si_liq ;
  double n_k_liq  = phi*c_k_liq ;
  double n_cl_liq = phi*c_cl_liq ;
  
  
  /* Back up */
  
  /* Liquid components */
  x[I_C_OH      ] = c_oh ;
  x[I_C_H       ] = c_h ;
  
  x[I_C_HS      ] = c_hs ;
  x[I_C_S       ] = c_s ;
  
  x[I_C_Ca      ] = c_ca ;
  
  x[I_C_H4SiO4  ] = c_h4sio4 ;
  x[I_C_H3SiO4  ] = c_h3sio4 ;
  x[I_C_H2SiO4  ] = c_h2sio4 ;
  
  x[I_C_CaH3SiO4] = c_cah3sio4 ;
  x[I_C_CaH2SiO4] = c_cah2sio4 ;
  
  x[I_C_CaHS    ] = c_cahs ;
  x[I_C_CaSaq   ] = c_casaq ;
  
  x[I_C_CaOH    ] = c_caoh ;
  
  x[I_S_CH      ] = s_ch ;


  /* Solid components */
  x[I_N_CH      ] = n_ch ;
  x[I_N_CaS     ] = n_cas ;
  x[I_V_Mat     ] = v_mat ;
  x[I_N_Si_sol  ] = n_si_sol ;
  
  /* Porosity */
  x[I_Phi       ] = phi ;
  
  /* Element contents */
  x[I_N_S       ] = n_s_liq  + n_s_sol ;
  x[I_N_Ca      ] = n_ca_liq + n_ca_sol ;
  x[I_N_Si      ] = n_si_liq + n_si_sol ;
  x[I_N_K       ] = n_k_liq  ;
  x[I_N_Cl      ] = n_cl_liq  ;
  /* Charge density */
  x[I_N_Q       ] = c_q ;

}
