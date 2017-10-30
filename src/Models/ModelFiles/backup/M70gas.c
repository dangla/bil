/*

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

/* Choose the numerical method */
#include "FVM.h"

#define TITLE "Hydrogen Sulfide (H2S) attack of concrete (May. 2012) + K+Cl" 
#define AUTHORS "Yuan"

#include "PredefinedMethods.h"

/* Macros */ 
#define NEQ     (6)
#define NVE     (31)/* g*/
#define NVE_TR  (31)/* g*/
#define NVI     (29)
#define NV0     (8)

#define E_S     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)
#define E_K     (4)
#define E_Cl    (5)

#define I_C_H2SF  (0)
#define I_PSI      (1)
#define I_ZN_Ca_S  (2)
#define I_ZN_Si_S  (3)
#define I_C_K      (4)
#define I_C_Cl     (5)

#define RHO     1
#define LOG_RHO 2
#define Ln10    2.302585093
#define U_H2SF   LOG_RHO
#define EXPLICIT  1
#define IMPLICIT  2
#define U_PHI     IMPLICIT

#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWN_n(n,i)   (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])

/*
#define UNKNOWN(n,i)     (Element_GetValueOfNodalUnknown(el,n,i))
#define UNKNOWN_n(n,i)   (Element_GetValueOfPreviousNodalUnknown(el,n,i))
*/


#if (U_H2SF == LOG_RHO)
  #define C_H2SF(n)   (exp(Ln10*UNKNOWN(n,I_C_H2SF)))
  #define C_H2SFn(n)  (exp(Ln10*UNKNOWN_n(n,I_C_H2SF)))
#else
  #define C_H2SF(n)   (UNKNOWN(n,I_C_H2SF))
  #define C_H2SFn(n)  (UNKNOWN_n(n,I_C_H2SF))
#endif
#define ZN_Ca_S(n)   (UNKNOWN(n,I_ZN_Ca_S))
#define ZN_Si_S(n)   (UNKNOWN(n,I_ZN_Si_S))
#define PSI(n)       (UNKNOWN(n,I_PSI))
#define C_K(n)       (UNKNOWN(n,I_C_K))
#define C_Cl(n)      (UNKNOWN(n,I_C_Cl))

#define ZN_Ca_Sn(n)  (UNKNOWN_n(n,I_ZN_Ca_S))
#define ZN_Si_Sn(n)  (UNKNOWN_n(n,I_ZN_Si_S))
#define PSIn(n)      (UNKNOWN_n(n,I_PSI))
#define C_Kn(n)      (UNKNOWN_n(n,I_C_K))
#define C_Cln(n)     (UNKNOWN_n(n,I_C_Cl))

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
#define A_H2S      (f[18])
#define N_CH(n)    (f[(19+n)])
#define N_CSF(n)  (f[(21+n)])
#define N_Si_S(n)  (f[(23+n)])
#define ZQ_CH(n)   (f[(25+n)])
#define PHI(n)     (f[(27+n)])

#define N_Sn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_Cln(n)   (f_n[(10+n)])
#define N_CHn(n)   (f_n[(19+n)])
#define N_CSFn(n) (f_n[(21+n)])
#define N_Si_Sn(n) (f_n[(23+n)])
#define ZQ_CHn(n)  (f_n[(25+n)])
#define PHIn(n)    (f_n[(27+n)])

#define KF_OH       (va[(0)])
#define KF_H        (va[(1)])
#define KF_CO2      (va[(2)])
#define KF_H2SF    (va[(3)])
#define KF_HSF     (va[(4)])
#define KF_SF      (va[(5)])
#define KF_Ca       (va[(6)])
#define KF_CaHSF   (va[(7)])
#define KF_CaH3SiO4 (va[(8)])
#define KF_H3SiO4   (va[(9)])
#define KF_H4SiO4   (va[(10)])
#define KF_H2SiO4   (va[(11)])
#define KF_CaH2SiO4 (va[(12)])
#define KF_CaSFaq  (va[(13)])
#define KF_CaOH     (va[(14)])
#define KF_K        (va[(15)])
#define KF_Cl       (va[(16)])

#define Kpsi_OH       (va[(17)])
#define Kpsi_H        (va[(18)])
#define Kpsi_HSF     (va[(19)])
#define Kpsi_SF      (va[(20)])
#define Kpsi_Ca       (va[(21)])
#define Kpsi_CaHSF   (va[(22)])
#define Kpsi_CaH3SiO4 (va[(23)])
#define Kpsi_H3SiO4   (va[(24)])
#define Kpsi_q        (va[(25)])
#define Kpsi_H2SiO4   (va[(26)])
#define Kpsi_CaOH     (va[(27)])
#define Kpsi_K        (va[(28)])
#define Kpsi_Cl       (va[(29)])
#define KD_CSF       (va[(30)])/* g*/

#define N_CH0(n)      (v0[(0+n)])
#define N_CSF0(n)    (v0[(2+n)])
#define N_Si_S0(n)    (v0[(4+n)])
#define ZQ_CH0(n)     (v0[(6+n)])

/*
  Solution aqueuse
*/

/* les valences */
#define z_ca       (2.)
#define z_h        (1.)
#define z_oh       (-1.)
#define z_hsf     (-1.)
#define z_sf      (-2.)
#define z_h3sio4   (-1.)
#define z_cahsf   (1.)
#define z_cah3sio4 (1.)
#define z_h2sio4   (-2.)
#define z_caoh     (1.)
#define z_k        (1.)
#define z_cl       (-1.)

/* volumes molaires partiels des ions (dm3/mole) from [Millero F J,Partial molar volum of ions in seawater]*/
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_h2o      (18.e-3)
#define v_co2      (32.81e-3)
#define v_h2sf    (34.9e-3)
#define v_hsf     (27.3e-3)   /* d'apres Lothenbach */
#define v_sf      (-0.2e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)
#define v_casfaq   (26.20e-3)    /* a modifier */
#define v_caoh     (26.20e-3)    /* a modifier */
#define v_k        (43.93e-3)    /* d'apres Antoine */
#define v_cl       (43.93e-3)    /*a modifier*/
#define v_koh      (27.44e-3)    /* d'apres Antoine */

/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s  */
#define d_oh       (1.22e-7)    /* 1.22e-7 (radius = 1.75e-10 m) */
#define d_h        (9.310e-7)    /* 4.76e-8 (radius = 4.5e-10 m) */
#define d_h2sf     (1.36e-9)   /* (radius = 1.5e-10 m) */
#define d_hsf      (1.73e-7)   /* (radius = 1.91e-10 m) */
#define d_sf       (1.13e-7)   /*? (radius = 1.89e-10 m) */
#define d_ca       (7.92e-8)
#define d_cahsf   (1.07e-7)    /* (radius = 2e-10 m) */
#define d_sioh4    (xxx)
#define d_h4sio4   (1.07e-7)     /* */
#define d_h3sio4   (1.07e-7)   /* (radius = 2e-10 m) */
#define d_h2sio4   (1.07e-7)    /*(radius = 2e-10 m) */
#define d_cah2sio4 (1.07e-7)    /* a modifier */
#define d_cah3sio4 (1.07e-7)     /*(radius = 2e-10 m) */
#define d_casfaq  (1.43e-7)    /* (radius = 1.5e-10 m) */
#define d_caoh     (1.07e-7)    /* (radius = 2e-10 m) */
#define d_k        (1.43e-7)   /* (radius = 1.5e-10 m) */
#define d_cl       (2.032e-7)   /* chemistry handbook */
#define d_koh      (1.43e-7)   /* (radius = 1.5e-10 m) */

/* viscosite (Pa.s) of csf*/
#define mu_csf       	(30)

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_h2o      (1.e-14)          /* autoprotolyse de l'eau */
#define K_henry    (2.5)           /* H2S KH=0.1; (KH*RT);cste de Henry du H2S / RT CO2 1.238*/

#define K_hsf      (8.9e-8)          /* H2S   = HS[-] + H[+] */
#define K_sf      (1.2e-13)        /* HS[-] = S[2-] + H[+] */

#define K_h2sio4   (4.68)            /* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define K_h4sio4   (6.45e9)          /* H3SiO4[-] + H[+] = H4SiO4 */
#define K_h3sio4   (1.55e-10)        /* H4SiO4    = H3SiO4[-] + H[+] */

#define K_cahsf    (1.276e+1)        /* ???Ca[2+] + HCO3[-]    = CaHCO3[+] */
#define K_casfaq   (3.5e+3)          /* Ca[2+] + S[2-]    = CaS[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-] = CaOH[+] */

#define KF_A_H2S      (2e-9)            /* mol/dm2 s */

/*
  Solides
  CH  = Portlandite
  CSF  = Calcuim Sulfide
  CSH = Hydrated Calcium Silicates
  SH  = Amorphous Silica
*/

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* 3.456e-9 6.456e-6 CH  = Ca[2+] + 2OH[-] */
#define K_SH       (1.93642e-3)      /* SHt = H4SiO4 + (t-2)H2O */
/* volumes molaires solides (dm3/mole) */
#define V_CH       (0)      /* (33.e-3) */

/* CSF */
/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CSF     (7.9e-7)         /* CSF = Ca[2+] + sf[2-] + 2H2O */
/* volumes molaires solides (dm3/mole) */
#define V_CSF     (0)      /* (28.e-3) */
/* Coefficient of compressibility (MPa) */
#define CC_CSF    (5.e+06)      /* (??) */
#define CC_Con     (5.e+9)      /* (??) */
#define CS_Con     (3.e+6)      /* (??) */
 

/* Curves for C-S-H */
/* C/S ratio */
#define X_CSH(q)    (Curve_ComputeValue(Element_GetCurve(el),q))
#define DX_CSH(q)   (Curve_ComputeDerivative(Element_GetCurve(el),q))
/* Molar Volume */
#define V_CSH(q)    (Curve_ComputeValue(Element_GetCurve(el) + 2,q))
#define DV_CSH(q)   (Curve_ComputeDerivative(Element_GetCurve(el) + 2,q))
/* Amorphous Silica Q_SH/K_SH */
#define Q_SH(q)     (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define DQ_SH(q)    (Curve_ComputeDerivative(Element_GetCurve(el) + 3,q))


#define C_H2SF_eq (K_h2o*K_h2o*K_CSF/(K_hsf*K_sf*K_CH))

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */
/*#define BOUN      (15)   /* F/RT (1/V) */

/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,elem_t*,double,double,double,double) ;
static double productionrate(double,elem_t *el, double, double);
static double absorptionrate(double,elem_t *el);
static double poly4(double,double,double,double,double) ;
static void   transfert(Element_t*,double**,double*) ;
static void   flux(Element_t*,double**) ;

static double poly41(double,double,double,double,double) ;
static double quartic_equation(double,double,double,double,double,double *,double *);
static double cubic_equation(double,double,double,double);
static void   quadratic_equation(double,double,double,double *);

static double* ComputeAqueousVariables(Element_t*,double**,int) ;
static void    ComputeAqueousParameters(Element_t*,double*) ;
static double* ComputeAqueousDerivatives(Element_t*,double*,double,int) ;
static double* ComputeAqueousDerivatives1(Element_t*,double*,double,int) ;

#define MIN(a,b)   ((a < b) ? a : b)
#define MAX(a,b)   ((a > b) ? a : b)
#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)
#define DNEGEXP(y) ((y < 0.) ? exp(y) : 0.)

/* Ion Activity Products */
#define IAP_CSF(zc_h2sf,zn_ca_s)       (K_CSF*NEGEXP(zn_ca_s)*MIN(zc_h2sf,1.))
#define IAP_CH(zc_h2sf,zn_ca_s)         (K_CH*NEGEXP(zn_ca_s)/MAX(zc_h2sf,1.))
#define IAP_SH(zq_ch,zn_si_s)            (K_SH*NEGEXP(zn_si_s)*Q_SH(zq_ch))
#define DIAP_SHSDQ_CH(zq_ch,zn_si_s)     (K_SH*NEGEXP(zn_si_s)*DQ_SH(zq_ch))
#define DIAP_SHSDZN_Si_S(zq_ch,zn_si_s)  (K_SH*DNEGEXP(zn_si_s)*Q_SH(zq_ch))


#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_S,"carbone") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_Si,"silicium") ;
  Model_CopyNameOfEquation(model,E_q,"charge") ;
  Model_CopyNameOfEquation(model,E_K,"potassium") ;
  Model_CopyNameOfEquation(model,E_Cl,"chloride") ;

#if (U_H2SF == LOG_RHO)
  Model_CopyNameOfUnknown(model,I_C_H2SF,"logc_h2sf") ;
#else
  Model_CopyNameOfUnknown(model,I_C_H2SF,"c_h2sf") ;
#endif
  Model_CopyNameOfUnknown(model,I_ZN_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,I_PSI,"psi") ;
  Model_CopyNameOfUnknown(model,I_ZN_Si_S,"z_si") ;
  Model_CopyNameOfUnknown(model,I_C_K,"c_k") ;
  Model_CopyNameOfUnknown(model,I_C_Cl,"c_cl") ;
  
  return(0) ;
}


/* Parametres */
static double phi0,c_h2sf_eq,t_ch,t_csf,p_s2h;
static double n_ca_ref,n_si_ref ;

#define NbOfAqueousVariables    (26)
static double var[NbOfAqueousVariables],dvar[NbOfAqueousVariables] ;

#define I_C_OH         (6)
#define I_C_H          (7)
#define I_C_HSF       (8)
#define I_C_SF        (9)
#define I_C_Ca         (10)
#define I_C_CaOH       (11)
#define I_C_CaHSF     (12)
#define I_C_CaSFaq    (13)
#define I_C_H2SiO4     (14)
#define I_C_H3SiO4     (15)
#define I_C_H4SiO4     (16)
#define I_C_CaH2SiO4   (17)
#define I_C_CaH3SiO4   (18)
#define I_ZQ_CH        (19)
#define I_N_Q          (20)
#define I_C_S_L        (21)
#define I_C_Ca_L       (22)
#define I_C_Si_L       (23)
#define I_C_K_L        (24)
#define I_C_Cl_L       (25)


int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"N_Si") == 0)     return (2) ;
  else if(strcmp(s,"C_H2SF_eq") == 0) return (3) ;
  else if(strcmp(s,"T_CH") == 0)     return (4) ;
  else if(strcmp(s,"T_CSF") == 0)     return (5) ;
  else if(strcmp(s,"P_H2S") == 0)    return (6) ;
  else return(-1) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 7 ;

  {
    /* initialisation automatique */
    double t_ch        = 600. ;
    double t_csf      = 600. ;
    double p_h2s      = 50. ;
    double n_ca_ref    = 1. ;
    double n_si_ref    = 1. ;

    Material_GetProperty(mat)[pm("N_CH")] = n_ca_ref ;
    Material_GetProperty(mat)[pm("N_Si")] = n_si_ref ;
    Material_GetProperty(mat)[pm("T_CH")] = t_ch ;
    Material_GetProperty(mat)[pm("T_CSF")] = t_csf ;
    Material_GetProperty(mat)[pm("P_H2S")] = p_h2s ;

    Material_ScanProperties(mat,datafile,pm) ;

    t_ch      = Material_GetProperty(mat)[pm("T_CH")] ;
    t_csf    = Material_GetProperty(mat)[pm("T_CSF")] ;
    p_h2s    = Material_GetProperty(mat)[pm("P_H2S")] ;

    if(t_csf  == 0.) Material_GetProperty(mat)[pm("T_CSF")]  = t_ch ;

    Material_GetProperty(mat)[pm("C_H2SF_eq")] = C_H2SF_eq ;
  }
  
  return(n_donnees) ;
}



int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;
  printf("Le systeme est forme de 6 equations:\n") ;
#if (U_H2SF == LOG_RHO)
  printf("\t- la conservation de la masse de C      (logc_h2sf)\n") ;
#else
  printf("\t- la conservation de la masse de C      (c_h2sf)\n") ;
#endif
  printf("\t- la conservation de la charge          (psi)\n") ;
  printf("\t- la conservation de la masse de Ca     (zn_ca_s)\n") ;
  printf("\t- la conservation de la masse de Si     (zn_si_s)\n") ;
  printf("\t- la conservation de la masse de K      (c_k)\n") ;
  printf("\t- la conservation de la masse de Cl     (c_cl)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n") ;

  printf("Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1       # Contenu en CH (moles/L)\n") ;
  fprintf(ficd,"N_Si  = 2.4       # contenu en Si solide (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4       # contenu en K (moles/L)\n") ;
  fprintf(ficd,"N_Cl   = 0.4       # contenu en Cl (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5      # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"T_CSF  = 1.e5      # Cinetique de dissolution de CSF (s)\n") ;
  fprintf(ficd,"courbes = solid   # Nom du fichier: q_ch X Y V f_S\n") ;

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
  FVM_t* fvm = FVM_GetInstance(el) ;
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
  double *u[MAX_NOEUDS] ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Donnees
  */
  phi0      = Element_GetProperty(el)[pm("porosite")] ;
  n_ca_ref  = Element_GetProperty(el)[pm("N_CH")] ;
  n_si_ref  = Element_GetProperty(el)[pm("N_Si")] ;
  c_h2sf_eq  = Element_GetProperty(el)[pm("C_H2SF_eq")] ;
  /*k_int       = Element_GetProperty(el)[pm("k_int")] ;*/
  /* Default initialization */
  for(i = 0 ; i < nn ; i++) {
    double c_h2sf      = (C_H2SF(i) > 0.) ? C_H2SF(i) : c_h2sf_eq ;

#if (U_H2SF == LOG_RHO)
    UNKNOWN(i,I_C_H2SF) = log(c_h2sf)/Ln10 ;
#else
    C_H2SF(i)   = c_h2sf ;
#endif
  }
  
  
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
    double *x = ComputeAqueousVariables(el,u,i) ;
    /*printf("c_cl         = %e\n",x[I_C_Cl]) ;*/
    double xx = Element_GetNodeCoordinate(el,i)[0] ;
    
    double c_h2sf    = x[I_C_H2SF] ;
    double zn_si_s    = x[I_ZN_Si_S] ;
    double zn_ca_s    = x[I_ZN_Ca_S] ;

    
    
    double zc_h2sf   = c_h2sf/c_h2sf_eq ;
    double zq_ch      = x[I_ZQ_CH] ;
    
    /* solid contents : CH, CSF, CSH */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_csf_eq  = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_ch       = (zc_h2sf <= 1) ? n_ch_eq  : 0 ;
    double n_csf     = (zc_h2sf >  1) ? n_csf_eq  : 0 ;
    double x_csh      = X_CSH(zq_ch) ;
    double n_si_s     = n_si_ref*MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_csf + x_csh*n_si_s ;
    double n_s_s      = n_csf ;


    /* porosity */
    double phi = phi0 ;

   

    /* liquid contents */
    double n_s_l      = phi*x[I_C_S_L] ;
    double n_ca_l     = phi*x[I_C_Ca_L] ;
    double n_si_l     = phi*x[I_C_Si_L] ;
    double n_k_l      = phi*x[I_C_K_L] ;
    double n_cl_l     = phi*x[I_C_Cl_L] ;
    
    N_S(i)  = n_s_l  + n_s_s ;
    N_Ca(i) = n_ca_l + n_ca_s ;
    N_Si(i) = n_si_l + n_si_s ;
    N_K(i)  = n_k_l ;
    N_Cl(i) = n_cl_l ;
    /* charge density */
    N_q(i)  = x[I_N_Q] ;

    
    N_CH(i)    = n_ch ;
    N_CSF(i)  = n_csf ;
    N_Si_S(i)  = n_si_s ;
    ZQ_CH(i)   = zq_ch ;
    PHI(i)     = phi ;

    N_CH0(i)   = n_ch ;
    N_CSF0(i) = n_csf ;
    N_Si_S0(i) = n_si_s ;
    ZQ_CH0(i)  = zq_ch ;
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  transfert(el,u,f) ;

  /* Flux */
  flux(el,u) ;
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Thermes explicites (va)  */
{
  double *f = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetPreviousNodalUnknown(el,i) ;
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Coefficients de transfert
  */
  transfert(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Donnees
  */
  phi0      = Element_GetProperty(el)[pm("porosite")] ;
  n_ca_ref  = Element_GetProperty(el)[pm("N_CH")] ;
  n_si_ref  = Element_GetProperty(el)[pm("N_Si")] ;
  c_h2sf_eq  = Element_GetProperty(el)[pm("C_H2SF_eq")] ;
  t_ch      = Element_GetProperty(el)[pm("T_CH")] ;
  t_csf    = Element_GetProperty(el)[pm("T_CSF")] ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
    double *x = ComputeAqueousVariables(el,u,i) ;
    double xx = Element_GetNodeCoordinate(el,i)[0] ;
   
    double c_h2sf    = x[I_C_H2SF] ;
    double zn_si_s    = x[I_ZN_Si_S] ;
    double zn_ca_s    = x[I_ZN_Ca_S] ;
    double zc_h2sf   = c_h2sf/c_h2sf_eq ;
    double zq_ch      = x[I_ZQ_CH] ;

    
    /* solid contents : CH, CSF, CSH */
    double n_chn      = N_CHn(i) ;
    double n_csfn    = N_CSFn(i) ;

    

    /*   kinetics */
    double n_ch_ci    = n_chn*pow(zc_h2sf,-dt/t_ch) ;  /* if zc_h2sf > 1 */
    double n_csf_ci  = n_csfn*pow(zc_h2sf,dt/t_csf) ;   /* if zc_h2sf < 1 */
    /*   equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) - n_csf_ci ;
    double n_csf_eq  = n_ca_ref*MAX(zn_ca_s,0.) - n_ch_ci ;
    double n_ch       = (zc_h2sf <= 1) ? n_ch_eq   : n_ch_ci ;
    double n_csf     = (zc_h2sf >  1) ? n_csf_eq : n_csf_ci ;
    double x_csh      = X_CSH(zq_ch) ;
    double n_si_s     = n_si_ref*MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_csf + x_csh*n_si_s ;
    double n_s_s      = n_csf ;

    /* porosity */
    double n_ch0      = N_CH0(i) ;
    double n_csf0    =  N_CSF0(i) ;
    double n_si_s0    =  N_Si_S0(i) ;
    double zq_ch0     = ZQ_CH0(i) ;
    double v_csh0     = V_CSH(zq_ch0) ;
    double v_csh      = V_CSH(zq_ch) ;
    double nv_agr    = 1 - phi0- V_CH*n_ch0 - v_csh0*n_si_s0 ;/*volume of sand*/
    double phi1       = phi0 + V_CH*(n_ch0 - n_ch) + V_CSF*(n_csf0 - n_csf) + v_csh0*n_si_s0 - v_csh*n_si_s ;
                      
    
#if (U_PHI == IMPLICIT)
    double phi        = phi1;
#else
    double phi        = PHIn(i) ;
#endif
    
    /*if (xx==1) {phi    = phi_CSF;}*/
    /* liquid contents */
    double n_s_l      = phi*x[I_C_S_L] ;
    double n_ca_l     = phi*x[I_C_Ca_L] ;
    double n_si_l     = phi*x[I_C_Si_L] ;
    double n_k_l      = phi*x[I_C_K_L] ;
    double n_cl_l     = phi*x[I_C_Cl_L] ;
    
    /*if (damage == 1 && n_ch ==0 && n_csf  ==0 )   
	{
		phi = 0.99;
		
	}*/ 

    
    /* Molar contents */
    N_S(i)  = n_s_l  + n_s_s ;
    N_Ca(i) = n_ca_l + n_ca_s ;
    N_Si(i) = n_si_l + n_si_s ;
    N_K(i)  = n_k_l ;
    N_Cl(i)  = n_cl_l ;
    /* charge density */
    N_q(i)  = x[I_N_Q] ;
    
    N_CH(i)    = n_ch ;
    N_CSF(i)  = n_csf ;
    N_Si_S(i)  = n_si_s ;
    ZQ_CH(i)   = zq_ch ;
    PHI(i)     = phi ;


    if(c_h2sf < 0. || n_ca_s < 0. || n_si_s < 0.) {
      
      double c_h3sio4 = x[I_C_H3SiO4] ;
      double c_oh     = x[I_C_OH] ;
      printf("x         = %e\n",xx) ;
      printf("c_h2sf   = %e\n",c_h2sf) ;
      printf("n_csf    = %e\n",n_csf) ;
      printf("n_ca_s    = %e\n",n_ca_s) ;
      printf("n_si_s    = %e\n",n_si_s) ;
      printf("zn_si_s   = %e\n",zn_si_s) ;
      printf("zn_ca_s   = %e\n",zn_ca_s) ;
      printf("c_h3sio4  = %e\n",c_h3sio4) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
 
    if(phi < 0.) {
	  
      double c_h3sio4 = x[I_C_H3SiO4] ;
      double c_oh     = x[I_C_OH] ;
      printf("phi = %e\n",phi) ;
      printf("CH = %e\n",n_ch) ;
      printf("CSF = %e\n",n_csf) ;
      printf("Si = %e\n",n_si_s) ;
      printf("c_h2sf   = %e\n",c_h2sf) ;
      printf("x         = %e\n",xx) ;
      return(-1) ;
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Flux */
  flux(el,u) ;

  return(0) ;
}


int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  Symmetry_t sym = Element_GetSymmetry(el) ;
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double c[2] ;

  double Dc_hSDc_h2sf[2]        ;
  double Dc_ohSDc_h2sf[2]       ;
  double Dc_hsfSDc_h2sf[2]     ;
  double Dc_sfSDc_h2sf[2]      ;
  double Dc_caSDc_h2sf[2]       ;
  double Dc_cahsfSDc_h2sf[2]   ;
  double Dc_h3sio4SDc_h2sf[2]   ;
  double Dc_h4sio4SDc_h2sf[2]   ;
  double Dc_cah3sio4SDc_h2sf[2] ;
  double Dc_h2sio4SDc_h2sf[2]   ;
  double Dc_cah2sio4SDc_h2sf[2] ;
  double Dc_casfaqSDc_h2sf[2]  ;  
  double Dc_caohSDc_h2sf[2]     ;
  double Dc_csfSDc_h2sf[2]     ;
  
  double Dc_hSDzn_si_s[2]        ;
  double Dc_ohSDzn_si_s[2]       ;
  double Dc_hsfSDzn_si_s[2]     ;
  double Dc_sfSDzn_si_s[2]      ;
  double Dc_caSDzn_si_s[2]       ;
  double Dc_cahsfSDzn_si_s[2]   ;
  double Dc_h3sio4SDzn_si_s[2]   ;
  double Dc_h4sio4SDzn_si_s[2]   ;
  double Dc_cah3sio4SDzn_si_s[2] ;
  double Dc_h2sio4SDzn_si_s[2]   ;
  double Dc_cah2sio4SDzn_si_s[2] ;
  double Dc_casfaqSDzn_si_s[2]  ;
  double Dc_caohSDzn_si_s[2]     ;

  double Dc_hSDzn_ca_s[2]        ;
  double Dc_ohSDzn_ca_s[2]       ;
  double Dc_hsfSDzn_ca_s[2]     ;
  double Dc_sfSDzn_ca_s[2]      ;
  double Dc_caSDzn_ca_s[2]       ;
  double Dc_cahsfSDzn_ca_s[2]   ;
  double Dc_h3sio4SDzn_ca_s[2]   ;
  double Dc_h4sio4SDzn_ca_s[2]   ;
  double Dc_cah3sio4SDzn_ca_s[2] ;
  double Dc_h2sio4SDzn_ca_s[2]   ;
  double Dc_cah2sio4SDzn_ca_s[2] ;
  double Dc_casfaqSDzn_ca_s[2]  ;  
  double Dc_caohSDzn_ca_s[2]     ;
  double Dc_csfSDzn_ca_s[2]     ;
  
  double Dc_hSDc_k[2]         ;
  double Dc_ohSDc_k[2]        ;
  double Dc_hsfSDc_k[2]      ;
  double Dc_sfSDc_k[2]       ;
  double Dc_caSDc_k[2]        ;
  double Dc_cahsfSDc_k[2]    ;
  double Dc_h3sio4SDc_k[2]    ;
  double Dc_cah3sio4SDc_k[2]  ;
  double Dc_h2sio4SDc_k[2]    ;
  double Dc_cah2sio4SDc_k[2]  ;
  double Dc_casfaqSDc_k[2]   ;  
  double Dc_caohSDc_k[2]      ;
  
  double Dc_hSDc_cl[2]         ;
  double Dc_ohSDc_cl[2]        ;
  double Dc_hsfSDc_cl[2]      ;
  double Dc_sfSDc_cl[2]       ;
  double Dc_caSDc_cl[2]        ;
  double Dc_cahsfSDc_cl[2]    ;
  double Dc_h3sio4SDc_cl[2]    ;
  double Dc_cah3sio4SDc_cl[2]  ;
  double Dc_h2sio4SDc_cl[2]    ;
  double Dc_cah2sio4SDc_cl[2]  ;
  double Dc_casfaqSDc_cl[2]   ;  
  double Dc_caohSDc_cl[2]      ;
  double *u[MAX_NOEUDS] ;
  double *u_n[MAX_NOEUDS] ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
    u_n[i] = Element_GetPreviousNodalUnknown(el,i) ;
  }
  
  /*
    Donnees
  */
  phi0      = Element_GetProperty(el)[pm("porosite")] ;
  n_ca_ref  = Element_GetProperty(el)[pm("N_CH")] ;
  n_si_ref  = Element_GetProperty(el)[pm("N_Si")] ;
  c_h2sf_eq  = Element_GetProperty(el)[pm("C_H2SF_eq")] ;
  t_ch      = Element_GetProperty(el)[pm("T_CH")] ;
  t_csf      = Element_GetProperty(el)[pm("T_CSF")] ;


  /*
    Initialisation 
  */
  for(i=0;i<nn*nn*NEQ*NEQ;i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    CALCUL DE volume ET DE surf 
  */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    dx = x1 - x0 ;
    xm = (x1 + x0)*0.5 ;
  }
  for(i=0;i<nn;i++) {
    double x = Element_GetNodeCoordinate(el,i)[0] ;
    volume[i] = fabs(dx)*0.5 ; 
    if(sym == AXIS) volume[i] *= M_PI*(x + xm) ; 
  }
  if(sym == AXIS) surf = 2*M_PI*xm ; else surf = 1 ;
  /*
    termes d'accumulation
  */
  for(i=0;i<nn;i++) {
    /* molarities */
    double *x = ComputeAqueousVariables(el,u,i) ;

    double c_h2sf    = x[I_C_H2SF] ;
    double zn_si_s    = x[I_ZN_Si_S];
    double zn_ca_s    = x[I_ZN_Ca_S] ;
    double zc_h2sf   = c_h2sf/c_h2sf_eq ;
    double zq_ch      = x[I_ZQ_CH] ;

    
    /* solid contents : CH, CSF, CSH */
    double n_csfn    = N_CSFn(i) ;
    double n_chn      = N_CHn(i) ;

    
    /*   kinetics */
    double n_ch_ci    = n_chn*pow(zc_h2sf,-dt/t_ch) ;
    double n_csf_ci  = n_csfn*pow(zc_h2sf,dt/t_csf) ;
    double n_si_s     = N_Si_S(i) ;
    double x_csh      = X_CSH(zq_ch) ;

    /* porosity */
    double v_csh      = V_CSH(zq_ch) ;
#if (U_PHI == IMPLICIT)
    double phi        = PHI(i) ;
#else
    double phi        = PHIn(i) ;
#endif



    /* Atom concentrations in the liquid phase */
    double c_s_l      = x[I_C_S_L] ;
    double c_ca_l     = x[I_C_Ca_L] ;
    double c_si_l     = x[I_C_Si_L] ;
    double c_k_l      = x[I_C_K_L] ;
    double c_cl_l      = x[I_C_Cl_L] ;
    

    int j = i*NEQ ;

    /* derivatives ... */
    /* ... with respect to c_h2sf */
    {
      double dc_h2sf     = 1.e-8*c_h2sf*((c_h2sf > C_H2SFn(i)) ? 1 : -1) ;
      double *dcsdc_h2sf = ComputeAqueousDerivatives(el,x,dc_h2sf,I_C_H2SF) ;
		                    
      /* solid */
      double dn_ch_cisdc_h2sf   = n_ch_ci*(-dt/t_ch)/c_h2sf ;
      double dn_csf_cisdc_h2sf = n_csf_ci*(dt/t_csf)/c_h2sf ;
      double dn_chsdc_h2sf      = (zc_h2sf <= 1) ? -dn_csf_cisdc_h2sf  : dn_ch_cisdc_h2sf ;
      double dn_csfsdc_h2sf    = (zc_h2sf >  1) ? -dn_ch_cisdc_h2sf    : dn_csf_cisdc_h2sf ;
      
      double dzq_chsdc_h2sf     = dcsdc_h2sf[I_ZQ_CH] ;
      double dx_cshsdc_h2sf     = DX_CSH(zq_ch)*dzq_chsdc_h2sf ;
      
      double dn_ca_ssdc_h2sf    = dn_chsdc_h2sf + dn_csfsdc_h2sf \
                                 + dx_cshsdc_h2sf*n_si_s;
      double dn_si_ssdc_h2sf    = 0. ;
      double dn_s_ssdc_h2sf     = dn_csfsdc_h2sf ;
      
      /* porosity */
#if (U_PHI == IMPLICIT)
      double dv_cshsdc_h2sf     = DV_CSH(zq_ch)*dzq_chsdc_h2sf ;
      double dphisdc_h2sf       = - V_CH*dn_chsdc_h2sf - V_CSF*dn_csfsdc_h2sf \
                                 - dv_cshsdc_h2sf*n_si_s ;
      double dphi_csdc_h2sf     = - V_CH*dn_chsdc_h2sf - dv_cshsdc_h2sf*n_si_s ;
#else
      double dphisdc_h2sf       = 0 ;
      double dphi_csdc_h2sf     = 0 ;
#endif
       
       
      /* liquid */
      double dc_s_lsdc_h2sf     = dcsdc_h2sf[I_C_S_L] ;
      double dc_ca_lsdc_h2sf    = dcsdc_h2sf[I_C_Ca_L] ;
      double dc_si_lsdc_h2sf    = dcsdc_h2sf[I_C_Si_L] ;
      double dc_k_lsdc_h2sf     = dcsdc_h2sf[I_C_K_L] ;
      double dc_cl_lsdc_h2sf    = dcsdc_h2sf[I_C_Cl_L] ;
      double dn_s_lsdc_h2sf     = phi*dc_s_lsdc_h2sf  + dphisdc_h2sf*c_s_l ;
      double dn_ca_lsdc_h2sf    = phi*dc_ca_lsdc_h2sf + dphisdc_h2sf*c_ca_l ;
      double dn_si_lsdc_h2sf    = phi*dc_si_lsdc_h2sf + dphisdc_h2sf*c_si_l ;
      double dn_k_lsdc_h2sf     = phi*dc_k_lsdc_h2sf  + dphisdc_h2sf*c_k_l ;
      double dn_cl_lsdc_h2sf    = phi*dc_cl_lsdc_h2sf  + dphisdc_h2sf*c_cl_l ;
      
      /* Balance of S (sulfur)  : (n_S1 - n_Sn) + dt * div(w_S) = 0 */
      K(E_S+j,I_C_H2SF+j)   += volume[i]*(dn_s_lsdc_h2sf + dn_s_ssdc_h2sf) ;
      /* Balance of Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0 */
      K(E_Ca+j,I_C_H2SF+j)  += volume[i]*(dn_ca_lsdc_h2sf + dn_ca_ssdc_h2sf) ;
      /* Balance of Si (silicon) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0 */
      K(E_Si+j,I_C_H2SF+j)  += volume[i]*(dn_si_lsdc_h2sf + dn_si_ssdc_h2sf) ;
      /* Balance of K (potassium): (n_K1 - n_Kn) + dt * div(w_K) = 0 */
      K(E_K+j,I_C_H2SF+j)   += volume[i]*(dn_k_lsdc_h2sf) ;
      /* Balance of Cl (chloride): (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0 */
      K(E_Cl+j,I_C_H2SF+j)  += volume[i]*(dn_cl_lsdc_h2sf) ;
    

      /* sauvegardes pour les termes de transport */
      Dc_hSDc_h2sf[i]        = dcsdc_h2sf[I_C_H] ;
      Dc_ohSDc_h2sf[i]       = dcsdc_h2sf[I_C_OH] ;
      Dc_hsfSDc_h2sf[i]     = dcsdc_h2sf[I_C_HSF] ;
      Dc_sfSDc_h2sf[i]      = dcsdc_h2sf[I_C_SF] ;
      Dc_caSDc_h2sf[i]       = dcsdc_h2sf[I_C_Ca] ;
      Dc_cahsfSDc_h2sf[i]   = dcsdc_h2sf[I_C_CaHSF] ;
      Dc_h3sio4SDc_h2sf[i]   = dcsdc_h2sf[I_C_H3SiO4] ;
      Dc_h4sio4SDc_h2sf[i]   = dcsdc_h2sf[I_C_H4SiO4] ;
      Dc_cah3sio4SDc_h2sf[i] = dcsdc_h2sf[I_C_CaH3SiO4] ;
      Dc_h2sio4SDc_h2sf[i]   = dcsdc_h2sf[I_C_H2SiO4] ;
      Dc_cah2sio4SDc_h2sf[i] = dcsdc_h2sf[I_C_CaH2SiO4] ;
      Dc_casfaqSDc_h2sf[i]  = dcsdc_h2sf[I_C_CaSFaq] ;
      Dc_caohSDc_h2sf[i]     = dcsdc_h2sf[I_C_CaOH] ;
      Dc_csfSDc_h2sf[i]     = dn_csfsdc_h2sf;
    }
			                    
    /* with respect to zn_si_s */
    {
      double dzn_si_s     = 1.e-8*((zn_si_s > ZN_Si_Sn(i)) ? 1 : -1) ; 
      double *dcsdzn_si_s = ComputeAqueousDerivatives(el,x,dzn_si_s,I_ZN_Si_S) ;
 
      /* solid */
      double dn_si_ssdzn_si_s    = (zn_si_s > 0) ? n_si_ref : 0. ;
      double dn_ca_ssdzn_si_s    = x_csh*dn_si_ssdzn_si_s ;
    
      /* porosity */
#if (U_PHI == IMPLICIT)
      double dphisdzn_si_s       = -v_csh*dn_si_ssdzn_si_s ;
#else
      double dphisdzn_si_s       = 0 ;
#endif
    
      /* liquid */
      double dc_s_lsdzn_si_s     = dcsdzn_si_s[I_C_S_L] ;
      double dc_ca_lsdzn_si_s    = dcsdzn_si_s[I_C_Ca_L] ;
      double dc_si_lsdzn_si_s    = dcsdzn_si_s[I_C_Si_L] ;
      double dc_k_lsdzn_si_s     = dcsdzn_si_s[I_C_K_L] ; 
      double dc_cl_lsdzn_si_s    = dcsdzn_si_s[I_C_Cl_L] ;      
      double dn_s_lsdzn_si_s     = phi*dc_s_lsdzn_si_s  + dphisdzn_si_s*c_s_l ;
      double dn_ca_lsdzn_si_s    = phi*dc_ca_lsdzn_si_s + dphisdzn_si_s*c_ca_l ;
      double dn_si_lsdzn_si_s    = phi*dc_si_lsdzn_si_s + dphisdzn_si_s*c_si_l ;
      double dn_k_lsdzn_si_s     = phi*dc_k_lsdzn_si_s  + dphisdzn_si_s*c_k_l ;
      double dn_cl_lsdzn_si_s    = phi*dc_cl_lsdzn_si_s + dphisdzn_si_s*c_cl_l ;
    
      K(E_S+j,I_ZN_Si_S+j)   += volume[i]*(dn_s_lsdzn_si_s) ;
      K(E_Ca+j,I_ZN_Si_S+j)  += volume[i]*(dn_ca_lsdzn_si_s + dn_ca_ssdzn_si_s) ;
      K(E_Si+j,I_ZN_Si_S+j)  += volume[i]*(dn_si_lsdzn_si_s + dn_si_ssdzn_si_s) ;
      K(E_K+j,I_ZN_Si_S+j)   += volume[i]*(dn_k_lsdzn_si_s) ;
      K(E_Cl+j,I_ZN_Si_S+j)  += volume[i]*(dn_cl_lsdzn_si_s) ;
    
      /* sauvegardes pour les termes de transport */
      Dc_hSDzn_si_s[i]        = dcsdzn_si_s[I_C_H] ;
      Dc_ohSDzn_si_s[i]       = dcsdzn_si_s[I_C_OH] ;
      Dc_hsfSDzn_si_s[i]     = dcsdzn_si_s[I_C_HSF] ;
      Dc_sfSDzn_si_s[i]      = dcsdzn_si_s[I_C_SF] ;   
      Dc_caSDzn_si_s[i]       = dcsdzn_si_s[I_C_Ca] ;
      Dc_cahsfSDzn_si_s[i]   = dcsdzn_si_s[I_C_CaHSF] ;
      Dc_h3sio4SDzn_si_s[i]   = dcsdzn_si_s[I_C_H3SiO4] ;
      Dc_h4sio4SDzn_si_s[i]   = dcsdzn_si_s[I_C_H4SiO4] ;
      Dc_cah3sio4SDzn_si_s[i] = dcsdzn_si_s[I_C_CaH3SiO4] ;
      Dc_h2sio4SDzn_si_s[i]   = dcsdzn_si_s[I_C_H2SiO4] ;
      Dc_cah2sio4SDzn_si_s[i] = dcsdzn_si_s[I_C_CaH2SiO4] ;
      Dc_casfaqSDzn_si_s[i]  = dcsdzn_si_s[I_C_CaSFaq] ;
      Dc_caohSDzn_si_s[i]     = dcsdzn_si_s[I_C_CaOH] ;  
    }
    
    /* with respect to zn_ca_s */
    {
      double dzn_ca_s            = 1.e-10*((zn_ca_s > ZN_Ca_Sn(i)) ? 1 : -1) ;
      double *dcsdzn_ca_s = ComputeAqueousDerivatives(el,x,dzn_ca_s,I_ZN_Ca_S) ;

      /* solid */
      double dzq_chsdzn_ca_s    = dcsdzn_ca_s[I_ZQ_CH] ;
      double dx_cshsdzn_ca_s    = DX_CSH(zq_ch)*dzq_chsdzn_ca_s ;     

      double dn_ch_eqsdzn_ca_s   = (zn_ca_s > 0) ? n_ca_ref : 0 ;
      double dn_csf_eqsdzn_ca_s = (zn_ca_s > 0) ? n_ca_ref : 0 ;
      double dn_chsdzn_ca_s      = (zc_h2sf <= 1) ? dn_ch_eqsdzn_ca_s    : 0 ;
      double dn_csfsdzn_ca_s    = (zc_h2sf >  1) ? dn_csf_eqsdzn_ca_s  : 0 ;
    
      double dn_si_ssdzn_ca_s    = 0. ;
      double dn_ca_ssdzn_ca_s    = dn_chsdzn_ca_s + dn_csfsdzn_ca_s \
                                 + dx_cshsdzn_ca_s*n_si_s ;
      double dn_s_ssdzn_ca_s     = dn_csfsdzn_ca_s ;

      /* porosity */
#if (U_PHI == IMPLICIT)
      double dv_cshsdzn_ca_s     = DV_CSH(zq_ch)*dzq_chsdzn_ca_s ;
      double dphisdzn_ca_s       = - V_CH*dn_chsdzn_ca_s \
                                 - V_CSF*dn_csfsdzn_ca_s \
                                 - dv_cshsdzn_ca_s*n_si_s ;
#else
      double dphisdzn_ca_s       = 0 ;
#endif

      /* liquid */
      double dc_s_lsdzn_ca_s     = dcsdzn_ca_s[I_C_S_L] ;
      double dc_ca_lsdzn_ca_s    = dcsdzn_ca_s[I_C_Ca_L] ;
      double dc_si_lsdzn_ca_s    = dcsdzn_ca_s[I_C_Si_L] ;
      double dc_k_lsdzn_ca_s     = dcsdzn_ca_s[I_C_K_L] ;   
      double dc_cl_lsdzn_ca_s    = dcsdzn_ca_s[I_C_Cl_L] ; 
      double dn_s_lsdzn_ca_s     = phi*dc_s_lsdzn_ca_s  + dphisdzn_ca_s*c_s_l ;
      double dn_ca_lsdzn_ca_s    = phi*dc_ca_lsdzn_ca_s + dphisdzn_ca_s*c_ca_l ;
      double dn_si_lsdzn_ca_s    = phi*dc_si_lsdzn_ca_s + dphisdzn_ca_s*c_si_l ;
      double dn_k_lsdzn_ca_s     = phi*dc_k_lsdzn_ca_s  + dphisdzn_ca_s*c_k_l ;
      double dn_cl_lsdzn_ca_s    = phi*dc_cl_lsdzn_ca_s + dphisdzn_ca_s*c_cl_l ;
    
      K(E_S+j,I_ZN_Ca_S+j)   += volume[i]*(dn_s_lsdzn_ca_s + dn_s_ssdzn_ca_s) ;
      K(E_Ca+j,I_ZN_Ca_S+j)  += volume[i]*(dn_ca_lsdzn_ca_s + dn_ca_ssdzn_ca_s) ;
      K(E_Si+j,I_ZN_Ca_S+j)  += volume[i]*(dn_si_lsdzn_ca_s + dn_si_ssdzn_ca_s) ;
      K(E_K+j,I_ZN_Ca_S+j)   += volume[i]*(dn_k_lsdzn_ca_s) ; 
      K(E_Cl+j,I_ZN_Ca_S+j)  += volume[i]*(dn_cl_lsdzn_ca_s) ;
      
      /* sauvegardes pour les termes de transport */
      Dc_hSDzn_ca_s[i]        = dcsdzn_ca_s[I_C_H] ;
      Dc_ohSDzn_ca_s[i]       = dcsdzn_ca_s[I_C_OH] ;
      Dc_hsfSDzn_ca_s[i]     = dcsdzn_ca_s[I_C_HSF] ;
      Dc_sfSDzn_ca_s[i]      = dcsdzn_ca_s[I_C_SF] ;
      Dc_caSDzn_ca_s[i]       = dcsdzn_ca_s[I_C_Ca] ;      
      Dc_cahsfSDzn_ca_s[i]   = dcsdzn_ca_s[I_C_CaHSF] ;
      Dc_h3sio4SDzn_ca_s[i]   = dcsdzn_ca_s[I_C_H3SiO4] ;
      Dc_h4sio4SDzn_ca_s[i]   = dcsdzn_ca_s[I_C_H4SiO4] ;
      Dc_cah3sio4SDzn_ca_s[i] = dcsdzn_ca_s[I_C_CaH3SiO4] ;
      Dc_h2sio4SDzn_ca_s[i]   = dcsdzn_ca_s[I_C_H2SiO4] ;
      Dc_cah2sio4SDzn_ca_s[i] = dcsdzn_ca_s[I_C_CaH2SiO4] ;
      Dc_casfaqSDzn_ca_s[i]  = dcsdzn_ca_s[I_C_CaSFaq] ;
      Dc_caohSDzn_ca_s[i]     = dcsdzn_ca_s[I_C_CaOH] ;
      Dc_csfSDzn_ca_s[i]     = dn_csfsdzn_ca_s ;
    }
    
    /* ... with respect to c_k */
    {
      double dc_k                = 1.e-6 ;
      double *dcsdc_k = ComputeAqueousDerivatives(el,x,dc_k,I_C_K) ;
      
      /* liquid */
      double dc_s_lsdc_k     = dcsdc_k[I_C_S_L] ;
      double dc_ca_lsdc_k    = dcsdc_k[I_C_Ca_L] ;
      double dc_si_lsdc_k    = dcsdc_k[I_C_Si_L] ;
      double dc_k_lsdc_k     = dcsdc_k[I_C_K_L] ;
      double dc_cl_lsdc_k    = dcsdc_k[I_C_Cl_L] ;
      double dn_s_lsdc_k     = phi*dc_s_lsdc_k ;
      double dn_ca_lsdc_k    = phi*dc_ca_lsdc_k ;
      double dn_si_lsdc_k    = phi*dc_si_lsdc_k ;
      double dn_k_lsdc_k     = phi*dc_k_lsdc_k ;
      double dn_cl_lsdc_k    = phi*dc_cl_lsdc_k ;
            
      K(E_S+j,I_C_K+j)     += volume[i]*(dn_s_lsdc_k) ; 
      K(E_Ca+j,I_C_K+j)    += volume[i]*(dn_ca_lsdc_k) ;
      K(E_Si+j,I_C_K+j)    += volume[i]*(dn_si_lsdc_k) ; 
      K(E_K+j,I_C_K+j)     += volume[i]*(dn_k_lsdc_k) ;
      K(E_Cl+j,I_C_K+j)    += volume[i]*(dn_cl_lsdc_k) ;
      
      /* sauvegardes pour les termes de transport */
      Dc_hSDc_k[i]        = dcsdc_k[I_C_H] ;
      Dc_ohSDc_k[i]       = dcsdc_k[I_C_OH] ;
      Dc_hsfSDc_k[i]     = dcsdc_k[I_C_HSF] ;
      Dc_sfSDc_k[i]      = dcsdc_k[I_C_SF] ;
      Dc_caSDc_k[i]       = dcsdc_k[I_C_Ca] ;      
      Dc_cahsfSDc_k[i]   = dcsdc_k[I_C_CaHSF] ;
      Dc_h3sio4SDc_k[i]   = dcsdc_k[I_C_H3SiO4] ;
      Dc_cah3sio4SDc_k[i] = dcsdc_k[I_C_CaH3SiO4] ;
      Dc_h2sio4SDc_k[i]   = dcsdc_k[I_C_H2SiO4] ;
      Dc_cah2sio4SDc_k[i] = dcsdc_k[I_C_CaH2SiO4] ;
      Dc_casfaqSDc_k[i]  = dcsdc_k[I_C_CaSFaq] ;
      Dc_caohSDc_k[i]     = dcsdc_k[I_C_CaOH] ;    
    }
    
        /* ... with respect to c_cl */
    {
      double dc_cl                = 1.e-6 ;
      double *dcsdc_cl = ComputeAqueousDerivatives(el,x,dc_cl,I_C_Cl) ;
      
      /* liquid */
      double dc_s_lsdc_cl     = dcsdc_cl[I_C_S_L] ;
      double dc_ca_lsdc_cl    = dcsdc_cl[I_C_Ca_L] ;
      double dc_si_lsdc_cl    = dcsdc_cl[I_C_Si_L] ;
      double dc_k_lsdc_cl     = dcsdc_cl[I_C_K_L] ;
      double dc_cl_lsdc_cl    = dcsdc_cl[I_C_Cl_L] ;
      double dn_s_lsdc_cl     = phi*dc_s_lsdc_cl ;
      double dn_ca_lsdc_cl    = phi*dc_ca_lsdc_cl ;
      double dn_si_lsdc_cl    = phi*dc_si_lsdc_cl ;
      double dn_k_lsdc_cl     = phi*dc_k_lsdc_cl ;
      double dn_cl_lsdc_cl    = phi*dc_cl_lsdc_cl ;
      
      K(E_S+j,I_C_Cl+j)     += volume[i]*(dn_s_lsdc_cl) ; 
      K(E_Ca+j,I_C_Cl+j)    += volume[i]*(dn_ca_lsdc_cl) ;
      K(E_Si+j,I_C_Cl+j)    += volume[i]*(dn_si_lsdc_cl) ; 
      K(E_K+j,I_C_Cl+j)     += volume[i]*(dn_k_lsdc_cl) ;
      K(E_Cl+j,I_C_Cl+j)    += volume[i]*(dn_cl_lsdc_cl) ;
    
      /* sauvegardes pour les termes de transport */
      Dc_hSDc_cl[i]        = dcsdc_cl[I_C_H] ;
      Dc_ohSDc_cl[i]       = dcsdc_cl[I_C_OH] ;
      Dc_hsfSDc_cl[i]     = dcsdc_cl[I_C_HSF] ;
      Dc_sfSDc_cl[i]      = dcsdc_cl[I_C_SF] ;
      Dc_caSDc_cl[i]       = dcsdc_cl[I_C_Ca] ;      
      Dc_cahsfSDc_cl[i]   = dcsdc_cl[I_C_CaHSF] ;
      Dc_h3sio4SDc_cl[i]   = dcsdc_cl[I_C_H3SiO4] ;
      Dc_cah3sio4SDc_cl[i] = dcsdc_cl[I_C_CaH3SiO4] ;
      Dc_h2sio4SDc_cl[i]   = dcsdc_cl[I_C_H2SiO4] ;
      Dc_cah2sio4SDc_cl[i] = dcsdc_cl[I_C_CaH2SiO4] ;
      Dc_casfaqSDc_cl[i]  = dcsdc_cl[I_C_CaSFaq] ;
      Dc_caohSDc_cl[i]     = dcsdc_cl[I_C_CaOH] ;    
    }

  }

  /* termes d'ecoulement */
  {
  double tr        = dt*surf/dx ;

  double trf_oh       = tr*KF_OH ;
  double trf_h        = tr*KF_H ;
  double trf_h2sf    = tr*KF_H2SF ;
  double trf_hsf     = tr*KF_HSF ;
  double trf_sf      = tr*KF_SF ;
  double trf_ca       = tr*KF_Ca ;
  double trf_cahsf   = tr*KF_CaHSF ;
  double trf_cah3sio4 = tr*KF_CaH3SiO4 ;
  double trf_h3sio4   = tr*KF_H3SiO4 ;
  double trf_h4sio4   = tr*KF_H4SiO4 ;
  double trf_h2sio4   = tr*KF_H2SiO4 ;
  double trf_cah2sio4 = tr*KF_CaH2SiO4 ;
  double trf_casfaq  = tr*KF_CaSFaq ;
  double trf_caoh     = tr*KF_CaOH ;
  double trf_k        = tr*KF_K ;
  double trf_cl       = tr*KF_Cl ;
  
  double tre_hsf     = tr*Kpsi_HSF ;
  double tre_sf      = tr*Kpsi_SF ;
  double tre_ca       = tr*Kpsi_Ca ;
  double tre_cahsf   = tr*Kpsi_CaHSF ;
  double tre_cah3sio4 = tr*Kpsi_CaH3SiO4 ;
  double tre_h3sio4   = tr*Kpsi_H3SiO4 ;
  double tre_h2sio4   = tr*Kpsi_H2SiO4 ;
  double tre_caoh     = tr*Kpsi_CaOH ;
  double tre_k        = tr*Kpsi_K ;
  double tre_cl       = tr*Kpsi_Cl ;

  double tre_q        = tr*Kpsi_q ;
  
  double trd_csf 	  = tr*KD_CSF ;
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h2sf + trf_hsf*Dc_hsfSDc_h2sf[i] + trf_sf*Dc_sfSDc_h2sf[i] + trf_cahsf*Dc_cahsfSDc_h2sf[i] \
           + trf_casfaq*Dc_casfaqSDc_h2sf[i] + trd_csf*Dc_csfSDc_h2sf[i]  ;
  }
  K(E_S,I_C_H2SF)          += + c[0] ;
  K(E_S,I_C_H2SF+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_C_H2SF)      += - c[0] ;
  K(E_S+NEQ,I_C_H2SF+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hsf*Dc_hsfSDzn_si_s[i] + trf_sf*Dc_sfSDzn_si_s[i] + trf_cahsf*Dc_cahsfSDzn_si_s[i] + trf_casfaq*Dc_casfaqSDzn_si_s[i] ;
  }
  K(E_S,I_ZN_Si_S)          += + c[0] ;
  K(E_S,I_ZN_Si_S+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_ZN_Si_S)      += - c[0] ;
  K(E_S+NEQ,I_ZN_Si_S+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hsf*Dc_hsfSDzn_ca_s[i] + trf_sf*Dc_sfSDzn_ca_s[i] + trf_cahsf*Dc_cahsfSDzn_ca_s[i] + trf_casfaq*Dc_casfaqSDzn_ca_s[i] + trd_csf*Dc_csfSDzn_ca_s[i]   ;
  }
  K(E_S,I_ZN_Ca_S)           += + c[0] ;
  K(E_S,I_ZN_Ca_S+NEQ)       += - c[1] ;
  K(E_S+NEQ,I_ZN_Ca_S)       += - c[0] ;
  K(E_S+NEQ,I_ZN_Ca_S+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_hsf*Dc_hsfSDc_k[i] + trf_sf*Dc_sfSDc_k[i] + trf_cahsf*Dc_cahsfSDc_k[i] + trf_casfaq*Dc_casfaqSDc_k[i] ;
  }
  K(E_S,I_C_K)          += + c[0] ;
  K(E_S,I_C_K+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_C_K)      += - c[0] ;
  K(E_S+NEQ,I_C_K+NEQ)  += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_hsf*Dc_hsfSDc_cl[i] + trf_sf*Dc_sfSDc_cl[i] + trf_cahsf*Dc_cahsfSDc_cl[i] + trf_casfaq*Dc_casfaqSDc_cl[i] ;
  }
  K(E_S,I_C_Cl)          += + c[0] ;
  K(E_S,I_C_Cl+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_C_Cl)      += - c[0] ;
  K(E_S+NEQ,I_C_Cl+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hsf + tre_sf + tre_cahsf ;
  }
  K(E_S,I_PSI)          += + c[0] ;
  K(E_S,I_PSI+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_PSI)      += - c[0] ;
  K(E_S+NEQ,I_PSI+NEQ)  += + c[1] ;

  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_h2sf[i] + trf_cahsf*Dc_cahsfSDc_h2sf[i] + trf_cah3sio4*Dc_cah3sio4SDc_h2sf[i] + trf_cah2sio4*Dc_cah2sio4SDc_h2sf[i] \
    + trf_caoh*Dc_caohSDc_h2sf[i] + trf_casfaq*Dc_casfaqSDc_h2sf[i] + trd_csf*Dc_csfSDc_h2sf[i]   ;
  }
  K(E_Ca,I_C_H2SF)         += + c[0] ;
  K(E_Ca,I_C_H2SF+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_C_H2SF)     += - c[0] ;
  K(E_Ca+NEQ,I_C_H2SF+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_si_s[i] + trf_cahsf*Dc_cahsfSDzn_si_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + trf_caoh*Dc_caohSDzn_si_s[i] + trf_casfaq*Dc_casfaqSDzn_si_s[i] ;
  }
  K(E_Ca,I_ZN_Si_S)         += + c[0] ;
  K(E_Ca,I_ZN_Si_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_ZN_Si_S)     += - c[0] ;
  K(E_Ca+NEQ,I_ZN_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_ca_s[i] + trf_cahsf*Dc_cahsfSDzn_ca_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + trf_caoh*Dc_caohSDzn_ca_s[i] + trf_casfaq*Dc_casfaqSDzn_ca_s[i] + trd_csf*Dc_csfSDzn_ca_s[i]   ;
  }
  K(E_Ca,I_ZN_Ca_S)         += + c[0] ;
  K(E_Ca,I_ZN_Ca_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_ZN_Ca_S)     += - c[0] ;
  K(E_Ca+NEQ,I_ZN_Ca_S+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_k[i] + trf_cahsf*Dc_cahsfSDc_k[i] + trf_cah3sio4*Dc_cah3sio4SDc_k[i] + trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + trf_caoh*Dc_caohSDc_k[i] + trf_casfaq*Dc_casfaqSDc_k[i] ;
  }
  K(E_Ca,I_C_K)         += + c[0] ;
  K(E_Ca,I_C_K+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_C_K)     += - c[0] ;
  K(E_Ca+NEQ,I_C_K+NEQ) += + c[1] ; 

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_cl[i] + trf_cahsf*Dc_cahsfSDc_cl[i] + trf_cah3sio4*Dc_cah3sio4SDc_cl[i] + trf_cah2sio4*Dc_cah2sio4SDc_cl[i] \
    + trf_caoh*Dc_caohSDc_cl[i] + trf_casfaq*Dc_casfaqSDc_cl[i] ;
  }
  K(E_Ca,I_C_Cl)         += + c[0] ;
  K(E_Ca,I_C_Cl+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_C_Cl)     += - c[0] ;
  K(E_Ca+NEQ,I_C_Cl+NEQ) += + c[1] ; 
  
  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_cahsf + tre_cah3sio4 + tre_caoh;
  }
  K(E_Ca,I_PSI)          += + c[0] ;
  K(E_Ca,I_PSI+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_PSI)      += - c[0] ;
  K(E_Ca+NEQ,I_PSI+NEQ)  += + c[1] ;
  /*
    Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_h2sf[i] + trf_h4sio4*Dc_h4sio4SDc_h2sf[i] + trf_cah3sio4*Dc_cah3sio4SDc_h2sf[i] + trf_cah2sio4*Dc_cah2sio4SDc_h2sf[i] \
    + trf_h2sio4*Dc_h2sio4SDc_h2sf[i] ;
  }
  K(E_Si,I_C_H2SF)         += + c[0] ;
  K(E_Si,I_C_H2SF+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_C_H2SF)     += - c[0] ;
  K(E_Si+NEQ,I_C_H2SF+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDzn_si_s[i] + trf_h4sio4*Dc_h4sio4SDzn_si_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i]  + trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + trf_h2sio4*Dc_h2sio4SDzn_si_s[i];
  }
  K(E_Si,I_ZN_Si_S)         += + c[0] ;
  K(E_Si,I_ZN_Si_S+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_ZN_Si_S)     += - c[0] ;
  K(E_Si+NEQ,I_ZN_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] + trf_h4sio4*Dc_h4sio4SDzn_ca_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i]  + trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + trf_h2sio4*Dc_h2sio4SDzn_ca_s[i];
  }
  K(E_Si,I_ZN_Ca_S)         += + c[0] ;
  K(E_Si,I_ZN_Ca_S+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_ZN_Ca_S)     += - c[0] ;
  K(E_Si+NEQ,I_ZN_Ca_S+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_k[i] + trf_cah3sio4*Dc_cah3sio4SDc_k[i]  + trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + trf_h2sio4*Dc_h2sio4SDc_k[i];
  }
  K(E_Si,I_C_K)         += + c[0] ;
  K(E_Si,I_C_K+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_C_K)     += - c[0] ;
  K(E_Si+NEQ,I_C_K+NEQ) += + c[1] ; 
  
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_cl[i] + trf_cah3sio4*Dc_cah3sio4SDc_cl[i]  + trf_cah2sio4*Dc_cah2sio4SDc_cl[i] \
    + trf_h2sio4*Dc_h2sio4SDc_cl[i];
  }
  K(E_Si,I_C_Cl)         += + c[0] ;
  K(E_Si,I_C_Cl+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_C_Cl)     += - c[0] ;
  K(E_Si+NEQ,I_C_Cl+NEQ) += + c[1] ; 

  for(i=0;i<2;i++){
    c[i] = tre_h3sio4 + tre_cah3sio4 + tre_h2sio4 ;
  }
  K(E_Si,I_PSI)          += + c[0] ;
  K(E_Si,I_PSI+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_PSI)      += - c[0] ;
  K(E_Si+NEQ,I_PSI+NEQ)  += + c[1] ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDc_h2sf[i] + z_oh*trf_oh*Dc_ohSDc_h2sf[i] + z_hsf*trf_hsf*Dc_hsfSDc_h2sf[i] + z_sf*trf_sf*Dc_sfSDc_h2sf[i] \
           + z_ca*trf_ca*Dc_caSDc_h2sf[i] + z_cahsf*trf_cahsf*Dc_cahsfSDc_h2sf[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_h2sf[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_h2sf[i] + z_caoh*trf_caoh*Dc_caohSDc_h2sf[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_h2sf[i];
  }
  K(E_q,I_C_H2SF)           += + c[0] ;
  K(E_q,I_C_H2SF+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_C_H2SF)       += - c[0] ;
  K(E_q+NEQ,I_C_H2SF+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_si_s[i] + z_oh*trf_oh*Dc_ohSDzn_si_s[i] + z_hsf*trf_hsf*Dc_hsfSDzn_si_s[i] + z_sf*trf_sf*Dc_sfSDzn_si_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_si_s[i] + z_cahsf*trf_cahsf*Dc_cahsfSDzn_si_s[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDzn_si_s[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + z_caoh*trf_caoh*Dc_caohSDzn_si_s[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDzn_si_s[i];
  }
  K(E_q,I_ZN_Si_S)           += + c[0] ;
  K(E_q,I_ZN_Si_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_ZN_Si_S)       += - c[0] ;
  K(E_q+NEQ,I_ZN_Si_S+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_ca_s[i] + z_oh*trf_oh*Dc_ohSDzn_ca_s[i] + z_hsf*trf_hsf*Dc_hsfSDzn_ca_s[i] + z_sf*trf_sf*Dc_sfSDzn_ca_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_ca_s[i] + z_cahsf*trf_cahsf*Dc_cahsfSDzn_ca_s[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + z_caoh*trf_caoh*Dc_caohSDzn_ca_s[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDzn_ca_s[i];
  }
  K(E_q,I_ZN_Ca_S)           += + c[0] ;
  K(E_q,I_ZN_Ca_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_ZN_Ca_S)       += - c[0] ;
  K(E_q+NEQ,I_ZN_Ca_S+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = z_k*trf_k + z_h*trf_h*Dc_hSDc_k[i] + z_oh*trf_oh*Dc_ohSDc_k[i] + z_hsf*trf_hsf*Dc_hsfSDc_k[i] + z_sf*trf_sf*Dc_sfSDc_k[i] \
           + z_ca*trf_ca*Dc_caSDc_k[i] + z_cahsf*trf_cahsf*Dc_cahsfSDc_k[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_k[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_k[i] + z_caoh*trf_caoh*Dc_caohSDc_k[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_k[i];
  }
  K(E_q,I_C_K)           += + c[0] ;
  K(E_q,I_C_K+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_C_K)       += - c[0] ;
  K(E_q+NEQ,I_C_K+NEQ)   += + c[1] ; 
  
  for(i=0;i<2;i++){
    c[i] = z_cl*trf_cl + z_h*trf_h*Dc_hSDc_cl[i] + z_oh*trf_oh*Dc_ohSDc_cl[i] + z_hsf*trf_hsf*Dc_hsfSDc_cl[i] + z_sf*trf_sf*Dc_sfSDc_cl[i] \
           + z_ca*trf_ca*Dc_caSDc_cl[i] + z_cahsf*trf_cahsf*Dc_cahsfSDc_cl[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_cl[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_cl[i] + z_caoh*trf_caoh*Dc_caohSDc_cl[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_cl[i];
  }
  K(E_q,I_C_Cl)           += + c[0] ;
  K(E_q,I_C_Cl+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_C_Cl)       += - c[0] ;
  K(E_q+NEQ,I_C_Cl+NEQ)   += + c[1] ; 

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_PSI)          += + c[0] ;
  K(E_q,I_PSI+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_PSI)      += - c[0] ;
  K(E_q+NEQ,I_PSI+NEQ)  += + c[1] ;
  
    /*
  Conservation de K (potassium)  : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  for(i=0;i<2;i++){
    c[i] = tre_k ;
  }
  K(E_K,I_PSI)          += + c[0] ;
  K(E_K,I_PSI+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_PSI)      += - c[0] ;
  K(E_K+NEQ,I_PSI+NEQ)  += + c[1] ;
    
   
  for(i=0;i<2;i++){
    c[i] = trf_k ;
  }
  K(E_K,I_C_K)          += + c[0] ;
  K(E_K,I_C_K+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_C_K)      += - c[0] ;
  K(E_K+NEQ,I_C_K+NEQ)  += + c[1] ; 
  
     /*
  Conservation de Cl (chloride)  : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  for(i=0;i<2;i++){
    c[i] = tre_cl ;
  }
  K(E_Cl,I_PSI)          += + c[0] ;
  K(E_Cl,I_PSI+NEQ)      += - c[1] ;
  K(E_Cl+NEQ,I_PSI)      += - c[0] ;
  K(E_Cl+NEQ,I_PSI+NEQ)  += + c[1] ;
    
   
  for(i=0;i<2;i++){
    c[i] = trf_cl ;
  }
  K(E_Cl,I_C_Cl)          += + c[0] ;
  K(E_Cl,I_C_Cl+NEQ)      += - c[1] ;
  K(E_Cl+NEQ,I_C_Cl)      += - c[0] ;
  K(E_Cl+NEQ,I_C_Cl+NEQ)  += + c[1] ; 
 }

#if (U_H2SF == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_C_H2SF)     *= Ln10*C_H2SF(0) ;
    K(i,I_C_H2SF+NEQ) *= Ln10*C_H2SF(1) ;
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
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  R(0,E_S) -= volume[0]*(N_S(0) - N_Sn(0)) + dt*surf*W_S + dt*surf*A_H2S;
  R(1,E_S) -= volume[1]*(N_S(1) - N_Sn(1)) - dt*surf*W_S - dt*surf*A_H2S;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;
  /*
      Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ; 
  
    /*
      Conservation de Cl (chloride) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*W_Cl ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*W_Cl ; 
  /*
    Conservation de la charge  : div(w_q) = 0
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
  double *u[Element_MaxNbOfNodes] ;
  int    i,j,nso ;
  double *h_s ;
  double zero = 0. ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }

  /* if(Element_IsSubmanifold(el)) return(0) ; */
  
  /*
    Donnees
  */
  c_h2sf_eq  = Element_GetProperty(el)[pm("C_H2SF_eq")] ;

  /* initialisation */
  nso = 31 ;
  for(i = 0 ; i < nso ; i++) for(j = 0 ; j < 9 ; j++) Result_GetValue(r+i)[j] = zero ;


  /* output quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    double *x = ComputeAqueousVariables(el,u,j) ;
    
    double c_h2sf    = x[I_C_H2SF] ;
    double zn_si_s    = x[I_ZN_Si_S] ;
    double zn_ca_s    = x[I_ZN_Ca_S] ;
    double c_k        = x[I_C_K] ;
    double c_cl       = x[I_C_Cl] ;
    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hsf     = x[I_C_HSF] ;
    double c_sf      = x[I_C_SF] ;
    double c_ca       = x[I_C_Ca] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_h4sio4   = x[I_C_H4SiO4] ;
    double c_cah2sio4 = x[I_C_CaH2SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_cahsf   = x[I_C_CaHSF] ;
    double c_casfaq  = x[I_C_CaSFaq] ;
    double c_caoh     = x[I_C_CaOH] ;
    double zq_ch      = x[I_ZQ_CH] ;
    double zc_h2sf   = c_h2sf/c_h2sf_eq ;
    
    /* charge density */
    double c_q = x[I_N_Q] ;
    /* solid contents */
    double n_ch       = N_CH(j) ;
    double n_csf     = N_CSF(j) ;
    double n_si_s     = N_Si_S(j)  ;
    double a_h2s      = A_H2S ;

    
    /* porosity */
    double v_csh      = V_CSH(zq_ch) ;
    double phi        = PHI(j) ;

    double psi        = PSI(j) ;
    double ph         = 14 + log(c_oh)/log(10.) ;
    double pk_ch      = zq_ch ;
    double x_csh      = X_CSH(zq_ch) ;
    double test       = W_S;
    /*test       = Curve_ComputeValue(Element_GetCurve(el) + 4,ph);*/


    i = 0 ;
    Result_Store(r + i++,&ph,"ph",1) ;
    Result_Store(r + i++,&c_h2sf,"c_h2sf",1) ;
    Result_Store(r + i++,&c_hsf,"c_hsf",1) ;
    Result_Store(r + i++,&c_sf,"c_sf",1) ;
    Result_Store(r + i++,&c_ca,"c_ca",1) ;
    Result_Store(r + i++,&c_caoh,"c_caoh",1) ;
    Result_Store(r + i++,&c_h2sio4,"c_h2sio4",1) ;
    Result_Store(r + i++,&c_h3sio4,"c_h3sio4",1) ;
    Result_Store(r + i++,&c_h4sio4,"c_h4sio4",1) ;
    Result_Store(r + i++,&c_cah2sio4,"c_cah2sio4",1) ;
    Result_Store(r + i++,&c_cah3sio4,"c_cah3sio4",1) ;
    Result_Store(r + i++,&c_casfaq,"c_casfaq",1) ;
    Result_Store(r + i++,&c_cahsf,"c_cahsf",1) ;
    Result_Store(r + i++,&c_k,"c_k",1) ;
    Result_Store(r + i++,&c_cl,"c_cl",1) ;
    Result_Store(r + i++,&c_oh,"c_oh",1) ;
    Result_Store(r + i++,&c_h,"c_h",1) ;
    Result_Store(r + i++,&n_ch,"n_ch",1) ;
    Result_Store(r + i++,&n_csf,"n_csf",1) ;
    Result_Store(r + i++,&n_si_s,"n_si_s",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&psi,"potentiel_electrique",1) ;
    Result_Store(r + i++,&c_q,"charge",1) ;
    Result_Store(r + i++,&zn_ca_s,"zn_ca_s",1) ;
    Result_Store(r + i++,&zc_h2sf,"zc_h2sf",1) ;
    Result_Store(r + i++,&pk_ch,"pk_ch",1) ;
    Result_Store(r + i++,&zn_si_s,"zn_si_s",1) ;
    Result_Store(r + i++,&v_csh,"V_CSH",1) ;
    Result_Store(r + i++,&x_csh,"C/S",1) ;
    Result_Store(r + i++,&a_h2s,"A_H2S",1) ;
    Result_Store(r + i++,&test,"TEST",1) ;
  }
  
  if(i != nso) arret("ComputeOutputs (M70b)") ;

  return(nso) ;
}


void transfert(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int i ;
  /*
    Donnees
  */
  phi0      = Element_GetProperty(el)[pm("porosite")] ;
  c_h2sf_eq  = Element_GetProperty(el)[pm("C_H2SF_eq")] ;
  /*k_int     = Element_GetProperty(el)[pm("k_int")] ;*/
  /* initialisation */
  for(i=0;i<NVE_TR;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<2;i++) {
    /* molarities */
    double *x = ComputeAqueousVariables(el,u,i) ;

    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hsf     = x[I_C_HSF] ;
    double c_sf      = x[I_C_SF] ;
    double c_ca       = x[I_C_Ca] ;
    double c_cahsf   = x[I_C_CaHSF] ;
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_caoh     = x[I_C_CaOH] ;
    double c_k        = x[I_C_K] ;
    double c_cl       = x[I_C_Cl] ;

    /* porosity */
    double phi        = PHI(i) ;
    
    /* Gypsum contents */
    double n_csf     = N_CSF(i) ;
    
    /* tortuosite liquide */
    /*double iff    = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;*/
    double iff    = 0.0002 ;
    /* permeabilite */
    /* double k_g    = (k_int/mu_csf)*pow(phi_c/phi0,3.)*pow(((1-phi0)/(1-phi_c)),2.) ; */
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_H2SF      += d_h2sf*iff;
    KF_HSF       += d_hsf*iff;
    KF_SF        += d_sf*iff;

    KF_Ca         += d_ca*iff ;
    KF_CaHSF     += d_cahsf*iff;
    KF_CaH3SiO4   += d_cah3sio4*iff;

    KF_H3SiO4     += d_h3sio4*iff ;
    KF_H4SiO4     += d_h4sio4*iff ;
    KF_H2SiO4     += d_h2sio4*iff ;
    KF_CaH2SiO4   += d_cah2sio4*iff ;
    KF_CaSFaq    += d_casfaq*iff;
    KF_CaOH       += d_caoh*iff ;
    
    KF_K          += d_k*iff;
    KF_Cl         += d_cl*iff;


    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HSF     += FsRT*KF_HSF*z_hsf*c_hsf ;
    Kpsi_SF      += FsRT*KF_SF*z_sf*c_sf ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaHSF   += FsRT*KF_CaHSF*z_cahsf*c_cahsf ;
    Kpsi_CaH3SiO4 += FsRT*KF_CaH3SiO4*z_cah3sio4*c_cah3sio4 ;
    Kpsi_H3SiO4   += FsRT*KF_H3SiO4*z_h3sio4*c_h3sio4 ;
    Kpsi_H2SiO4   += FsRT*KF_H2SiO4*z_h2sio4*c_h2sio4 ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    
    Kpsi_K        += FsRT*KF_K*z_k*c_k ;
    Kpsi_Cl       += FsRT*KF_Cl*z_cl*c_cl ;
    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hsf*Kpsi_HSF + z_sf*Kpsi_SF + z_ca*Kpsi_Ca + z_cahsf*Kpsi_CaHSF + z_h3sio4*Kpsi_H3SiO4 \
                   + z_cah3sio4*Kpsi_CaH3SiO4 + z_caoh*Kpsi_CaOH + z_h2sio4*Kpsi_H2SiO4 + z_k*Kpsi_K + z_cl*Kpsi_Cl ;
    /*KD_CSF       += n_csf*k_g;*/
  }
  
  /* moyenne */
  for(i=0;i<NVE_TR;i++) va[i] *= 0.5 ;
}


void flux(Element_t *el,double **u)
/* Les flux (f) */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  double r_h[2],r_oh[2] ;
  double r_h2sf[2],r_hsf[2],r_sf[2] ;
  double r_h3sio4[2],r_h4sio4[2],r_h2sio4[2] ;
  double r_ca[2],r_caoh[2] ;
  double r_cahsf[2],r_casfaq[2],r_cah3sio4[2],r_cah2sio4[2] ;
  double r_k[2] ;
  double r_cl[2] ;
  double r_csf[2],r_p_h2s[2] ;
  double n_ch0=0;
  int    i ;
  double p_h2s      = Element_GetProperty(el)[pm("P_H2S")] ;

  for(i=0;i<2;i++) {
    /* molarities */
    double *x = ComputeAqueousVariables(el,u,i) ;
    double xx = Element_GetNodeCoordinate(el,i)[0] ;
    n_ch0 = N_CH0(i);
    
    r_oh[i]       = x[I_C_OH] ;
    r_h[i]        = x[I_C_H] ;
    r_h2sf[i]    = x[I_C_H2SF] ;
    r_hsf[i]     = x[I_C_HSF] ;
    r_sf[i]      = x[I_C_SF] ;
    r_ca[i]       = x[I_C_Ca] ;
    r_cahsf[i]   = x[I_C_CaHSF] ;
    r_h3sio4[i]   = x[I_C_H3SiO4] ;
    r_h4sio4[i]   = x[I_C_H4SiO4] ;
    r_cah3sio4[i] = x[I_C_CaH3SiO4] ;
    r_h2sio4[i]   = x[I_C_H2SiO4] ;
    r_cah2sio4[i] = x[I_C_CaH2SiO4] ;
    r_caoh[i]     = x[I_C_CaOH] ;
    r_casfaq[i]  = x[I_C_CaSFaq] ;
    r_k[i]        = x[I_C_K] ;
    r_cl[i]       = x[I_C_Cl] ;
    r_p_h2s[i]    = (xx==1)?1:0 ;
    
  }

  /* Gradients */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    double xx = Element_GetNodeCoordinate(el,i)[0] ;
    double grd_h        = (r_h[1]        - r_h[0]       )/dx ;
    double grd_oh       = (r_oh[1]       - r_oh[0]      )/dx ;
    double grd_h2sf    = (r_h2sf[1]    - r_h2sf[0]   )/dx ;
    double grd_hsf     = (r_hsf[1]     - r_hsf[0]    )/dx ;
    double grd_sf      = (r_sf[1]      - r_sf[0]     )/dx ;
    double grd_cahsf   = (r_cahsf[1]   - r_cahsf[0]  )/dx ;
    double grd_ca       = (r_ca[1]       - r_ca[0]      )/dx ;
    double grd_cah3sio4 = (r_cah3sio4[1] - r_cah3sio4[0])/dx ;
    double grd_h3sio4   = (r_h3sio4[1]   - r_h3sio4[0]  )/dx ;
    double grd_h4sio4   = (r_h4sio4[1]   - r_h4sio4[0]  )/dx ;
    double grd_h2sio4   = (r_h2sio4[1]   - r_h2sio4[0]  )/dx ;
    double grd_cah2sio4 = (r_cah2sio4[1] - r_cah2sio4[0])/dx ;
    double grd_casfaq  = (r_casfaq[1]  - r_casfaq[0] )/dx ;
    double grd_caoh     = (r_caoh[1]     - r_caoh[0]    )/dx ;
    double grd_k        = (r_k[1]        - r_k[0]       )/dx ;
    double grd_cl       = (r_cl[1]       - r_cl[0]      )/dx ;
        
    double grd_psi      = (PSI(1)        - PSI(0)       )/dx ;
    double grd_csf     = (r_csf[1]     - r_csf[0]    )/dx ;
    double grd_p_h2s     = (r_p_h2s[1]     - r_p_h2s[0]    );
    
    /*   producing H2SF */
    /*a_h2s = absorptionrate(p_h2s,el);*/
    /*double nk    = (p_h2s > 25)? 0.7:0.5;*/
    double a_h2s = KF_A_H2S*pow(p_h2s,0.55);
    
    /* Flux */
    
    double w_h2sf    =  - KF_H2SF*grd_h2sf;
     
    A_H2S = - a_h2s;
    
    double w_hsf     = - KF_HSF*grd_hsf          - Kpsi_HSF*grd_psi ;
    double w_sf      = - KF_SF*grd_sf            - Kpsi_SF*grd_psi  ;
    double w_cahsf   = - KF_CaHSF*grd_cahsf      - Kpsi_CaHSF*grd_psi ;
    double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
    double w_cah3sio4 = - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
    double w_h3sio4   = - KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi ;
    double w_h2sio4   = - KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi ;
    double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi ;
    double w_h4sio4   = - KF_H4SiO4*grd_h4sio4 ;
    double w_cah2sio4 = - KF_CaH2SiO4*grd_cah2sio4 ;
    double w_casfaq  = - KF_CaSFaq*grd_casfaq ;    
    double w_k        = - KF_K*grd_k                - Kpsi_K*grd_psi ;
    double w_cl       = - KF_Cl*grd_cl              - Kpsi_Cl*grd_psi ; 
    
    double w_q        = - z_h*KF_H*grd_h		      \
                        - z_oh*KF_OH*grd_oh		      \
                        - z_hsf*KF_HSF*grd_hsf             \
                        - z_sf*KF_SF*grd_sf		      \
                        - z_ca*KF_Ca*grd_ca		      \
                        - z_cahsf*KF_CaHSF*grd_cahsf	      \
                        - z_h3sio4*KF_H3SiO4*grd_h3sio4	      \
                        - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                        - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                        - z_caoh*KF_CaOH*grd_caoh \
                        - z_k*KF_K*grd_k \
                        - z_cl*KF_Cl*grd_cl \
                        - Kpsi_q*grd_psi ;
   
    /*A_H2S  = a_h2s ;*/
    W_S     = w_h2sf + w_hsf + w_sf + w_cahsf + w_casfaq ;/* + w_csf  ;*/
    W_Ca    = w_ca + w_cahsf + w_cah3sio4 + w_casfaq + w_caoh + w_cah2sio4; /*w_csf  ;*/
    W_Si    = w_h3sio4 + w_h4sio4 + w_cah3sio4 + w_cah2sio4 + w_h2sio4 ;
    W_q     = w_q ;
    W_K     = w_k ;
    W_Cl    = w_cl;
    
   } 
}


double concentration_oh(double c_h2sf,elem_t *el,double zn_ca_s,double c_k,double c_cl,double zn_si_s)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_h2sf, zn_si_s et zn_ca_s sont fixes */
  double c_h2sf_eq   = C_H2SF_eq ;
  double zc_h2sf     = c_h2sf/c_h2sf_eq ;
  /* les produits sont donc aussi fixes */
  double Q_CSF     = IAP_CSF(zc_h2sf,zn_ca_s) ;
  double Q_CH       = IAP_CH(zc_h2sf,zn_ca_s) ;
  double zq_ch      = Q_CH/K_CH ;
  double c_h4sio4   = IAP_SH(zq_ch,zn_si_s) ;

  /*
  rappel des expressions c_i = A_i*(c_h)**n   : n
     c_h        = K_h2o/c_oh                     : +1
     c_hsf     = K_hsf*c_h2sf/c_h             : -1
     c_sf      = K_sf*c_hsf/c_h               : -2
     c_ca       = Q_CSF/c_sf                     : +2
     c_cahsf   = K_cahsf*c_ca*c_hsf           : +1
     c_h4sio4   = IAP_SH                         :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)        : -1
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4       : +1
     c_casfaq  = K_casfaq*c_ca*c_sf ;         :  0      
     c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;       : -2       
     c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;     :  0      
     c_caoh     = K_caoh*c_ca*c_oh ;             : +1       
  */
  double A_hsf     = K_hsf*c_h2sf ;
  double A_sf      = K_sf*A_hsf ;
  double A_ca       = Q_CSF/A_sf ;
  double A_cahsf   = K_cahsf*A_ca*A_hsf ;
  double A_h3sio4   = c_h4sio4/K_h4sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;
  double A_h2sio4   = K_h2sio4*A_h3sio4*K_h2o ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahsf*A_cahsf + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_k*c_k + z_cl*c_cl;
  double d = z_oh*K_h2o + z_hsf*A_hsf + z_h3sio4*A_h3sio4 ;
  double e = z_sf*A_sf + z_h2sio4*A_h2sio4 ;

  /* a > 0 ; b > 0 ; c > 0 ; d < 0 ; e < 0 */
  double c_h = poly41(a,b,c,d,e) ;
 
  return(K_h2o/c_h) ;
}

double productionrate(double c_h2sf,elem_t *el, double c_oh,double j_h2sf)
{
	double ph         = 14 + log(c_oh)/log(10.) ;
	double c_h         = K_h2o/c_oh ;
	double c_h_i = 0;
	double c_h2sf_i = 0;
	
	/*if (ph>4.5) {return c_h2sf_i = (0.003/86400)*log(10)*c_h;}
	if ((ph>1.5)&&(ph<=4.5)) { return c_h2sf_i = (0.014/86400)*log(10)*c_h;}
	if ((ph>0.5)&&(ph<=1.5)) { return c_h2sf_i = (0.0006/86400)*log(10)*(0.0316/c_h);}
	else return c_h2sf_i = (0.007/86400)*log(10)/100;
    /*if (ph<0) {return c_h2sf_i = 5e-7;}*/
    return c_h2sf_i = (0.0006/86400)*log(10)*100;
	
	}
	
	double absorptionrate(double p_h2s,elem_t *el)
{
	double a_h2s = 0;
	return a_h2s = KF_A_H2S*pow(p_h2s,0.55);
	
	}

double poly4(double a,double b,double c,double d,double e)
/* on resout ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double x0 = pow(-d/a,1./3) ;/*ana_abcd(a,b,c,d);*/
  double x  = x0 ;
  int    i  = 0 ;
  /* 
     on neglige a et c
     solution de bx^3 + dx + e = 0 
     que l'on met sous la forme x^3 + px + q = 0
  double p = d/b,q = e/b ;
  double r = sqrt(q*q/4 + p*p*p/27) ;
  double uu = pow(-q/2. + r,1/3.), vv = pow(-q/2. - r,1/3.) ;
  double x = uu + vv ; *//* formule de Cardan */
  
  /* Newton */
  do {
    double x2 = x*x,x3 = x2*x,x4 = x2*x2 ;
    double f  = a*x4 + b*x3 + c*x2 + d*x + e ;
    double df = 4*a*x3 + 3*b*x2 + 2*c*x + d ;
    double dx = -f/df ;
    err = fabs(dx/x) ;
    x += dx ;
    if(i++ > 40) {
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      arret("poly4 : non convergence") ;
    }
  } while(err > tol) ;
  return(x) ;
}



double* ComputeAqueousVariables(Element_t *el,double **u,int n)
{
  double *x = var ;
  
  /*x[I_C_H2SF]  = (C_H2SF(n)>1.e-30)? C_H2SF(n):1.e-30 ;*/
  x[I_C_H2SF]  = C_H2SF(n) ;
  x[I_ZN_Ca_S ] = ZN_Ca_S(n) ;
  x[I_ZN_Si_S ] = ZN_Si_S(n) ;
  x[I_C_K    ]  = C_K(n) ;
  x[I_C_Cl   ]  = C_Cl(n) ;  
  
  /*printf("c_cl         = %e\n",x[I_C_Cl   ]) ;
  printf("c_k         = %e\n",C_K(n)) ;
  printf("c_oh         = %e\n",x[I_C_OH      ]) ;
  /*printf("n            = %e\n",n) ;*/
  
  ComputeAqueousParameters(el,x) ;
  return(x) ;
}


double* ComputeAqueousDerivatives(Element_t *el,double *x,double dxi,int i)
{
  double *dx = dvar ;
  int j ;
  
  dx[I_C_H2SF] = x[I_C_H2SF] ;
  dx[I_ZN_Ca_S ] = x[I_ZN_Ca_S ] ;
  dx[I_ZN_Si_S ] = x[I_ZN_Si_S ] ;
  dx[I_C_K ] = x[I_C_K] ;
  dx[I_C_Cl] = x[I_C_Cl] ;
  
  dx[i] += dxi ;
  
  ComputeAqueousParameters(el,dx) ;
  
  for(j = 0 ; j < NbOfAqueousVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
} 



void  ComputeAqueousParameters(Element_t *el,double *x)
{
  double c_h2sf    = x[I_C_H2SF] ;
  double zn_si_s    = x[I_ZN_Si_S] ;
  double zn_ca_s    = x[I_ZN_Ca_S] ;
  double c_k        = x[I_C_K] ;
  double c_cl       = x[I_C_Cl] ;      
    
  double zc_h2sf   = c_h2sf/c_h2sf_eq ;

  double Q_CSF     = IAP_CSF(zc_h2sf,zn_ca_s) ;
  double Q_CH       = IAP_CH(zc_h2sf,zn_ca_s) ;
  double zq_ch      = Q_CH/K_CH ;
  double c_h4sio4   = IAP_SH(zq_ch,zn_si_s) ;

  double c_oh       = concentration_oh(c_h2sf,el,zn_ca_s,c_k,c_cl,zn_si_s) ;
  double c_h        = K_h2o/c_oh ;
  double c_hsf     = K_hsf*c_h2sf/c_h ;
  double c_sf      = K_sf*c_hsf/c_h ;
  double c_ca       = Q_CSF/c_sf ;
  double c_cahsf   = K_cahsf*c_ca*c_hsf ;
  double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
  double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
  double c_casfaq  = K_casfaq*c_ca*c_sf ;
  double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
  double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
  double c_caoh     = K_caoh*c_ca*c_oh ;
  
  x[I_C_OH      ] = c_oh ;
  x[I_C_H       ] = c_h ;
  
  x[I_C_HSF    ] = c_hsf ;
  x[I_C_SF     ] = c_sf ;
  
  x[I_C_Ca      ] = c_ca ;
  
  x[I_C_H4SiO4  ] = c_h4sio4 ;
  x[I_C_H3SiO4  ] = c_h3sio4 ;
  x[I_C_H2SiO4  ] = c_h2sio4 ;
  
  x[I_C_CaH3SiO4] = c_cah3sio4 ;
  x[I_C_CaH2SiO4] = c_cah2sio4 ;
  
  x[I_C_CaHSF  ] = c_cahsf ;
  x[I_C_CaSFaq ] = c_casfaq ;
  
  x[I_C_CaOH    ] = c_caoh ;
  
  x[I_ZQ_CH   ] = zq_ch ;

  x[I_C_S_L     ] = c_h2sf + c_hsf + c_sf + c_cahsf + c_casfaq ;
  x[I_C_Ca_L    ] = c_ca + c_cahsf + c_cah3sio4 + c_cah2sio4 + c_casfaq + c_caoh ;
  x[I_C_Si_L    ] = c_h3sio4 + c_h4sio4 + c_cah3sio4 + c_h2sio4 + c_cah2sio4 ;
  x[I_C_K_L     ] = c_k ;
  x[I_C_Cl_L    ] = c_cl ;

  /*printf("c_cl1         = %e\n",c_cl) ;
  printf("c_oh         = %e\n",c_oh) ;*/
  
  /* charge density */
  x[I_N_Q     ]  = z_h*c_h + z_oh*c_oh \
                 + z_hsf*c_hsf + z_sf*c_sf \
                 + z_ca*c_ca + z_caoh*c_caoh \
                 + z_h3sio4*c_h3sio4 + z_h2sio4*c_h2sio4 \
                 + z_cah3sio4*c_cah3sio4 + z_cahsf*c_cahsf \
                 + z_k*c_k + z_cl*c_cl;
}


double* ComputeAqueousDerivatives1(Element_t *el,double *x,double dxi,int i)
{
  double c_h2sf    = x[I_C_H2SF] ;
  double zn_si_s    = x[I_ZN_Si_S] ;
  double zn_ca_s    = x[I_ZN_Ca_S] ;
  double c_k        = x[I_C_K] ;
  double c_cl       = x[I_C_Cl] ;
  double c_oh       = x[I_C_OH] ;
  double c_h        = x[I_C_H] ;
  double c_hsf     = x[I_C_HSF] ;
  double c_sf      = x[I_C_SF] ;
  double c_ca       = x[I_C_Ca] ;
  double c_h2sio4   = x[I_C_H2SiO4] ;
  double c_h3sio4   = x[I_C_H3SiO4] ;
  double c_h4sio4   = x[I_C_H4SiO4] ;
  double zq_ch      = x[I_ZQ_CH] ;
  double zc_h2sf   = c_h2sf/c_h2sf_eq ;
  double Q_CSF     = IAP_CSF(zc_h2sf,zn_ca_s) ;
  double Q_CH       = IAP_CH(zc_h2sf,zn_ca_s) ;

  double dc_oh       = 0 ;
  double dc_h        = 0 ;
  double dc_h2sf    = 0 ;
  double dc_hsf     = 0 ;
  double dc_sf      = 0 ;
  double dc_ca       = 0 ;
  double dc_h4sio4   = 0 ;
  double dc_h3sio4   = 0 ;
  double dc_h2sio4   = 0 ;
  double dc_cah3sio4 = 0 ;
  double dc_cah2sio4 = 0 ;
  double dc_casfaq  = 0 ;
  double dc_cahsf   = 0 ;
  double dc_caoh     = 0 ;
  double dc_k        = 0 ;
  double dc_cl       = 0 ;  
  double dzq_ch      = 0 ;
  
  double *dx = dvar ;
  int j ;
  
  for(j = 0 ; j < NbOfAqueousVariables ; j++) dx[j] = 0 ;
  
  if(i == I_C_H2SF) {
    double dc          = dxi ;
    double cc_h2sf    = c_h2sf + dc ;
    double cc_oh       = concentration_oh(cc_h2sf,el,zn_ca_s,c_k,c_cl,zn_si_s) ;
			                     
    double dQ_CSF     = (zc_h2sf < 1.) ? Q_CSF/c_h2sf : 0. ;
    double dQ_CH       = (zc_h2sf < 1.) ? 0. : -Q_CH/c_h2sf ;
    
    dc_h2sf    = 1 ;
    dc_oh       = (cc_oh - c_oh)/dc ;
    dc_h        = - c_h*dc_oh/c_oh ;
    dc_hsf     = c_hsf*(dc_h2sf/c_h2sf + dc_oh/c_oh) ;
    dc_sf      = c_sf*(dc_hsf/c_hsf + dc_oh/c_oh) ;
    
    dc_ca       = (dQ_CSF - c_ca*dc_sf)/c_sf ;
    dc_cahsf   = K_cahsf*(dc_ca*c_hsf + c_ca*dc_hsf) ;

    dzq_ch      = dQ_CH/K_CH ;
    dc_h4sio4   = DIAP_SHSDQ_CH(zq_ch,zn_si_s)*dzq_ch  ;
    dc_h3sio4   = c_h3sio4*(dc_h4sio4/c_h4sio4 - dc_h/c_h) ;
    dc_cah3sio4 = K_cah3sio4*(dc_ca*c_h3sio4 + c_ca*dc_h3sio4) ;
    dc_h2sio4   = K_h2sio4*(dc_h3sio4*c_oh	+ dc_oh*c_h3sio4) ;
    dc_cah2sio4 = K_cah2sio4*(dc_h2sio4*c_ca	+ dc_ca*c_h2sio4) ;
    dc_casfaq  = K_casfaq*(dc_sf*c_ca + dc_ca*c_sf) ;
    dc_caoh     = K_caoh*(dc_ca*c_oh + dc_oh*c_ca) ;
    dc_k        = 0 ;
    dc_cl        = 0 ;
    
  } else if(i == I_ZN_Ca_S) {
    double dz          = dxi ;
    double zzn_ca_s    = zn_ca_s + dz ;
    double cc_oh       = concentration_oh(c_h2sf,el,zzn_ca_s,c_k,c_cl,zn_si_s) ;

    double dQ_CSF     = (zn_ca_s < 0) ? Q_CSF : 0 ;
    double dQ_CH       = (zn_ca_s < 0) ? Q_CH : 0 ;
    
    dc_oh       = (cc_oh - c_oh)/dz ;
    dc_h        = - c_h*dc_oh/c_oh ;
    dc_hsf     = c_hsf*(dc_oh/c_oh) ;
    dc_sf      = c_sf*(dc_hsf/c_hsf + dc_oh/c_oh) ;
    
    dc_ca       = (dQ_CSF - c_ca*dc_sf)/c_sf ;
    dc_cahsf   = K_cahsf*(dc_ca*c_hsf + c_ca*dc_hsf) ;
    
    dzq_ch      = dQ_CH/K_CH ;
    dc_h4sio4   = DIAP_SHSDQ_CH(zq_ch,zn_si_s)*dzq_ch  ;
    dc_h3sio4   = c_h3sio4*(dc_h4sio4/c_h4sio4 - dc_h/c_h) ;
    dc_cah3sio4 = K_cah3sio4*(dc_ca*c_h3sio4 + c_ca*dc_h3sio4) ;
    dc_h2sio4   = K_h2sio4*(dc_h3sio4*c_oh + dc_oh*c_h3sio4) ;
    dc_cah2sio4 = K_cah2sio4*(dc_h2sio4*c_ca	+ dc_ca*c_h2sio4) ;
    dc_casfaq  = K_casfaq*(dc_sf*c_ca + dc_ca*c_sf) ;
    dc_caoh     = K_caoh*(dc_ca*c_oh + dc_oh*c_ca) ;	
    dc_k        = 0 ;
    dc_cl        = 0 ;
    dc_h2sf    = 0 ;
    
  } else if(i == I_ZN_Si_S) {
    double dz          = dxi ;
    double zzn_si_s    = zn_si_s + dz ;
    double cc_oh       = concentration_oh(c_h2sf,el,zn_ca_s,c_k,c_cl,zzn_si_s) ;
    
    dc_oh       = (cc_oh - c_oh)/dz ;
    dc_h        = - c_h*dc_oh/c_oh ;
    dc_hsf     = c_hsf*(dc_oh/c_oh) ;
    dc_sf      = c_sf*(dc_hsf/c_hsf + dc_oh/c_oh) ;
    dc_ca       = - c_ca*dc_sf/c_sf ;
    dc_cahsf   = K_cahsf*(dc_ca*c_hsf + c_ca*dc_hsf) ;
    dc_h4sio4   = DIAP_SHSDZN_Si_S(zq_ch,zn_si_s) ;
    dc_h3sio4   = c_h3sio4*(dc_h4sio4/c_h4sio4 - dc_h/c_h) ;
    dc_cah3sio4 = K_cah3sio4*(dc_ca*c_h3sio4 + c_ca*dc_h3sio4) ;
    dc_h2sio4   = K_h2sio4*(dc_h3sio4*c_oh	+ dc_oh*c_h3sio4) ;
    dc_cah2sio4 = K_cah2sio4*(dc_h2sio4*c_ca + dc_ca*c_h2sio4) ;
    dc_casfaq  = K_casfaq*(dc_sf*c_ca + dc_ca*c_sf) ;
    dc_caoh     = K_caoh*(dc_ca*c_oh + dc_oh*c_ca) ;
    dc_k        = 0 ;
    dc_cl        = 0 ;
    dc_h2sf    = 0 ;
    dzq_ch      = 0 ;
    
  } else if(i == I_C_K) {
    double dc          = dxi ;
    double cc_k        = c_k + dc ;
    double cc_oh       = concentration_oh(c_h2sf,el,zn_ca_s,cc_k,c_cl,zn_si_s) ;
    
    dc_oh       = (cc_oh - c_oh)/dc ;
    dc_h        = - c_h*dc_oh/c_oh ;
    dc_hsf     = c_hsf*(dc_oh/c_oh) ;
    dc_sf      = c_sf*(dc_hsf/c_hsf + dc_oh/c_oh) ; 
	                           
    dc_ca       = - c_ca*dc_sf/c_sf ;
    dc_cahsf   = K_cahsf*(dc_ca*c_hsf + c_ca*dc_hsf) ;

    dc_h3sio4   = - c_h3sio4*dc_h/c_h ;
    dc_cah3sio4 = K_cah3sio4*(dc_ca*c_h3sio4 + c_ca*dc_h3sio4) ;
    dc_h2sio4   = K_h2sio4*(dc_h3sio4*c_oh + dc_oh*c_h3sio4) ;
    dc_cah2sio4 = K_cah2sio4*(dc_h2sio4*c_ca + dc_ca*c_h2sio4) ;
    dc_casfaq  = K_casfaq*(dc_sf*c_ca + dc_ca*c_sf) ;
    dc_caoh     = K_caoh*(dc_ca*c_oh + dc_oh*c_ca) ;
    dc_h4sio4   = 0 ;
    dc_k        = 1 ;
    dc_cl       = 0 ;
    dc_h2sf    = 0 ;
    dzq_ch      = 0 ;
    
   } else if(i == I_C_Cl) {
    double dc          = dxi ;
    double cc_cl       = c_cl + dc ;
    double cc_oh       = concentration_oh(c_h2sf,el,zn_ca_s,c_k,cc_cl,zn_si_s) ;
    
    dc_oh       = (cc_oh - c_oh)/dc ;
    dc_h        = - c_h*dc_oh/c_oh ;
    dc_hsf     = c_hsf*(dc_oh/c_oh) ;
    dc_sf      = c_sf*(dc_hsf/c_hsf + dc_oh/c_oh) ; 
	                           
    dc_ca       = - c_ca*dc_sf/c_sf ;
    dc_cahsf   = K_cahsf*(dc_ca*c_hsf + c_ca*dc_hsf) ;

    dc_h3sio4   = - c_h3sio4*dc_h/c_h ;
    dc_cah3sio4 = K_cah3sio4*(dc_ca*c_h3sio4 + c_ca*dc_h3sio4) ;
    dc_h2sio4   = K_h2sio4*(dc_h3sio4*c_oh + dc_oh*c_h3sio4) ;
    dc_cah2sio4 = K_cah2sio4*(dc_h2sio4*c_ca + dc_ca*c_h2sio4) ;
    dc_casfaq  = K_casfaq*(dc_sf*c_ca + dc_ca*c_sf) ;
    dc_caoh     = K_caoh*(dc_ca*c_oh + dc_oh*c_ca) ;
    dc_h4sio4   = 0 ;
    dc_k        = 0 ;
    dc_cl       = 1 ;
    dc_h2sf    = 0 ;
    dzq_ch      = 0 ;   
  } else {
    arret("ComputeAqueousDerivatives1") ;
  }
  
    
  dx[I_C_OH      ] = dc_oh ;
  dx[I_C_H       ] = dc_h ;
  dx[I_C_H2SF   ] = dc_h2sf ;
  dx[I_C_HSF    ] = dc_hsf ;
  dx[I_C_SF     ] = dc_sf ;
  dx[I_C_Ca      ] = dc_ca ;
  dx[I_C_H4SiO4  ] = dc_h4sio4 ;
  dx[I_C_H3SiO4  ] = dc_h3sio4 ;
  dx[I_C_H2SiO4  ] = dc_h2sio4 ;
  dx[I_C_CaH3SiO4] = dc_cah3sio4 ;
  dx[I_C_CaH2SiO4] = dc_cah2sio4 ;
  dx[I_C_CaHSF  ] = dc_cahsf ;
  dx[I_C_CaSFaq ] = dc_casfaq ;
  dx[I_C_CaOH    ] = dc_caoh ;
  
  dx[I_ZQ_CH   ] = dzq_ch ;

  dx[I_C_S_L     ] = dc_h2sf + dc_hsf + dc_sf + dc_cahsf + dc_casfaq ;
  dx[I_C_Ca_L    ] = dc_ca + dc_cahsf + dc_cah3sio4 + dc_cah2sio4 + dc_casfaq + dc_caoh ;
  dx[I_C_Si_L    ] = dc_h3sio4 + dc_h4sio4 + dc_cah3sio4 + dc_h2sio4 + dc_cah2sio4 ;
  dx[I_C_K_L     ] = dc_k ;
  dx[I_C_Cl_L    ] = dc_cl ;  

  dx[I_N_Q     ] = z_h*dc_h + z_oh*dc_oh \
                 + z_hsf*dc_hsf + z_sf*dc_sf \
                 + z_ca*dc_ca + z_caoh*dc_caoh \
                 + z_h3sio4*dc_h3sio4 + z_h2sio4*dc_h2sio4 \
                 + z_cah3sio4*dc_cah3sio4 + z_cahsf*dc_cahsf \
                 + z_k*dc_k + z_cl*dc_cl;

  return(dx) ;
} 

double poly41(double a,double b,double c,double d,double e)
/* on resout ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double iks[4],zed[2];
  double x0 = quartic_equation(a,b,c,d,e,iks,zed); /*quartic_equation(a,b,c,d,e,iks,zed)*/
  double x  = MAX(x0,0.) ;
  int    i = 0 ;
  /* 
     on neglige a et c
     solution de bx^3 + dx + e = 0 
     que l'on met sous la forme x^3 + px + q = 0
  double p = d/b,q = e/b ;
  double r = sqrt(q*q/4 + p*p*p/27) ;
  double uu = pow(-q/2. + r,1/3.), vv = pow(-q/2. - r,1/3.) ;
  double x = uu + vv ; *//* formule de Cardan */
  
  /* Newton */
  do {
    double x2 = x*x,x3 = x2*x,x4 = x2*x2 ;
    double f  = a*x4 + b*x3 + c*x2 + d*x + e ;
    double df = 4*a*x3 + 3*b*x2 + 2*c*x + d ;
    double dx = -f/df ;
    err = fabs(dx/x) ;
    x += dx ;
    if(i++ > 40) {
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      printf("c  = %e\n",c) ;
      arret("poly_sansc : non convergence") ;
    }
  } while(err > tol) ;
  return(x) ;
}

double cubic_equation(double a,double b,double c,double d)
{ /* solution de ax^3 + bx^2+ cx + d = 0 */ 
  double x,y1,y2,y3,p,q,r,v,uu,vv,x1,x2,x3,k;  

  p = c/a-b*b/(3.0*a*a);
  q = d/a + 2.0*pow(b,3.0)/(27.0*pow(a,3.0)) - b*c/(3.0*a*a);
  v = q*q/4.0 + p*p*p/27.0 ;
  if (v>=0){
  
  r = sqrt(v) ;
  uu = pow(-q/2. + r,1/3.), vv = ((-q/2. - r)>0.) ? pow(-q/2. - r,1/3.) : -pow(q/2. + r,1/3.);
  
  y1 = uu + vv ;
  x1 = y1 - b/(3.0*a) ;
  x2 = 0. ;
  x3 = 0.;
  }  
  else {
  k = acos(-1.5*q*sqrt(-3.*p)/(p*p)); 
  y1 = 2.*sqrt(-3.*p)*cos(k/3.)/3.;
  y2 = -sqrt(-3.*p) * (cos(k/3.) + sqrt(3.)*sin(k/3.)) /3.;
  y3 = -sqrt(-3.*p) * (cos(k/3.) - sqrt(3.)*sin(k/3.)) /3.;
  x1 = y1 - b/(3.*a) ;
  x2 = y2 - b/(3.*a) ;
  x3 = y3 - b/(3.*a) ;
  }
  x = MAX(x1,x2) ;
  return(MAX(x,x3)) ;
}

void quadratic_equation(double a,double b,double c,double *zed)
{ /* solution de ax^2+ bx + c = 0 */ 
  double delta,temp1,temp2;  
  delta = b*b - 4.*a*c;
  if (delta>=0.){
  temp1 = -b/(2.*a);
  temp2 = sqrt(delta)/(2.*a);
  zed[0] = temp1 + temp2;
  zed[1] = temp1 - temp2; 
  } 
  else{
  zed[0] = 0.;
  zed[1] = 0.;
  }     
}

double quartic_equation(double a,double b,double c,double d,double e,double *iks,double *zed)
{ /* solution de ax^4 + bx^3+ cx^2 + dx + e = 0 */ 
  double y,a2,b2,c2,a4,b4,c4,a3,b3,c3,d3,a1,b1,c1,d1,e1,x1,x2,x3,x4,x5,x6,x;
  
  a1 = a/a;
  b1 = b/a;
  c1 = c/a;
  d1 = d/a;
  e1 = e/a; 
  if(b==0. && c==0. && d==0. && e==0.){
  iks[0]=0.;
  iks[1]=0.;
  iks[2]=0.;
  iks[3]=0.;
  } 
  else if(b==0. && d==0. && e==0.){
  quadratic_equation(a,0.,c,zed);
  iks[0] = zed[0];
  iks[1] = zed[1];
  iks[2] = 0.;
  iks[3] = 0.;
  }   
  else{
  a3 = 8.;
  b3 = -4.*c/a;
  c3 = 2.0*b*d/(a*a)- 8.*e/a;
  d3 = e*(4.*c/a-b*b/(a*a))/a - d*d/(a*a);
  y = cubic_equation(a3,b3,c3,d3);
  a2 = 1.;
  b2 = b1/2. - sqrt(8.*y + b1*b1 - 4.*c1)/2.;
  c2 = y - (b1*y - d1)/sqrt(8.*y + b1*b1 - 4.*c1);
  quadratic_equation(a2,b2,c2,zed);
  iks[0]=zed[0];
  iks[1]=zed[1];
  a4 = 1.;
  b4 = b1/2. + sqrt(8.*y + b1*b1 - 4.*c1)/2.;
  c4 = y + (b1*y - d1)/sqrt(8.*y + b1*b1 - 4.*c1);
  quadratic_equation(a4,b4,c4,zed);
  iks[2]=zed[0];
  iks[3]=zed[1];
  }
  x1 = (iks[0]>0. && iks[0]<1.)? iks[0] : 0. ;
  x2 = (iks[0]>0. && iks[0]<1.)? iks[1] : 0. ;
  x3 = (iks[0]>0. && iks[0]<1.)? iks[2] : 0. ;
  x4 = (iks[0]>0. && iks[0]<1.)? iks[3] : 0. ;
  x5 = MAX(x1,x2);
  x6 = MAX(x3,x4);
  x  = MAX(x5,x6);
  return(x);
  
}
