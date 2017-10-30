/*

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "FVM.h"

#define TITLE "sc-Carbonation of saturated CBM (05/2011)"
#define AUTHORS "Shen"

#include "PredefinedMethods.h"


/* Macros */
#define NEQ     (5)
#define NVE     (28)
#define NVI     (21)
#define NV0     (8)

#define E_C     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)
#define E_K     (4)

#define I_H2CO3   (0)
#define I_psi   (1)
#define I_Ca_S  (2)
#define I_Si_S  (3)
#define I_K     (4)

#define RHO     1
#define LOG_RHO 2
#define Ln10    2.302585093
#define U_H2CO3   LOG_RHO


#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWN_n(n,i)   (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])

/*
#define UNKNOWN(n,i)     (Element_GetValueOfNodalUnknown(el,n,i))
#define UNKNOWN_n(n,i)   (Element_GetValueOfPreviousNodalUnknown(el,n,i))
*/


#if (U_H2CO3 == LOG_RHO)
  #define C_H2CO3(n)   (exp(Ln10*UNKNOWN(n,I_H2CO3)))
  #define C_H2CO3n(n)  (exp(Ln10*UNKNOWN_n(n,I_H2CO3)))
#else
  #define C_H2CO3(n)   (UNKNOWN(n,I_H2CO3))
  #define C_H2CO3n(n)  (UNKNOWN_n(n,I_H2CO3))
#endif
#define ZN_Ca_S(n)   (UNKNOWN(n,I_Ca_S))
#define ZN_Si_S(n)   (UNKNOWN(n,I_Si_S))
#define PSI(n)       (UNKNOWN(n,I_psi))
#define C_K(n)       (UNKNOWN(n,I_K))

#define ZN_Ca_Sn(n)  (UNKNOWN_n(n,I_Ca_S))
#define ZN_Si_Sn(n)  (UNKNOWN_n(n,I_Si_S))
#define PSIn(n)      (UNKNOWN_n(n,I_psi))
#define C_Kn(n)      (UNKNOWN_n(n,I_K))

#define N_C(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define N_Si(n)    (f[(6+n)])
#define N_K(n)     (f[(8+n)])
#define W_C        (f[10])
#define W_q        (f[11])
#define W_Ca       (f[12])
#define W_Si       (f[13])
#define W_K        (f[14])
#define N_CH(n)    (f[(15+n)])
#define N_CC(n)    (f[(17+n)])
#define N_Si_S(n)  (f[(19+n)])

#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_CHn(n)   (f_n[(15+n)])
#define N_CCn(n)   (f_n[(17+n)])
#define N_Si_Sn(n) (f_n[(19+n)])


#define KF_OH       (va[(0)])
#define KF_H        (va[(1)])
#define KF_CO2      (va[(2)])
#define KF_H2CO3    (va[(3)])
#define KF_HCO3     (va[(4)])
#define KF_CO3      (va[(5)])
#define KF_Ca       (va[(6)])
#define KF_CaHCO3   (va[(7)])
#define KF_CaH3SiO4 (va[(8)])
#define KF_H3SiO4   (va[(9)])
#define KF_H4SiO4   (va[(10)])
#define KF_H2SiO4   (va[(11)])
#define KF_CaH2SiO4 (va[(12)])
#define KF_CaCO3aq  (va[(13)])
#define KF_CaOH     (va[(14)])
#define KF_K        (va[(15)])

#define Kpsi_OH       (va[(16)])
#define Kpsi_H        (va[(17)])
#define Kpsi_HCO3     (va[(18)])
#define Kpsi_CO3      (va[(19)])
#define Kpsi_Ca       (va[(20)])
#define Kpsi_CaHCO3   (va[(21)])
#define Kpsi_CaH3SiO4 (va[(22)])
#define Kpsi_H3SiO4   (va[(23)])
#define Kpsi_q        (va[(24)])
#define Kpsi_H2SiO4   (va[(25)])
#define Kpsi_CaOH     (va[(26)])
#define Kpsi_K        (va[(27)])

#define N_CH0(n)      (v0[(0+n)])
#define N_CC0(n)      (v0[(2+n)])
#define N_Si_S0(n)    (v0[(4+n)])
#define ZP_CH0(n)     (v0[(6+n)])

/*
  Solution aqueuse
*/

/* les valences */
#define z_ca       (2.)
#define z_h        (1.)
#define z_oh       (-1.)
#define z_hco3     (-1.)
#define z_co3      (-2.)
#define z_h3sio4   (-1.)
#define z_cahco3   (1.)
#define z_cah3sio4 (1.)
#define z_h2sio4   (-2.)
#define z_caoh     (1.)
#define z_k        (1.)

/* volumes molaires partiels des ions (dm3/mole) */
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_h2o      (18.e-3)
#define v_co2      (32.81e-3)
#define v_h2co3    (50.e-3)
#define v_hco3     (24.21.e-3)   /* d'apres Lothenbach */
#define v_co3      (-6.06e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)
#define v_caco3aq  (26.20e-3)    /* a modifier */
#define v_caoh     (26.20e-3)    /* a modifier */
#define v_k        (43.93e-3)    /* d'apres Antoine */
#define v_koh      (27.44e-3)    /* d'apres Antoine */

/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s  */
#define d_oh       (1.22e-7)    /* 1.22e-7 (radius = 1.75e-10 m) */
#define d_h        (9.310e-7)    /* 4.76e-8 (radius = 4.5e-10 m) */
#define d_co2      (1.91e-7)     /* 1.43e-7 (radius = 1.5e-10 m) */
#define d_h2co3    d_co2         /* (7.2e-8) */
#define d_hco3     (1.18e-7)    /* 1.07e-7 (radius = 2e-10 m) */
#define d_co3      (9.55e-8)    /* 9.53e-8 (radius = 2.25e-10 m) */
#define d_ca       (7.92e-8)
#define d_cahco3   (1.07e-7)    /* (radius = 2e-10 m) */
#define d_sioh4    (xxx)
#define d_h4sio4   (1.07e-7)     /* */
#define d_h3sio4   (1.07e-7)   /* (radius = 2e-10 m) */
#define d_h2sio4   (1.07e-7)    /*(radius = 2e-10 m) */
#define d_cah2sio4 (1.07e-7)    /* a modifier */
#define d_cah3sio4 (1.07e-7)     /*(radius = 2e-10 m) */
#define d_caco3aq  (1.43e-7)    /* (radius = 1.5e-10 m) */
#define d_caoh     (1.07e-7)    /* (radius = 2e-10 m) */
#define d_k        (1.43e-7)   /* (radius = 1.5e-10 m) */
#define d_koh      (1.43e-7)   /* (radius = 1.5e-10 m) */

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_h2o      (1.e-14)          /* autoprotolyse de l'eau */
#define K_henry    (1.238)           /* cste de Henry du CO2 / RT */

#define K_hco3     (4.25e-7)         /* H2CO3   = HCO3[-] + H[+] */
#define K_co3      (4.69e-11)        /* HCO3[-] = CO3[2-] + H[+] */

#define K_h2sio4   (4.68)            /* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define K_h4sio4   (6.45e9)          /* H3SiO4[-] + H[+] = H4SiO4 */
#define K_h3sio4   (1.55e-10)        /* H4SiO4    = H3SiO4[-] + H[+] */

#define K_cahco3   (1.276e+1)        /* Ca[2+] + HCO3[-]    = CaHCO3[+] */
#define K_caco3aq  (1.4e+3)          /* Ca[2+] + CO3[2-]    = CaCO3[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-] = CaOH[+] */

/*
  Solides
  CH  = Portlandite
  CC  = Calcite
  CSH = Hydrated Calcium Silicates
  SH  = Amorphous Silica
*/

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* CH  = Ca[2+] + 2OH[-] */
#define K_CC       (3.89e-9)         /* CC  = Ca[2+] + CO3[2-]  */ 
#define K_SH       (1.93642e-3)      /* SHt = H4SiO4 + (t-2)H2O */
/* volumes molaires solides (dm3/mole) */
#define V_CH       (33.e-3)      /* (33.e-3) */
#define V_CC       (37.e-3)      /* (37.e-3) */


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


#define C_H2CO3_eq   (K_h2o*K_h2o*K_CC/(K_hco3*K_co3*K_CH))

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */


/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,elem_t*,double,double,double) ;
static double poly4(double,double,double,double,double) ;
static void   transfert(Element_t*,double**,double*) ;
static void   flux(Element_t*,double**) ;

#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)
#define DNEGEXP(y) ((y < 0.) ? exp(y) : 0.)

/* Ion Activity Products */
#define IAP_CC(zc_h2co3,zn_ca_s)           (K_CC*NEGEXP(zn_ca_s)*MIN(zc_h2co3,1.))
#define IAP_CH(zc_h2co3,zn_ca_s)           (K_CH*NEGEXP(zn_ca_s)/MAX(zc_h2co3,1.))
#define IAP_SH(zp_ch,zn_si_s)            (K_SH*NEGEXP(zn_si_s)*Q_SH(zp_ch))
#define DIAP_SHSDQ_CH(zp_ch,zn_si_s)     (K_SH*NEGEXP(zn_si_s)*DQ_SH(zp_ch))
#define DIAP_SHSDZN_Si_S(zp_ch,zn_si_s)  (K_SH*DNEGEXP(zn_si_s)*Q_SH(zp_ch))


#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_C,"carbone") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_Si,"silicium") ;
  Model_CopyNameOfEquation(model,E_q,"charge") ;
  Model_CopyNameOfEquation(model,E_K,"potassium") ;

#if (U_H2CO3 == LOG_RHO)
  Model_CopyNameOfUnknown(model,I_H2CO3,"logc_h2co3") ;
#else
  Model_CopyNameOfUnknown(model,I_H2CO3,"c_h2co3") ;
#endif
  Model_CopyNameOfUnknown(model,I_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,I_psi,"psi") ;
  Model_CopyNameOfUnknown(model,I_Si_S,"z_si") ;
  Model_CopyNameOfUnknown(model,I_K,"c_k") ;
  
  return(0) ;
}


/* Parametres */
static double phi0,c_h2co3_eq,t_ch,t_cc ;
static double n_ca_ref,n_si_ref ;


int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"N_Si") == 0)     return (2) ;
  else if(strcmp(s,"C_H2CO3_eq") == 0) return (3) ;
  else if(strcmp(s,"T_CH") == 0)     return (4) ;
  else if(strcmp(s,"T_CC") == 0)     return (5) ;
  else if(strcmp(s,"courbes") == 0)  return (6) ;
  else return(-1) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  FILE *ficd = DataFile_GetFileStream(datafile) ;
  int  n_donnees = 7 ;

  {
    /* initialisation automatique */
    double t_ch        = 600. ;
    double t_cc        = 0. ;
    double n_ca_ref    = 1. ;
    double n_si_ref    = 1. ;

    Material_GetProperty(mat)[pm("N_CH")] = n_ca_ref ;
    Material_GetProperty(mat)[pm("N_Si")] = n_si_ref ;
    Material_GetProperty(mat)[pm("T_CH")] = t_ch ;
    Material_GetProperty(mat)[pm("T_CC")] = t_cc ;

    dmat(mat,ficd,pm,n_donnees) ;

    t_ch      = Material_GetProperty(mat)[pm("T_CH")] ;
    t_cc      = Material_GetProperty(mat)[pm("T_CC")] ;

    if(t_cc  == 0.) Material_GetProperty(mat)[pm("T_CC")]  = t_ch ;

    Material_GetProperty(mat)[pm("C_H2CO3_eq")] = C_H2CO3_eq ;
  }
  
  return(n_donnees) ;
}



int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;

  printf("Description:\n") ;

  printf("Carbonation of saturated CBM. We take into account the following properties:\n ") ;

  printf("\t - alcalis\n") ;
  printf("\t - general approach for CSH\n") ;
  printf("\t - no water release\n") ;
  
  printf("\n\n") ;
  printf("The set of 5 equations is:\n") ;
#if (U_H2CO3 == LOG_RHO)
  printf("\t- mass balance of C      (logc_h2co3)\n") ;
#else
  printf("\t- mass balance of C      (c_h2co3)\n") ;
#endif
  printf("\t- charge balance         (psi)\n") ;
  printf("\t- mass balance of Ca     (zn_ca_s)\n") ;
  printf("\t- mass balance of Si     (zn_si_s)\n") ;
  printf("\t- mass balance of K      (c_k)\n") ;

  printf("\n\
ATTENTION to units : \n\
\t length : dm !\n\
\t time   : s !\n") ;

  printf("Example of input data\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1       # Contenu en CH (moles/L)\n") ;
  fprintf(ficd,"N_Si  = 2.4       # contenu en Si solide (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4       # contenu en K (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5      # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"T_CC  = 1.e5      # Cinetique de dissolution de CC (s)\n") ;
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
  int    i ;

  FVM_ComputeSurfaceLoadResidu(el,cg,t,dt,r) ;
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r[i] ;
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
  c_h2co3_eq  = Element_GetProperty(el)[pm("C_H2CO3_eq")] ;
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
    double c_h2co3      = (C_H2CO3(i) > 0.) ? C_H2CO3(i) : c_h2co3_eq ;
    double zc_h2co3     = c_h2co3/c_h2co3_eq ;
    double zn_si_s    = ZN_Si_S(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double c_k        = C_K(i) ;
    
    double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;
    
    double c_h4sio4   = P_Sam ;
    double c_oh       = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ; 
    
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    
    /* equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_cc_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    
    /* solid contents : CH, CC, CSH */
    double n_ch       = (zc_h2co3 <= 1) ? n_ch_eq  : 0 ;
    double n_cc       = (zc_h2co3 >  1) ? n_cc_eq  : 0 ;
    double x_csh      = X_CSH(zp_ch) ;
    double n_si_s     = n_si_ref*MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;

    /* porosity */
    double phi = phi0 ;

    /* molar contents */
    double n_h2co3    = phi*c_h2co3 ;
    double n_hco3     = phi*c_hco3 ;
    double n_co3      = phi*c_co3 ;
    double n_ca       = phi*c_ca ;
    double n_cahco3   = phi*c_cahco3 ;
    double n_h3sio4   = phi*c_h3sio4 ;
    double n_h4sio4   = phi*c_h4sio4 ;
    double n_cah3sio4 = phi*c_cah3sio4 ;
    double n_h2sio4   = phi*c_h2sio4 ;
    double n_cah2sio4 = phi*c_cah2sio4 ;
    double n_caco3aq  = phi*c_caco3aq ;
    double n_caoh     = phi*c_caoh ;
    double n_k        = phi*c_k ;
    
    N_C(i)  = n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc + n_caco3aq ;
    N_Ca(i) = n_ca + n_cahco3 + n_cah3sio4 + n_cah2sio4 + n_caco3aq + n_caoh + n_ca_s ;
    N_Si(i) = n_h3sio4 + n_h4sio4 + n_cah3sio4 + n_h2sio4 + n_cah2sio4 + n_si_s ;
    N_K(i)  = n_k ;

    /* charge density */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 \
    + z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k ;

#if (U_H2CO3 == LOG_RHO)
    UNKNOWN(i,I_H2CO3) = log(c_h2co3)/Ln10 ;
#else
    C_H2CO3(i)   = c_h2co3 ;
#endif
    ZN_Si_S(i) = zn_si_s ;
    ZN_Ca_S(i) = zn_ca_s ;
    C_K(i)     = c_k;

    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_Si_S(i)  = n_si_s ;

    N_CH0(i)   = n_ch ;
    N_CC0(i)   = n_cc ;
    N_Si_S0(i) = n_si_s ;
    ZP_CH0(i)  = zp_ch ;
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
  c_h2co3_eq  = Element_GetProperty(el)[pm("C_H2CO3_eq")] ;
  t_ch      = Element_GetProperty(el)[pm("T_CH")] ;
  t_cc      = Element_GetProperty(el)[pm("T_CC")] ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
    double c_h2co3      = C_H2CO3(i) ;
    double zn_si_s    = ZN_Si_S(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2co3     = c_h2co3/c_h2co3_eq ;
    double c_k        = C_K(i) ;

    double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;

    double c_h4sio4   = P_Sam ;
    double c_oh       = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* kinetics */
    double n_chn      = N_CHn(i) ;
    double n_ccn      = N_CCn(i) ;
    double n_ch_ci    = n_chn*pow(zc_h2co3,-dt/t_ch) ;  /* if zc_h2co3 > 1 */
    double n_cc_ci    = n_ccn*pow(zc_h2co3,dt/t_cc) ;   /* if zc_h2co3 < 1 */

    /* equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_cc_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    
    /* solid contents : CH, CC, CSH */
    double n_ch       = (zc_h2co3 <= 1) ? n_ch_eq  : n_ch_ci ;
    double n_cc       = (zc_h2co3 >  1) ? n_cc_eq  : n_cc_ci ;
    double x_csh      = X_CSH(zp_ch) ;
    double n_si_s     = n_si_ref*MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;

    /* porosity */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double n_si_s0    = N_Si_S0(i) ;
    double zp_ch0     = ZP_CH0(i) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double v_csh      = V_CSH(zp_ch) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) \
                      + v_csh0*n_si_s0 - v_csh*n_si_s ;

    /* molar contents */
    double n_h2co3    = phi*c_h2co3 ;
    double n_hco3     = phi*c_hco3 ;
    double n_co3      = phi*c_co3 ;
    double n_ca       = phi*c_ca ;
    double n_cahco3   = phi*c_cahco3 ;
    double n_h3sio4   = phi*c_h3sio4 ;
    double n_h4sio4   = phi*c_h4sio4 ;
    double n_cah3sio4 = phi*c_cah3sio4 ;
    double n_h2sio4   = phi*c_h2sio4 ;
    double n_cah2sio4 = phi*c_cah2sio4 ;
    double n_caco3aq  = phi*c_caco3aq ;
    double n_caoh     = phi*c_caoh ;
    double n_k        = phi*c_k ;

    N_C(i)  = n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc + n_caco3aq ;
    N_Ca(i) = n_ca + n_cahco3 + n_cah3sio4 + n_cah2sio4 + n_caco3aq + n_caoh + n_ca_s ;
    N_Si(i) = n_h3sio4 + n_h4sio4 + n_cah3sio4 + n_h2sio4 + n_cah2sio4 + n_si_s ;
    N_K(i)  = n_k ;

    /* charge density */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 + \
    z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k ;

    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_Si_S(i)  = n_si_s ;

    if(c_h2co3 < 0. || n_ca_s < 0. || n_si_s < 0.) {
      double x = Element_GetNodeCoordinate(el,i)[0] ;
      printf("x         = %e\n",x) ;
      printf("c_h2co3     = %e\n",c_h2co3) ;
      printf("n_cc      = %e\n",n_cc) ;
      printf("n_ca_s    = %e\n",n_ca_s) ;
      printf("n_si_s    = %e\n",n_si_s) ;
      printf("zn_si_s   = %e\n",zn_si_s) ;
      printf("zn_ca_s   = %e\n",zn_ca_s) ;
      printf("c_h3sio4  = %e\n",c_h3sio4) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
    if(phi < 0.) {
      printf("phi = %e\n",phi) ;
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
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double c[2] ;

  double Dc_hSDc_h2co3[2]        ;
  double Dc_ohSDc_h2co3[2]       ;
  double Dc_hco3SDc_h2co3[2]     ;
  double Dc_co3SDc_h2co3[2]      ;
  double Dc_caSDc_h2co3[2]       ;
  double Dc_cahco3SDc_h2co3[2]   ;
  double Dc_h3sio4SDc_h2co3[2]   ;
  double Dc_h4sio4SDc_h2co3[2]   ;
  double Dc_cah3sio4SDc_h2co3[2] ;
  double Dc_h2sio4SDc_h2co3[2]   ;
  double Dc_cah2sio4SDc_h2co3[2] ;
  double Dc_caco3aqSDc_h2co3[2]  ;  
  double Dc_caohSDc_h2co3[2]     ;
  
  double Dc_hSDzn_si_s[2]        ;
  double Dc_ohSDzn_si_s[2]       ;
  double Dc_hco3SDzn_si_s[2]     ;
  double Dc_co3SDzn_si_s[2]      ;
  double Dc_caSDzn_si_s[2]       ;
  double Dc_cahco3SDzn_si_s[2]   ;
  double Dc_h3sio4SDzn_si_s[2]   ;
  double Dc_h4sio4SDzn_si_s[2]   ;
  double Dc_cah3sio4SDzn_si_s[2] ;
  double Dc_h2sio4SDzn_si_s[2]   ;
  double Dc_cah2sio4SDzn_si_s[2] ;
  double Dc_caco3aqSDzn_si_s[2]  ;
  double Dc_caohSDzn_si_s[2]     ;

  double Dc_hSDzn_ca_s[2]        ;
  double Dc_ohSDzn_ca_s[2]       ;
  double Dc_hco3SDzn_ca_s[2]     ;
  double Dc_co3SDzn_ca_s[2]      ;
  double Dc_caSDzn_ca_s[2]       ;
  double Dc_cahco3SDzn_ca_s[2]   ;
  double Dc_h3sio4SDzn_ca_s[2]   ;
  double Dc_h4sio4SDzn_ca_s[2]   ;
  double Dc_cah3sio4SDzn_ca_s[2] ;
  double Dc_h2sio4SDzn_ca_s[2]   ;
  double Dc_cah2sio4SDzn_ca_s[2] ;
  double Dc_caco3aqSDzn_ca_s[2]  ;  
  double Dc_caohSDzn_ca_s[2]     ;
  
  double Dc_hSDc_k[2]         ;
  double Dc_ohSDc_k[2]        ;
  double Dc_hco3SDc_k[2]      ;
  double Dc_co3SDc_k[2]       ;
  double Dc_caSDc_k[2]        ;
  double Dc_cahco3SDc_k[2]    ;
  double Dc_h3sio4SDc_k[2]    ;
  double Dc_cah3sio4SDc_k[2]  ;
  double Dc_h2sio4SDc_k[2]    ;
  double Dc_cah2sio4SDc_k[2]  ;
  double Dc_caco3aqSDc_k[2]   ;  
  double Dc_caohSDc_k[2]      ;
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
  c_h2co3_eq  = Element_GetProperty(el)[pm("C_H2CO3_eq")] ;
  t_ch      = Element_GetProperty(el)[pm("T_CH")] ;
  t_cc      = Element_GetProperty(el)[pm("T_CC")] ;


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
    /* molarites */
    double c_h2co3      = C_H2CO3(i) ;
    double zn_si_s    = ZN_Si_S(i);
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2co3     = c_h2co3/c_h2co3_eq ;
    double c_k        = C_K(i) ;

    double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;
 
    double c_h4sio4   = P_Sam ;
    double c_oh       = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* solid contents : CH, CC, CSH */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_si_s      = N_Si_S(i) ;
    double x_csh      = X_CSH(zp_ch) ;

    /* kinetics */
    double n_ccn      = N_CCn(i) ;
    double n_chn      = N_CHn(i) ;
    double n_ch_ci    = n_chn*pow(zc_h2co3,-dt/t_ch) ;
    double n_cc_ci    = n_ccn*pow(zc_h2co3,dt/t_cc) ;

    /* porosity */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double n_si_s0    = N_Si_S0(i) ;
    double zp_ch0     = ZP_CH0(i) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double v_csh      = V_CSH(zp_ch) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) \
                      + v_csh0*n_si_s0 - v_csh*n_si_s ;

    /* derivatives ... */
    /* ... with respect to c_h2co3 */
    double dc_h2co3             = 1.e-4*c_h2co3*((c_h2co3 > C_H2CO3n(i)) ? 1 : -1) ;
    double c_h2co32             = c_h2co3 + dc_h2co3 ;
    double c_oh2              = concentration_oh(c_h2co32,el,zn_ca_s,c_k,zn_si_s) ;   
    double dc_ohsdc_h2co3       = (c_oh2 - c_oh)/dc_h2co3 ;
    double dc_hsdc_h2co3        = - c_h*dc_ohsdc_h2co3/c_oh ;
    double dc_hco3sdc_h2co3     = c_hco3*(1/c_h2co3	\
                              + dc_ohsdc_h2co3/c_oh) ;
    double dc_co3sdc_h2co3      = c_co3*(dc_hco3sdc_h2co3/c_hco3 \
                              + dc_ohsdc_h2co3/c_oh) ;
			                     
    double dP_CCsdc_h2co3       = (zc_h2co3 < 1.) ? P_CC/c_h2co3 : 0. ;
    double dP_CHsdc_h2co3       = (zc_h2co3 < 1.) ? 0. : -P_CH/c_h2co3 ;
    double dP_Samsdc_h2co3      = DIAP_SHSDQ_CH(zp_ch,zn_si_s)*dP_CHsdc_h2co3/K_CH  ;
    
    double dc_casdc_h2co3       = (dP_CCsdc_h2co3 - c_ca*dc_co3sdc_h2co3)/c_co3 ;
    double dc_cahco3sdc_h2co3   = K_cahco3*(dc_casdc_h2co3*c_hco3 \
                              + c_ca*dc_hco3sdc_h2co3) ;

    double dx_cshsdc_h2co3      = DX_CSH(zp_ch)*dP_CHsdc_h2co3/K_CH ;

    double dc_h4sio4sdc_h2co3   = dP_Samsdc_h2co3  ;
    double dc_h3sio4sdc_h2co3   = c_h3sio4*(dc_h4sio4sdc_h2co3/c_h4sio4 \
                              - dc_hsdc_h2co3/c_h) ;
    double dc_cah3sio4sdc_h2co3 = K_cah3sio4*(dc_casdc_h2co3*c_h3sio4 \
                              + c_ca*dc_h3sio4sdc_h2co3) ;
    double dc_h2sio4sdc_h2co3   = K_h2sio4*(dc_h3sio4sdc_h2co3*c_oh	\
                              + dc_ohsdc_h2co3*c_h3sio4) ;
    double dc_cah2sio4sdc_h2co3 = K_cah2sio4*(dc_h2sio4sdc_h2co3*c_ca	\
                              + dc_casdc_h2co3*c_h2sio4) ;
    double dc_caco3aqsdc_h2co3  = K_caco3aq*(dc_co3sdc_h2co3*c_ca \
                              + dc_casdc_h2co3*c_co3) ;
    double dc_caohsdc_h2co3     = K_caoh*(dc_casdc_h2co3*c_oh \
                              + dc_ohsdc_h2co3*c_ca) ;
		                    
    double dn_ch_cisdc_h2co3    = n_ch_ci*(-dt/t_ch)/c_h2co3 ;
    double dn_cc_cisdc_h2co3    = n_cc_ci*(dt/t_cc)/c_h2co3 ;

    double dn_ch_eqsdc_h2co3    = 0 ;
    double dn_cc_eqsdc_h2co3    = 0 ;


    double dn_chsdc_h2co3       = (zc_h2co3 <= 1) ? dn_ch_eqsdc_h2co3  : dn_ch_cisdc_h2co3 ;
    double dn_ccsdc_h2co3       = (zc_h2co3 >  1) ? dn_cc_eqsdc_h2co3  : dn_cc_cisdc_h2co3 ;
    
    double dn_ca_ssdc_h2co3     = dn_chsdc_h2co3 + dn_ccsdc_h2co3 \
                              + dx_cshsdc_h2co3*n_si_s ;
    double dn_si_ssdc_h2co3     = 0. ;

    double dv_cshsdc_h2co3      = DV_CSH(zp_ch)*dP_CHsdc_h2co3/K_CH ;
    double dphisdc_h2co3        = - V_CH*dn_chsdc_h2co3 - V_CC*dn_ccsdc_h2co3 \
                              - dv_cshsdc_h2co3*n_si_s ;

    double dn_h2co3sdc_h2co3    = phi + dphisdc_h2co3*c_h2co3 ;
    double dn_hco3sdc_h2co3     = phi*dc_hco3sdc_h2co3 + dphisdc_h2co3*c_hco3 ;
    double dn_co3sdc_h2co3      = phi*dc_co3sdc_h2co3 + dphisdc_h2co3*c_co3 ;
    double dn_casdc_h2co3       = phi*dc_casdc_h2co3 + dphisdc_h2co3*c_ca ;
    double dn_cahco3sdc_h2co3   = phi*dc_cahco3sdc_h2co3 + dphisdc_h2co3*c_cahco3 ;
    double dn_cah3sio4sdc_h2co3 = phi*dc_cah3sio4sdc_h2co3 + dphisdc_h2co3*c_cah3sio4 ;
    double dn_h3sio4sdc_h2co3   = phi*dc_h3sio4sdc_h2co3 + dphisdc_h2co3*c_h3sio4 ;
    double dn_h4sio4sdc_h2co3   = phi*dc_h4sio4sdc_h2co3 + dphisdc_h2co3*c_h4sio4 ;
    double dn_h2sio4sdc_h2co3   = phi*dc_h2sio4sdc_h2co3 + dphisdc_h2co3*c_h2sio4 ;
    double dn_cah2sio4sdc_h2co3 = phi*dc_cah2sio4sdc_h2co3 + dphisdc_h2co3*c_cah2sio4 ;
    double dn_caco3aqsdc_h2co3  = phi*dc_caco3aqsdc_h2co3 + dphisdc_h2co3*c_caco3aq ;
    double dn_caohsdc_h2co3     = phi*dc_caohsdc_h2co3 + dphisdc_h2co3*c_caoh ;
    double dn_ksdc_h2co3        = dphisdc_h2co3*c_k ;
			                    
    /* with respect to zn_si_s */
    double dzn_si_s            = ((zn_si_s > ZN_Si_Sn(i)) ? 1 : -1)*1.e-4; 
    double zn_si_s2            = zn_si_s + dzn_si_s ;
    
    double c_oh1               = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s2) ;
    double dc_ohsdzn_si_s      = (c_oh1 - c_oh)/dzn_si_s ;
    double dc_hsdzn_si_s       = - c_h*dc_ohsdzn_si_s/c_oh ;
    double dc_hco3sdzn_si_s    = c_hco3*(dc_ohsdzn_si_s/c_oh) ;
    double dc_co3sdzn_si_s     = c_co3*(dc_hco3sdzn_si_s/c_hco3 \
                               + dc_ohsdzn_si_s/c_oh) ;
                
    double dP_Samsdzn_si_s     = DIAP_SHSDZN_Si_S(zp_ch,zn_si_s) ;
    
    double dc_casdzn_si_s      = - c_ca*dc_co3sdzn_si_s/c_co3 ;
    double dc_cahco3sdzn_si_s  = K_cahco3*(dc_casdzn_si_s*c_hco3 \
                               + c_ca*dc_hco3sdzn_si_s) ;

    double dc_h4sio4sdzn_si_s  = dP_Samsdzn_si_s ;
    double dc_h3sio4sdzn_si_s  = c_h3sio4*(dc_h4sio4sdzn_si_s/c_h4sio4 \
                               - dc_hsdzn_si_s/c_h) ;
    double dc_cah3sio4sdzn_si_s = K_cah3sio4*(dc_casdzn_si_s*c_h3sio4 \
                               + c_ca*dc_h3sio4sdzn_si_s) ;
    double dc_h2sio4sdzn_si_s  = K_h2sio4*(dc_h3sio4sdzn_si_s*c_oh	\
                               + dc_ohsdzn_si_s*c_h3sio4) ;
    double dc_cah2sio4sdzn_si_s = K_cah2sio4*(dc_h2sio4sdzn_si_s*c_ca	\
                               + dc_casdzn_si_s*c_h2sio4) ;
    double dc_caco3aqsdzn_si_s = K_caco3aq*(dc_co3sdzn_si_s*c_ca \
                               + dc_casdzn_si_s*c_co3) ;
    double dc_caohsdzn_si_s    = K_caoh*(dc_casdzn_si_s*c_oh \
                               + dc_ohsdzn_si_s*c_ca) ;				      
 
    double dn_si_ssdzn_si_s    = (zn_si_s > 0) ? n_si_ref : 0. ;
    double dn_ca_ssdzn_si_s    = x_csh*dn_si_ssdzn_si_s ;
    
    double dphisdzn_si_s       = -v_csh*dn_si_ssdzn_si_s ;

    double dn_h2co3sdzn_si_s   = dphisdzn_si_s*c_h2co3 ;
    double dn_hco3sdzn_si_s    = phi*dc_hco3sdzn_si_s + dphisdzn_si_s*c_hco3 ;
    double dn_co3sdzn_si_s     = phi*dc_co3sdzn_si_s + dphisdzn_si_s*c_co3 ;
    double dn_casdzn_si_s      = phi*dc_casdzn_si_s + dphisdzn_si_s*c_ca ;
    double dn_cahco3sdzn_si_s  = phi*dc_cahco3sdzn_si_s + dphisdzn_si_s*c_cahco3 ;
    double dn_cah3sio4sdzn_si_s = phi*dc_cah3sio4sdzn_si_s + dphisdzn_si_s*c_cah3sio4 ;
    double dn_h3sio4sdzn_si_s  = phi*dc_h3sio4sdzn_si_s + dphisdzn_si_s*c_h3sio4 ;
    double dn_h4sio4sdzn_si_s  = phi*dc_h4sio4sdzn_si_s + dphisdzn_si_s*c_h4sio4 ;
    double dn_h2sio4sdzn_si_s  = phi*dc_h2sio4sdzn_si_s + dphisdzn_si_s*c_h2sio4 ;
    double dn_cah2sio4sdzn_si_s = phi*dc_cah2sio4sdzn_si_s + dphisdzn_si_s*c_cah2sio4 ;
    double dn_caco3aqsdzn_si_s = phi*dc_caco3aqsdzn_si_s + dphisdzn_si_s*c_caco3aq ;
    double dn_caohsdzn_si_s    = phi*dc_caohsdzn_si_s + dphisdzn_si_s*c_caoh ;
    double dn_ksdzn_si_s       = dphisdzn_si_s*c_k ;
    
    /* with respect to zn_ca_s */
    /* double dzn_ca_s             = ((zn_ca_s > 0.) ? 1 : -1)*1.e-2 ; */
    double dzn_ca_s            = 1.e-6*((zn_ca_s > ZN_Ca_Sn(i)) ? 1 : -1) ;
    double zn_ca_s2            = zn_ca_s + dzn_ca_s ;
    double c_oh3               = concentration_oh(c_h2co3,el,zn_ca_s2,c_k,zn_si_s) ;
    double dc_ohsdzn_ca_s      = (c_oh3 - c_oh)/dzn_ca_s ;
    double dc_hsdzn_ca_s       = - c_h*dc_ohsdzn_ca_s/c_oh ;
    double dc_hco3sdzn_ca_s    = c_hco3*(dc_ohsdzn_ca_s/c_oh) ;
    double dc_co3sdzn_ca_s     = c_co3*(dc_hco3sdzn_ca_s/c_hco3 \
                               + dc_ohsdzn_ca_s/c_oh) ;

    double dP_CCsdzn_ca_s      = (zn_ca_s < 0) ? P_CC : 0 ;
    double dP_CHsdzn_ca_s      = (zn_ca_s < 0) ? P_CH : 0 ;
    double dP_Samsdzn_ca_s     = DIAP_SHSDQ_CH(zp_ch,zn_si_s)*dP_CHsdzn_ca_s/K_CH ;
    
    double dc_casdzn_ca_s      = (dP_CCsdzn_ca_s - c_ca*dc_co3sdzn_ca_s)/c_co3 ;
    double dc_cahco3sdzn_ca_s  = K_cahco3*(dc_casdzn_ca_s*c_hco3 \
                               + c_ca*dc_hco3sdzn_ca_s) ;

    double dx_cshsdzn_ca_s     = DX_CSH(zp_ch)*dP_CHsdzn_ca_s/K_CH ;
    double dv_cshsdzn_ca_s     = DV_CSH(zp_ch)*dP_CHsdzn_ca_s/K_CH ;
    
    double dc_h4sio4sdzn_ca_s  = dP_Samsdzn_ca_s ;
    double dc_h3sio4sdzn_ca_s  = c_h3sio4*(dc_h4sio4sdzn_ca_s/c_h4sio4 \
                               - dc_hsdzn_ca_s/c_h) ;
    double dc_cah3sio4sdzn_ca_s = K_cah3sio4*(dc_casdzn_ca_s*c_h3sio4 \
                               + c_ca*dc_h3sio4sdzn_ca_s) ;
    double dc_h2sio4sdzn_ca_s  = K_h2sio4*(dc_h3sio4sdzn_ca_s*c_oh	\
                               + dc_ohsdzn_ca_s*c_h3sio4) ;
    double dc_cah2sio4sdzn_ca_s = K_cah2sio4*(dc_h2sio4sdzn_ca_s*c_ca	\
                               + dc_casdzn_ca_s*c_h2sio4) ;
    double dc_caco3aqsdzn_ca_s = K_caco3aq*(dc_co3sdzn_ca_s*c_ca \
                               + dc_casdzn_ca_s*c_co3) ;
    double dc_caohsdzn_ca_s    = K_caoh*(dc_casdzn_ca_s*c_oh \
                               + dc_ohsdzn_ca_s*c_ca) ;					       

    double dn_ch_cisdzn_ca_s   = 0 ;
    double dn_cc_cisdzn_ca_s   = 0 ;


    double dn_ch_eqsdzn_ca_s   = (zn_ca_s > 0) ? n_ca_ref : 0 ;
    double dn_cc_eqsdzn_ca_s   = (zn_ca_s > 0) ? n_ca_ref : 0 ;
     
    double dn_chsdzn_ca_s      = (zc_h2co3 <= 1) ? dn_ch_eqsdzn_ca_s  : dn_ch_cisdzn_ca_s ;
    double dn_ccsdzn_ca_s      = (zc_h2co3 >  1) ? dn_cc_eqsdzn_ca_s  : dn_cc_cisdzn_ca_s ;
    
    double dn_si_ssdzn_ca_s    = 0. ;
    double dn_ca_ssdzn_ca_s    = dn_chsdzn_ca_s + dn_ccsdzn_ca_s \
                               + dx_cshsdzn_ca_s*n_si_s ;

    double dphisdzn_ca_s       = - V_CH*dn_chsdzn_ca_s \
                               - V_CC*dn_ccsdzn_ca_s \
                               - dv_cshsdzn_ca_s*n_si_s ;
              
    double dn_h2co3sdzn_ca_s   = dphisdzn_ca_s*c_h2co3 ;
    double dn_hco3sdzn_ca_s    = phi*dc_hco3sdzn_ca_s + dphisdzn_ca_s*c_hco3 ;
    double dn_co3sdzn_ca_s     = phi*dc_co3sdzn_ca_s + dphisdzn_ca_s*c_co3 ;
    double dn_casdzn_ca_s      = phi*dc_casdzn_ca_s + dphisdzn_ca_s*c_ca ;
    double dn_cahco3sdzn_ca_s  = phi*dc_cahco3sdzn_ca_s + dphisdzn_ca_s*c_cahco3 ;
    double dn_cah3sio4sdzn_ca_s = phi*dc_cah3sio4sdzn_ca_s + dphisdzn_ca_s*c_cah3sio4 ;
    double dn_h3sio4sdzn_ca_s  = phi*dc_h3sio4sdzn_ca_s + dphisdzn_ca_s*c_h3sio4 ;
    double dn_h4sio4sdzn_ca_s  = phi*dc_h4sio4sdzn_ca_s + dphisdzn_ca_s*c_h4sio4 ;
    double dn_h2sio4sdzn_ca_s  = phi*dc_h2sio4sdzn_ca_s + dphisdzn_ca_s*c_h2sio4 ;
    double dn_cah2sio4sdzn_ca_s = phi*dc_cah2sio4sdzn_ca_s + dphisdzn_ca_s*c_cah2sio4 ;
    double dn_caco3aqsdzn_ca_s = phi*dc_caco3aqsdzn_ca_s + dphisdzn_ca_s*c_caco3aq ;
    double dn_caohsdzn_ca_s    = phi*dc_caohsdzn_ca_s + dphisdzn_ca_s*c_caoh ;
    double dn_ksdzn_ca_s       = dphisdzn_ca_s*c_k;
    
    /* ... with respect to c_k */
    double dc_k                = 0.4e-4 ;/*c_k*((c_k > C_Kn(i)) ? 1 : -1) ;*/
    double c_k2                = c_k + dc_k ;
    double c_oh4               = concentration_oh(c_h2co3,el,zn_ca_s,c_k2,zn_si_s) ;
    double dc_ohsdc_k          = (c_oh4 - c_oh)/dc_k ;
    double dc_hsdc_k           = - c_h*dc_ohsdc_k/c_oh ;
    double dc_hco3sdc_k        = c_hco3*(dc_ohsdc_k/c_oh) ;
    double dc_co3sdc_k         = c_co3*(dc_hco3sdc_k/c_hco3 + dc_ohsdc_k/c_oh) ; 
	                           
    double dc_casdc_k          = - c_ca*dc_co3sdc_k/c_co3 ;
    double dc_cahco3sdc_k      = K_cahco3*(dc_casdc_k*c_hco3 + c_ca*dc_hco3sdc_k) ;

    double dc_h3sio4sdc_k      = - c_h3sio4*dc_hsdc_k/c_h ;
    double dc_cah3sio4sdc_k    = K_cah3sio4*(dc_casdc_k*c_h3sio4 + c_ca*dc_h3sio4sdc_k) ;
    double dc_h2sio4sdc_k      = K_h2sio4*(dc_h3sio4sdc_k*c_oh + dc_ohsdc_k*c_h3sio4) ;
    double dc_cah2sio4sdc_k    = K_cah2sio4*(dc_h2sio4sdc_k*c_ca + dc_casdc_k*c_h2sio4) ;
    double dc_caco3aqsdc_k     = K_caco3aq*(dc_co3sdc_k*c_ca + dc_casdc_k*c_co3) ;				                  
    double dc_caohsdc_k        = K_caoh*(dc_casdc_k*c_oh + dc_ohsdc_k*c_ca) ;					       

    double dn_hco3sdc_k        = phi*dc_hco3sdc_k ;
    double dn_co3sdc_k         = phi*dc_co3sdc_k ;
    double dn_casdc_k          = phi*dc_casdc_k ;
    double dn_cahco3sdc_k      = phi*dc_cahco3sdc_k ;
    double dn_cah3sio4sdc_k    = phi*dc_cah3sio4sdc_k ;
    double dn_h3sio4sdc_k      = phi*dc_h3sio4sdc_k ;
    double dn_h2sio4sdc_k      = phi*dc_h2sio4sdc_k ;
    double dn_cah2sio4sdc_k    = phi*dc_cah2sio4sdc_k ;
    double dn_caco3aqsdc_k     = phi*dc_caco3aqsdc_k ;
    double dn_caohsdc_k        = phi*dc_caohsdc_k ;

    j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_H2CO3+j)   += volume[i]*(dn_h2co3sdc_h2co3 + dn_hco3sdc_h2co3 + dn_co3sdc_h2co3 + dn_cahco3sdc_h2co3 \
                        + dn_caco3aqsdc_h2co3 + dn_ccsdc_h2co3) ;
    K(E_C+j,I_Ca_S+j)   += volume[i]*(dn_h2co3sdzn_ca_s + dn_hco3sdzn_ca_s + dn_co3sdzn_ca_s + dn_cahco3sdzn_ca_s \
                        + dn_caco3aqsdzn_ca_s + dn_ccsdzn_ca_s) ;
    K(E_C+j,I_Si_S+j)   += volume[i]*(dn_h2co3sdzn_si_s + dn_hco3sdzn_si_s + dn_co3sdzn_si_s + dn_cahco3sdzn_si_s \
                        + dn_caco3aqsdzn_si_s) ;
    K(E_C+j,I_K+j)     += volume[i]*(dn_hco3sdc_k + dn_co3sdc_k + dn_cahco3sdc_k + dn_caco3aqsdc_k) ;                          
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_H2CO3+j)  += volume[i]*(dn_casdc_h2co3 + dn_cahco3sdc_h2co3 + dn_cah3sio4sdc_h2co3 + dn_ca_ssdc_h2co3 \
                        + dn_cah2sio4sdc_h2co3 + dn_caco3aqsdc_h2co3 +dn_caohsdc_h2co3) ;
    K(E_Ca+j,I_Ca_S+j)  += volume[i]*(dn_casdzn_ca_s + dn_cahco3sdzn_ca_s + dn_cah3sio4sdzn_ca_s + dn_ca_ssdzn_ca_s \
                        + dn_cah2sio4sdzn_ca_s + dn_caco3aqsdzn_ca_s +dn_caohsdzn_ca_s ) ;
    K(E_Ca+j,I_Si_S+j)  += volume[i]*(dn_casdzn_si_s + dn_cahco3sdzn_si_s + dn_cah3sio4sdzn_si_s + dn_ca_ssdzn_si_s \
                        + dn_cah2sio4sdzn_si_s + dn_caco3aqsdzn_si_s +dn_caohsdzn_si_s ) ;
    K(E_Ca+j,I_K+j)    += volume[i]*(dn_casdc_k + dn_cahco3sdc_k + dn_cah3sio4sdc_k \
                        + dn_cah2sio4sdc_k + dn_caco3aqsdc_k +dn_caohsdc_k ) ;
    /*
      Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
    */
    K(E_Si+j,I_H2CO3+j)  += volume[i]*(dn_h3sio4sdc_h2co3 + dn_h4sio4sdc_h2co3 + dn_cah3sio4sdc_h2co3 + dn_si_ssdc_h2co3\
                        + dn_h2sio4sdc_h2co3 + dn_cah2sio4sdc_h2co3) ;
    K(E_Si+j,I_Ca_S+j)  += volume[i]*(dn_h3sio4sdzn_ca_s + dn_h4sio4sdzn_ca_s + dn_cah3sio4sdzn_ca_s + dn_si_ssdzn_ca_s \
                        + dn_h2sio4sdzn_ca_s + dn_cah2sio4sdzn_ca_s) ;
    K(E_Si+j,I_Si_S+j)  += volume[i]*(dn_h3sio4sdzn_si_s + dn_h4sio4sdzn_si_s + dn_cah3sio4sdzn_si_s + dn_si_ssdzn_si_s \
                        + dn_h2sio4sdzn_si_s + dn_cah2sio4sdzn_si_s) ;
    K(E_Si+j,I_K+j)    += volume[i]*(dn_h3sio4sdc_k + dn_cah3sio4sdc_k + dn_h2sio4sdc_k + dn_cah2sio4sdc_k) ;  
    /*
      Conservation de la charge  : div(w_q) = 0
    */
    /*
      Conservation de K (potassium)  : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
    K(E_K+j,I_H2CO3+j)   += volume[i]*(dn_ksdc_h2co3) ; 
    K(E_K+j,I_Ca_S+j)   += volume[i]*(dn_ksdzn_ca_s) ; 
    K(E_K+j,I_Si_S+j)   += volume[i]*(dn_ksdzn_si_s) ;    
    K(E_K+j,I_K+j)     += volume[i]*(phi*1.) ;

    /* sauvegardes pour les termes de transport */
    Dc_hSDc_h2co3[i]        = dc_hsdc_h2co3 ;
    Dc_ohSDc_h2co3[i]       = dc_ohsdc_h2co3 ;
    Dc_hco3SDc_h2co3[i]     = dc_hco3sdc_h2co3 ;
    Dc_co3SDc_h2co3[i]      = dc_co3sdc_h2co3 ;
    Dc_caSDc_h2co3[i]       = dc_casdc_h2co3 ;
    Dc_cahco3SDc_h2co3[i]   = dc_cahco3sdc_h2co3 ;
    Dc_h3sio4SDc_h2co3[i]   = dc_h3sio4sdc_h2co3 ;
    Dc_h4sio4SDc_h2co3[i]   = dc_h4sio4sdc_h2co3 ;
    Dc_cah3sio4SDc_h2co3[i] = dc_cah3sio4sdc_h2co3 ;
    Dc_h2sio4SDc_h2co3[i]   = dc_h2sio4sdc_h2co3 ;
    Dc_cah2sio4SDc_h2co3[i] = dc_cah2sio4sdc_h2co3 ;
    Dc_caco3aqSDc_h2co3[i]  = dc_caco3aqsdc_h2co3 ;
    Dc_caohSDc_h2co3[i]     = dc_caohsdc_h2co3 ;
    
    Dc_hSDzn_si_s[i]        = dc_hsdzn_si_s ;
    Dc_ohSDzn_si_s[i]       = dc_ohsdzn_si_s ;
    Dc_hco3SDzn_si_s[i]     = dc_hco3sdzn_si_s ;
    Dc_co3SDzn_si_s[i]      = dc_co3sdzn_si_s ;   
    Dc_caSDzn_si_s[i]       = dc_casdzn_si_s ;
    Dc_cahco3SDzn_si_s[i]   = dc_cahco3sdzn_si_s ;
    Dc_h3sio4SDzn_si_s[i]   = dc_h3sio4sdzn_si_s ;
    Dc_h4sio4SDzn_si_s[i]   = dc_h4sio4sdzn_si_s ;
    Dc_cah3sio4SDzn_si_s[i] = dc_cah3sio4sdzn_si_s ;
    Dc_h2sio4SDzn_si_s[i]   = dc_h2sio4sdzn_si_s ;
    Dc_cah2sio4SDzn_si_s[i] = dc_cah2sio4sdzn_si_s ;
    Dc_caco3aqSDzn_si_s[i]  = dc_caco3aqsdzn_si_s ;
    Dc_caohSDzn_si_s[i]     = dc_caohsdzn_si_s ;  

    Dc_hSDzn_ca_s[i]        = dc_hsdzn_ca_s ;
    Dc_ohSDzn_ca_s[i]       = dc_ohsdzn_ca_s ;
    Dc_hco3SDzn_ca_s[i]     = dc_hco3sdzn_ca_s ;
    Dc_co3SDzn_ca_s[i]      = dc_co3sdzn_ca_s ;
    Dc_caSDzn_ca_s[i]       = dc_casdzn_ca_s ;      
    Dc_cahco3SDzn_ca_s[i]   = dc_cahco3sdzn_ca_s ;
    Dc_h3sio4SDzn_ca_s[i]   = dc_h3sio4sdzn_ca_s ;
    Dc_h4sio4SDzn_ca_s[i]   = dc_h4sio4sdzn_ca_s ;
    Dc_cah3sio4SDzn_ca_s[i] = dc_cah3sio4sdzn_ca_s ;
    Dc_h2sio4SDzn_ca_s[i]   = dc_h2sio4sdzn_ca_s ;
    Dc_cah2sio4SDzn_ca_s[i] = dc_cah2sio4sdzn_ca_s ;
    Dc_caco3aqSDzn_ca_s[i]  = dc_caco3aqsdzn_ca_s ;
    Dc_caohSDzn_ca_s[i]     = dc_caohsdzn_ca_s ;
    
    Dc_hSDc_k[i]        = dc_hsdc_k ;
    Dc_ohSDc_k[i]       = dc_ohsdc_k ;
    Dc_hco3SDc_k[i]     = dc_hco3sdc_k ;
    Dc_co3SDc_k[i]      = dc_co3sdc_k ;
    Dc_caSDc_k[i]       = dc_casdc_k ;      
    Dc_cahco3SDc_k[i]   = dc_cahco3sdc_k ;
    Dc_h3sio4SDc_k[i]   = dc_h3sio4sdc_k ;
    Dc_cah3sio4SDc_k[i] = dc_cah3sio4sdc_k ;
    Dc_h2sio4SDc_k[i]   = dc_h2sio4sdc_k ;
    Dc_cah2sio4SDc_k[i] = dc_cah2sio4sdc_k ;
    Dc_caco3aqSDc_k[i]  = dc_caco3aqsdc_k ;
    Dc_caohSDc_k[i]     = dc_caohsdc_k ;    
  }

  /* termes d'ecoulement */
  {
  double tr        = dt*surf/dx ;

  double trf_oh       = tr*KF_OH ;
  double trf_h        = tr*KF_H ;
  double trf_h2co3    = tr*KF_H2CO3 ;
  double trf_hco3     = tr*KF_HCO3 ;
  double trf_co3      = tr*KF_CO3 ;
  double trf_ca       = tr*KF_Ca ;
  double trf_cahco3   = tr*KF_CaHCO3 ;
  double trf_cah3sio4 = tr*KF_CaH3SiO4 ;
  double trf_h3sio4   = tr*KF_H3SiO4 ;
  double trf_h4sio4   = tr*KF_H4SiO4 ;
  double trf_h2sio4   = tr*KF_H2SiO4 ;
  double trf_cah2sio4 = tr*KF_CaH2SiO4 ;
  double trf_caco3aq  = tr*KF_CaCO3aq ;
  double trf_caoh     = tr*KF_CaOH ;
  double trf_k        = tr*KF_K ;                   /* ajouter2 */  
  
  double tre_hco3     = tr*Kpsi_HCO3 ;
  double tre_co3      = tr*Kpsi_CO3 ;
  double tre_ca       = tr*Kpsi_Ca ;
  double tre_cahco3   = tr*Kpsi_CaHCO3 ;
  double tre_cah3sio4 = tr*Kpsi_CaH3SiO4 ;
  double tre_h3sio4   = tr*Kpsi_H3SiO4 ;
  double tre_h2sio4   = tr*Kpsi_H2SiO4 ;
  double tre_caoh     = tr*Kpsi_CaOH ;
  double tre_k        = tr*Kpsi_K ;                  /* ajouter2 */    

  double tre_q        = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h2co3 + trf_hco3*Dc_hco3SDc_h2co3[i] + trf_co3*Dc_co3SDc_h2co3[i] + trf_cahco3*Dc_cahco3SDc_h2co3[i] \
           + trf_caco3aq*Dc_caco3aqSDc_h2co3[i] ;
  }
  K(E_C,I_H2CO3)          += + c[0] ;
  K(E_C,I_H2CO3+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_H2CO3)      += - c[0] ;
  K(E_C+NEQ,I_H2CO3+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDzn_si_s[i] + trf_co3*Dc_co3SDzn_si_s[i] + trf_cahco3*Dc_cahco3SDzn_si_s[i] + trf_caco3aq*Dc_caco3aqSDzn_si_s[i] ;
  }
  K(E_C,I_Si_S)          += + c[0] ;
  K(E_C,I_Si_S+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_Si_S)      += - c[0] ;
  K(E_C+NEQ,I_Si_S+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDzn_ca_s[i] + trf_co3*Dc_co3SDzn_ca_s[i] + trf_cahco3*Dc_cahco3SDzn_ca_s[i] + trf_caco3aq*Dc_caco3aqSDzn_ca_s[i] ;
  }
  K(E_C,I_Ca_S)           += + c[0] ;
  K(E_C,I_Ca_S+NEQ)       += - c[1] ;
  K(E_C+NEQ,I_Ca_S)       += - c[0] ;
  K(E_C+NEQ,I_Ca_S+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDc_k[i] + trf_co3*Dc_co3SDc_k[i] + trf_cahco3*Dc_cahco3SDc_k[i] + trf_caco3aq*Dc_caco3aqSDc_k[i] ;
  }
  K(E_C,I_K)          += + c[0] ;
  K(E_C,I_K+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_K)      += - c[0] ;
  K(E_C+NEQ,I_K+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 + tre_cahco3 ;
  }
  K(E_C,I_psi)          += + c[0] ;
  K(E_C,I_psi+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_psi)      += - c[0] ;
  K(E_C+NEQ,I_psi+NEQ)  += + c[1] ;

  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_h2co3[i] + trf_cahco3*Dc_cahco3SDc_h2co3[i] + trf_cah3sio4*Dc_cah3sio4SDc_h2co3[i] + trf_cah2sio4*Dc_cah2sio4SDc_h2co3[i] \
    + trf_caoh*Dc_caohSDc_h2co3[i] + trf_caco3aq*Dc_caco3aqSDc_h2co3[i] ;
  }
  K(E_Ca,I_H2CO3)         += + c[0] ;
  K(E_Ca,I_H2CO3+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_H2CO3)     += - c[0] ;
  K(E_Ca+NEQ,I_H2CO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_si_s[i] + trf_cahco3*Dc_cahco3SDzn_si_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + trf_caoh*Dc_caohSDzn_si_s[i] + trf_caco3aq*Dc_caco3aqSDzn_si_s[i] ;
  }
  K(E_Ca,I_Si_S)         += + c[0] ;
  K(E_Ca,I_Si_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Si_S)     += - c[0] ;
  K(E_Ca+NEQ,I_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_ca_s[i] + trf_cahco3*Dc_cahco3SDzn_ca_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + trf_caoh*Dc_caohSDzn_ca_s[i] + trf_caco3aq*Dc_caco3aqSDzn_ca_s[i] ;
  }
  K(E_Ca,I_Ca_S)         += + c[0] ;
  K(E_Ca,I_Ca_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Ca_S)     += - c[0] ;
  K(E_Ca+NEQ,I_Ca_S+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_k[i] + trf_cahco3*Dc_cahco3SDc_k[i] + trf_cah3sio4*Dc_cah3sio4SDc_k[i] + trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + trf_caoh*Dc_caohSDc_k[i] + trf_caco3aq*Dc_caco3aqSDc_k[i] ;
  }
  K(E_Ca,I_K)         += + c[0] ;
  K(E_Ca,I_K+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_K)     += - c[0] ;
  K(E_Ca+NEQ,I_K+NEQ) += + c[1] ; 

  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_cahco3 + tre_cah3sio4 + tre_caoh;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_h2co3[i] + trf_h4sio4*Dc_h4sio4SDc_h2co3[i] + trf_cah3sio4*Dc_cah3sio4SDc_h2co3[i] + trf_cah2sio4*Dc_cah2sio4SDc_h2co3[i] \
    + trf_h2sio4*Dc_h2sio4SDc_h2co3[i] ;
  }
  K(E_Si,I_H2CO3)         += + c[0] ;
  K(E_Si,I_H2CO3+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_H2CO3)     += - c[0] ;
  K(E_Si+NEQ,I_H2CO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDzn_si_s[i] + trf_h4sio4*Dc_h4sio4SDzn_si_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i]  + trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + trf_h2sio4*Dc_h2sio4SDzn_si_s[i];
  }
  K(E_Si,I_Si_S)         += + c[0] ;
  K(E_Si,I_Si_S+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_Si_S)     += - c[0] ;
  K(E_Si+NEQ,I_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] + trf_h4sio4*Dc_h4sio4SDzn_ca_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i]  + trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + trf_h2sio4*Dc_h2sio4SDzn_ca_s[i];
  }
  K(E_Si,I_Ca_S)         += + c[0] ;
  K(E_Si,I_Ca_S+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_Ca_S)     += - c[0] ;
  K(E_Si+NEQ,I_Ca_S+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_k[i] + trf_cah3sio4*Dc_cah3sio4SDc_k[i]  + trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + trf_h2sio4*Dc_h2sio4SDc_k[i];
  }
  K(E_Si,I_K)         += + c[0] ;
  K(E_Si,I_K+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_K)     += - c[0] ;
  K(E_Si+NEQ,I_K+NEQ) += + c[1] ; 

  for(i=0;i<2;i++){
    c[i] = tre_h3sio4 + tre_cah3sio4 + tre_h2sio4 ;
  }
  K(E_Si,I_psi)          += + c[0] ;
  K(E_Si,I_psi+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_psi)      += - c[0] ;
  K(E_Si+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDc_h2co3[i] + z_oh*trf_oh*Dc_ohSDc_h2co3[i] + z_hco3*trf_hco3*Dc_hco3SDc_h2co3[i] + z_co3*trf_co3*Dc_co3SDc_h2co3[i] \
           + z_ca*trf_ca*Dc_caSDc_h2co3[i] + z_cahco3*trf_cahco3*Dc_cahco3SDc_h2co3[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_h2co3[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_h2co3[i] + z_caoh*trf_caoh*Dc_caohSDc_h2co3[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_h2co3[i];
  }
  K(E_q,I_H2CO3)           += + c[0] ;
  K(E_q,I_H2CO3+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_H2CO3)       += - c[0] ;
  K(E_q+NEQ,I_H2CO3+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_si_s[i] + z_oh*trf_oh*Dc_ohSDzn_si_s[i] + z_hco3*trf_hco3*Dc_hco3SDzn_si_s[i] + z_co3*trf_co3*Dc_co3SDzn_si_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_si_s[i] + z_cahco3*trf_cahco3*Dc_cahco3SDzn_si_s[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDzn_si_s[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + z_caoh*trf_caoh*Dc_caohSDzn_si_s[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDzn_si_s[i];
  }
  K(E_q,I_Si_S)           += + c[0] ;
  K(E_q,I_Si_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Si_S)       += - c[0] ;
  K(E_q+NEQ,I_Si_S+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_ca_s[i] + z_oh*trf_oh*Dc_ohSDzn_ca_s[i] + z_hco3*trf_hco3*Dc_hco3SDzn_ca_s[i] + z_co3*trf_co3*Dc_co3SDzn_ca_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_ca_s[i] + z_cahco3*trf_cahco3*Dc_cahco3SDzn_ca_s[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + z_caoh*trf_caoh*Dc_caohSDzn_ca_s[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDzn_ca_s[i];
  }
  K(E_q,I_Ca_S)           += + c[0] ;
  K(E_q,I_Ca_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Ca_S)       += - c[0] ;
  K(E_q+NEQ,I_Ca_S+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = z_k*trf_k + z_h*trf_h*Dc_hSDc_k[i] + z_oh*trf_oh*Dc_ohSDc_k[i] + z_hco3*trf_hco3*Dc_hco3SDc_k[i] + z_co3*trf_co3*Dc_co3SDc_k[i] \
           + z_ca*trf_ca*Dc_caSDc_k[i] + z_cahco3*trf_cahco3*Dc_cahco3SDc_k[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_k[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_k[i] + z_caoh*trf_caoh*Dc_caohSDc_k[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_k[i];
  }
  K(E_q,I_K)           += + c[0] ;
  K(E_q,I_K+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_K)       += - c[0] ;
  K(E_q+NEQ,I_K+NEQ)   += + c[1] ; 

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_psi)          += + c[0] ;
  K(E_q,I_psi+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_psi)      += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)  += + c[1] ;
  
    /*
  Conservation de K (potassium)  : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  for(i=0;i<2;i++){
    c[i] = tre_k ;
  }
  K(E_K,I_psi)          += + c[0] ;
  K(E_K,I_psi+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_psi)      += - c[0] ;
  K(E_K+NEQ,I_psi+NEQ)  += + c[1] ;
    
   
  for(i=0;i<2;i++){
    c[i] = trf_k ;
  }
  K(E_K,I_K)          += + c[0] ;
  K(E_K,I_K+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_K)      += - c[0] ;
  K(E_K+NEQ,I_K+NEQ)  += + c[1] ; 
  
 }

#if (U_H2CO3 == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_H2CO3)     *= Ln10*C_H2CO3(0) ;
    K(i,I_H2CO3+NEQ) *= Ln10*C_H2CO3(1) ;
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
  R(0,E_C) -= volume[0]*(N_C(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C(1) - N_Cn(1)) - dt*surf*W_C ;
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
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[Element_MaxNbOfNodes] ;
  int    i,nso ;
  double zero = 0. ;
#define UNKNOWN_s(i)   UNKNOWN(j,i)
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }

  /* if(el.dim < dim) return(0) ; */
  
  /*
    Donnees
  */
  phi0      = Element_GetProperty(el)[pm("porosite")] ;
  c_h2co3_eq  = Element_GetProperty(el)[pm("C_H2CO3_eq")] ;

  /* initialisation */
  nso = 28 ;
  for(i = 0 ; i < nso ; i++) {
    int j ;
    for(j = 0 ; j < 9 ; j++) Result_GetValue(r+i)[j] = zero ;
  }


  /* output quantities */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double xm = (x0 + x1)*0.5 ;
    int j = (Element_IsSubmanifold(el)) ? 0 : ((s[0] < xm) ? 0 : 1) ;
    /* molarities */
#if (U_H2CO3 == LOG_RHO)
    double c_h2co3      = exp(Ln10*UNKNOWN_s(I_H2CO3)) ;
#else
    double c_h2co3      =  UNKNOWN_s(I_H2CO3) ;
#endif
    double zc_h2co3     = c_h2co3/c_h2co3_eq ;
    double zn_ca_s    =  UNKNOWN_s(I_Ca_S) ;
    double c_k        =  UNKNOWN_s(I_K) ;
    double zn_si_s    =  UNKNOWN_s(I_Si_S) ;
    
    double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;
    
    double c_h4sio4   = P_Sam ;
    double c_oh       = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;    
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* charge density */
    double c_q = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 \
    + z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k ;
    /* solid contents */
    double n_ch       = N_CH(j) ;
    double n_cc       = N_CC(j) ;
    double n_si_s     = N_Si_S(j)  ;
    
    /* porosity */
    double n_ch0      = N_CH0(j) ;
    double n_cc0      = N_CC0(j) ;
    double n_si_s0    = N_Si_S0(j) ;
    double zp_ch0     = ZP_CH0(j) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double v_csh      = V_CSH(zp_ch) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) \
                      + v_csh0*n_si_s0 - v_csh*n_si_s ;

    double psi        = UNKNOWN_s(I_psi) ;
    double ph         = 14 + log(c_oh)/log(10.) ;
    double pk_ch      = log10(zp_ch) ;
    double x_csh      = X_CSH(zp_ch) ;

    i = 0 ;
    Result_Store(r + i++,&c_h2co3,"c_h2co3",1) ;
    Result_Store(r + i++,&ph,"ph",1) ;
    Result_Store(r + i++,&c_h2co3,"c_h2co3",1) ;
    Result_Store(r + i++,&c_hco3,"c_hco3",1) ;
    Result_Store(r + i++,&c_co3,"c_co3",1) ;
    Result_Store(r + i++,&c_ca,"c_ca",1) ;
    Result_Store(r + i++,&c_cahco3,"c_cahco3",1) ;
    Result_Store(r + i++,&c_cah3sio4,"c_cah3sio4",1) ;
    Result_Store(r + i++,&c_h3sio4,"c_h3sio4",1) ;
    Result_Store(r + i++,&c_h4sio4,"c_h4sio4",1) ;
    Result_Store(r + i++,&n_ch,"n_ch",1) ;
    Result_Store(r + i++,&n_cc,"n_cc",1) ;
    Result_Store(r + i++,&n_si_s,"n_si_s",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&c_oh,"c_oh",1) ;
    Result_Store(r + i++,&psi,"potentiel_electrique",1) ;
    Result_Store(r + i++,&c_q,"charge",1) ;
    Result_Store(r + i++,&zn_ca_s,"zn_ca_s",1) ;
    Result_Store(r + i++,&zc_h2co3,"zc_h2co3",1) ;
    Result_Store(r + i++,&c_caco3aq,"c_caco3aq",1) ;
    Result_Store(r + i++,&c_caoh,"c_caoh",1) ;
    Result_Store(r + i++,&c_h2sio4,"c_h2sio4",1) ;
    Result_Store(r + i++,&c_cah2sio4,"c_cah2sio4",1) ;
    Result_Store(r + i++,&c_k,"c_k",1) ;
    Result_Store(r + i++,&pk_ch,"pk_ch",1) ;
    Result_Store(r + i++,&zn_si_s,"zn_si_s",1) ;
    Result_Store(r + i++,&v_csh,"V_CSH",1) ;
    Result_Store(r + i++,&x_csh,"C/S",1) ;
  }
  
  if(i != nso) arret("so53") ;

  return(nso) ;
}


void transfert(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int    i ;
  /*
    Donnees
  */
  phi0      = Element_GetProperty(el)[pm("porosite")] ;
  c_h2co3_eq  = Element_GetProperty(el)[pm("C_H2CO3_eq")] ;

  /* initialisation */
  for(i=0;i<NVE;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<2;i++) {
    /* molarites */
    double c_h2co3      = C_H2CO3(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2co3     = c_h2co3/c_h2co3_eq ;
    double c_k        = C_K(i) ;
    double zn_si_s    = ZN_Si_S(i); 

    double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ; 
    double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;

    double c_h4sio4   = P_Sam ;
    double c_oh       = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* solid contents */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_si_s     = N_Si_S(i) ;

    /* porosity */
    double n_cc0      = N_CC0(i) ;
    double n_ch0      = N_CH0(i) ;
    double n_si_s0    = N_Si_S0(i) ;
    double zp_ch0     = ZP_CH0(i) ;

    double v_csh0     = V_CSH(zp_ch0) ;
    double v_csh      = V_CSH(zp_ch) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) \
                      + v_csh0*n_si_s0 - v_csh*n_si_s ;

    /* tortuosite liquide */
    double iff    = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_H2CO3      += d_h2co3*iff ;
    KF_HCO3       += d_hco3*iff ;
    KF_CO3        += d_co3*iff ;

    KF_Ca         += d_ca*iff ;
    KF_CaHCO3     += d_cahco3*iff ;
    KF_CaH3SiO4   += d_cah3sio4*iff ;

    KF_H3SiO4     += d_h3sio4*iff ;
    KF_H4SiO4     += d_h4sio4*iff ;
    KF_H2SiO4     += d_h2sio4*iff ;
    KF_CaH2SiO4   += d_cah2sio4*iff ;
    KF_CaCO3aq    += d_caco3aq*iff ;
    KF_CaOH       += d_caoh*iff ;
    
    KF_K          += d_k*iff;


    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HCO3     += FsRT*KF_HCO3*z_hco3*c_hco3 ;
    Kpsi_CO3      += FsRT*KF_CO3*z_co3*c_co3 ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaHCO3   += FsRT*KF_CaHCO3*z_cahco3*c_cahco3 ;
    Kpsi_CaH3SiO4 += FsRT*KF_CaH3SiO4*z_cah3sio4*c_cah3sio4 ;
    Kpsi_H3SiO4   += FsRT*KF_H3SiO4*z_h3sio4*c_h3sio4 ;
    Kpsi_H2SiO4   += FsRT*KF_H2SiO4*z_h2sio4*c_h2sio4 ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    
    Kpsi_K        += FsRT*KF_K*z_k*c_k ;

    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca + z_cahco3*Kpsi_CaHCO3 + z_h3sio4*Kpsi_H3SiO4 \
                   + z_cah3sio4*Kpsi_CaH3SiO4 + z_caoh*Kpsi_CaOH + z_h2sio4*Kpsi_H2SiO4 + z_k*Kpsi_K ;
  }
  
  /* moyenne */
  for(i=0;i<NVE;i++) va[i] *= 0.5 ;
}


void flux(Element_t *el,double **u)
/* Les flux (f) */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double r_h[2],r_oh[2] ;
  double r_h2co3[2],r_hco3[2],r_co3[2],r_caco3aq[2],r_caoh[2] ;
  double r_ca[2],r_cahco3[2],r_h3sio4[2],r_h4sio4[2],r_cah3sio4[2],r_cah2sio4[2],r_h2sio4[2] ;
  double r_k[2] ;

  int    i ;

  for(i=0;i<2;i++) {
    double c_h2co3      = C_H2CO3(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2co3     = c_h2co3/c_h2co3_eq ;
    double c_k        = C_K(i);
    double zn_si_s    = ZN_Si_S(i); 

    double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ;    
    double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;
 
    double c_h4sio4   = P_Sam ;
    double c_oh       = concentration_oh(c_h2co3,el,zn_ca_s,c_k,zn_si_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ; 
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    r_oh[i]       = c_oh ;
    r_h[i]        = c_h ;
    r_h2co3[i]    = c_h2co3 ;
    r_hco3[i]     = c_hco3 ;
    r_co3[i]      = c_co3 ;
    r_ca[i]       = c_ca ;
    r_cahco3[i]   = c_cahco3 ;
    r_h3sio4[i]   = c_h3sio4 ;
    r_h4sio4[i]   = c_h4sio4 ;
    r_cah3sio4[i] = c_cah3sio4 ;
    r_h2sio4[i]   = c_h2sio4 ;
    r_cah2sio4[i] = c_cah2sio4 ;
    r_caoh[i]     = c_caoh ;
    r_caco3aq[i]  = c_caco3aq ;
    r_k[i]        = c_k ;
  }

  /* Gradients */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    double grd_h        = (r_h[1]        - r_h[0]       )/dx ;
    double grd_oh       = (r_oh[1]       - r_oh[0]      )/dx ;
    double grd_h2co3    = (r_h2co3[1]    - r_h2co3[0]   )/dx ;
    double grd_hco3     = (r_hco3[1]     - r_hco3[0]    )/dx ;
    double grd_co3      = (r_co3[1]      - r_co3[0]     )/dx ;
    double grd_cahco3   = (r_cahco3[1]   - r_cahco3[0]  )/dx ;
    double grd_ca       = (r_ca[1]       - r_ca[0]      )/dx ;
    double grd_cah3sio4 = (r_cah3sio4[1] - r_cah3sio4[0])/dx ;
    double grd_h3sio4   = (r_h3sio4[1]   - r_h3sio4[0]  )/dx ;
    double grd_h4sio4   = (r_h4sio4[1]   - r_h4sio4[0]  )/dx ;
    double grd_h2sio4   = (r_h2sio4[1]   - r_h2sio4[0]  )/dx ;
    double grd_cah2sio4 = (r_cah2sio4[1] - r_cah2sio4[0])/dx ;
    double grd_caco3aq  = (r_caco3aq[1]  - r_caco3aq[0] )/dx ;
    double grd_caoh     = (r_caoh[1]     - r_caoh[0]    )/dx ;
    double grd_k        = (r_k[1]        - r_k[0]       )/dx ;
    
    double grd_psi      = (PSI(1)        - PSI(0)       )/dx ;
    
    /* Flux */
    double w_h2co3    = - KF_H2CO3*grd_h2co3   ;
    double w_hco3     = - KF_HCO3*grd_hco3          - Kpsi_HCO3*grd_psi ;
    double w_co3      = - KF_CO3*grd_co3            - Kpsi_CO3*grd_psi      ;
    double w_cahco3   = - KF_CaHCO3*grd_cahco3      - Kpsi_CaHCO3*grd_psi ;
    double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
    double w_cah3sio4 = - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
    double w_h3sio4   = - KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi ;
    double w_h2sio4   = - KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi ;
    double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi ;
    double w_h4sio4   = - KF_H4SiO4*grd_h4sio4 ;
    double w_cah2sio4 = - KF_CaH2SiO4*grd_cah2sio4 ;
    double w_caco3aq  = - KF_CaCO3aq*grd_caco3aq ;    
    double w_k        = - KF_K*grd_k                - Kpsi_K*grd_psi ; 
    
    double w_q        = - z_h*KF_H*grd_h		      \
                        - z_oh*KF_OH*grd_oh		      \
                        - z_hco3*KF_HCO3*grd_hco3             \
                        - z_co3*KF_CO3*grd_co3		      \
                        - z_ca*KF_Ca*grd_ca		      \
                        - z_cahco3*KF_CaHCO3*grd_cahco3	      \
                        - z_h3sio4*KF_H3SiO4*grd_h3sio4	      \
                        - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                        - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                        - z_caoh*KF_CaOH*grd_caoh \
                        - z_k*KF_K*grd_k \
                        - Kpsi_q*grd_psi ;

    W_C     = w_h2co3 + w_hco3 + w_co3 + w_cahco3 + w_caco3aq ;
    W_Ca    = w_ca + w_cahco3 + w_cah3sio4 + w_caco3aq + w_caoh + w_cah2sio4 ;
    W_Si    = w_h3sio4 + w_h4sio4 + w_cah3sio4 + w_cah2sio4 + w_h2sio4 ;
    W_q     = w_q ;
    W_K     = w_k ;
  }
}


double concentration_oh(double c_h2co3,elem_t *el,double zn_ca_s,double c_k,double zn_si_s)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_h2co3, zn_si_s et zn_ca_s sont fixes */
  double c_h2co3_eq   = C_H2CO3_eq ;
  double zc_h2co3     = c_h2co3/c_h2co3_eq ;
  /* les produits sont donc aussi fixes */
  double P_CC       = IAP_CC(zc_h2co3,zn_ca_s) ;
  double P_CH       = IAP_CH(zc_h2co3,zn_ca_s) ;
  double zp_ch      = P_CH/K_CH ;
  double P_Sam      = IAP_SH(zp_ch,zn_si_s) ;
  
  double c_h4sio4   = P_Sam ;

  /*
  rappel des expressions c_i = A_i*(c_h)**n   : n
     c_h        = K_h2o/c_oh                     : +1
     c_hco3     = K_hco3*c_h2co3/c_h             : -1
     c_co3      = K_co3*c_hco3/c_h               : -2
     c_ca       = P_CC/c_co3                     : +2
     c_cahco3   = K_cahco3*c_ca*c_hco3           : +1
     c_h4sio4   = P_Sam                          :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)        : -1
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4       : +1
     c_caco3aq  = K_caco3aq*c_ca*c_co3 ;         :  0      
     c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;       : -2       
     c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;     :  0      
     c_caoh     = K_caoh*c_ca*c_oh ;             : +1       
  */
  double A_hco3     = K_hco3*c_h2co3 ;
  double A_co3      = K_co3*A_hco3 ;
  double A_ca       = P_CC/A_co3 ;
  double A_cahco3   = K_cahco3*A_ca*A_hco3 ;
  double A_h3sio4   = c_h4sio4/K_h4sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;
  double A_h2sio4   = K_h2sio4*A_h3sio4*K_h2o ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahco3*A_cahco3 + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_k*c_k ;
  double d = z_oh*K_h2o + z_hco3*A_hco3 + z_h3sio4*A_h3sio4 ;
  double e = z_co3*A_co3 + z_h2sio4*A_h2sio4 ;

  /* a > 0 ; b > 0 ; c > 0 ; d < 0 ; e < 0 */
  double c_h = poly4(a,b,c,d,e) ;
 
  return(K_h2o/c_h) ;
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
    if(i++ > 20) {
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      arret("poly4 : non convergence") ;
    }
  } while(err > tol) ;
  return(x) ;
}
