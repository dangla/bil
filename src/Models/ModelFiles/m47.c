/*
  Trois reactions avec cinetique :
   1. acido-basique : HCO3- + H20 <-> H2CO3 + OH-    (loi de Danckwerts)
   2. dissolution de portlandite : Ca(OH)2 <-> Ca2+ + 2OH-    (loi en log)
   3. Carbo CSH  : CSH + 2H2CO3 <-> 3CaCO3 + 2SiO2.3H2O + 3H2O (loi macro)
  Cristaux de portlandite spherique, 
  Diffusion gaz MT,
  Transferts avec diffusion ionique
  Electroneutralite
  Ajout des alcalins
  Ajout d'espèces CaOH, CaHCO3, CaCO3aq, CaOH2aq
  Prise en commpte de la quantitée totale en alcalins
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  47
#define TITLE "Carbonatation du beton (2009)"
#define AUTHORS "Thiery"

#include "OldMethods.h"

/* Macros */
#define NEQ     (8)
#define NVE     (49)
#define NVI     (27)

#define E_C     (0)
#define E_q     (1)
#define E_mass  (2)
#define E_Ca    (3)
#define E_k     (4)
#define E_el    (5)
#define E_Na    (6)
#define E_K     (7)

#define I_CO2   (0)
#define I_OH    (5)
#define I_P_l   (2)
#define I_CaCO3 (3)
#define I_HCO3  (4)
#define I_psi   (1)
#define I_Na    (6)
#define I_K     (7)


#define X_CO2(n)   (u[(n)][I_CO2])
#define X_OH(n)    (u[(n)][I_OH])
#define N_CaCO3(n) (u[(n)][I_CaCO3])
#define P_l(n)     (u[(n)][I_P_l])
#define X_HCO3(n)  (u[(n)][I_HCO3])
#define PSI(n)     (u[(n)][I_psi])
#define X_Na(n)    (u[(n)][I_Na])
#define X_K(n)     (u[(n)][I_K])

#define N_C(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define M(n)       (f[(4+n)])
#define N_Ca(n)    (f[(6+n)])
#define N_k(n)     (f[(8+n)])
#define N_Na(n)    (f[(10+n)])
#define N_K(n)     (f[(12+n)])
#define XI(n)      (f[(14+n)])
#define W_C        (f[16])
#define W_q        (f[17])
#define W_m        (f[18])
#define W_Ca       (f[19])
#define W_k        (f[20])
#define W_Na       (f[21])
#define W_K        (f[22])
#define N_CaOH2(n) (f[(23+n)])
#define N_CSH(n)   (f[(25+n)])

#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define M_n(n)     (f_n[(4+n)])
#define N_Can(n)   (f_n[(6+n)])
#define N_kn(n)    (f_n[(8+n)])
#define N_Nan(n)   (f_n[(10+n)])
#define N_Kn(n)    (f_n[(12+n)])

#define N_CaOH2n(n) (f_n[(23+n)])
#define N_CSHn(n)   (f_n[(25+n)])


#define KD_Ca      (va[(0)])
#define KD_OH      (va[(1)])
#define KD_H       (va[(2)])
#define KD_H2CO3   (va[(3)])
#define KD_HCO3    (va[(4)])
#define KD_CO3     (va[(5)])
#define KD_Na      (va[(6)])
#define KD_NaOH    (va[(7)])
#define KD_NaHCO3  (va[(8)])
#define KD_NaCO3   (va[(9)])
#define KD_m       (va[(10)])
#define KD_K       (va[(11)])
#define KD_KOH     (va[(12)])
#define KD_CaOH    (va[(13)])
#define KD_CaHCO3    (va[(14)])
#define KD_CaCO3aq (va[(15)])
#define KD_CaOH2aq (va[(16)])

#define KF_CO2     (va[(17)])
#define KF_Ca      (va[(18)])
#define KF_OH      (va[(19)])
#define KF_H       (va[(20)])
#define KF_H2CO3   (va[(21)])
#define KF_HCO3    (va[(22)])
#define KF_CO3     (va[(23)])
#define KF_Na      (va[(24)])
#define KF_NaOH    (va[(25)])
#define KF_NaHCO3  (va[(26)])
#define KF_NaCO3   (va[(27)])
#define KF_K       (va[(28)])
#define KF_KOH     (va[(29)])
#define KF_CaOH    (va[(30)])
#define KF_CaHCO3    (va[(31)])
#define KF_CaCO3aq (va[(32)])
#define KF_CaOH2aq (va[(33)])


#define Kpsi_Ca    (va[(34)])
#define Kpsi_OH    (va[(35)])
#define Kpsi_H     (va[(36)])
#define Kpsi_HCO3  (va[(37)])
#define Kpsi_CO3   (va[(38)])
#define Kpsi_Na    (va[(39)])
#define Kpsi_NaCO3 (va[(40)])
#define Kpsi_q     (va[(41)])
#define Kpsi_K     (va[(42)])
#define Kpsi_CaOH  (va[(43)])
#define Kpsi_CaHCO3 (va[(44)])

#define DN_CSHSDT(n)   (va[(45+n)])
#define DN_CaOH2SDT(n) (va[(47+n)])


/* les valences */
#define z_ca    (2.)
#define z_h     (1.)
#define z_oh    (-1.)
#define z_hco3  (-1.)
#define z_co3   (-2.)
#define z_na    (1.)
#define z_naco3 (-1.)
#define z_caoh   (1.)
#define z_cahco3 (1.)
#define z_k      (1.)
#define z_caoh   (1.)
#define z_cahco3 (1.)

/* volumes molaires partiels des ions (dm3/mole) */
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (23.50e-3)    /* (23.50e-3)  d'apres TQN */
#define v_h2o      (18.e-3)
#define v_h2co3    (50.e-3)
#define v_hco3     (50.e-3)
#define v_co3      (-2.3e-3)     /* (-2.3e-3) d'apres 3w */
#define v_ca       (-18.7e-3)    /* (-18.7e-3 d'apres */
#define v_na       (22.4e-3)
#define v_naoh     (22.4e-3)
#define v_nahco3   (22.4e-3) 
#define v_naco3    (22.4e-3) 

#define v_k        (43.93e-3)
#define v_koh      (27.44e-3)

#define v_caoh     (26.20e-3)  /*a modifier*/
#define v_cahco3   (26.20e-3)  /*a modifier*/
#define v_caco3aq  (26.20e-3)  /*a modifier*/
#define v_caoh2aq  (26.20e-3)  /*a modifier*/

/* volumes molaires solides (dm3/mole) */
#define v_caoh2    (33.e-3)
#define v_caco3    (37.e-3)
#define v_csh      (72.e-3)
#define v_caco3csh  v_caco3

/* Masses molaires (unite arbitraire = M_H) */
#define M_Ca       (40.1)
#define M_H2CO3    (62.)
#define M_HCO3     (61.)
#define M_CO3      (60.)
#define M_OH       (17.)
#define M_H        (1.)
#define M_H2O      (18.)
#define M_Na	   (23.)
#define M_NaOH	   (40.)
#define M_NaHCO3   (84.)
#define M_NaCO3	   (83.)

#define M_CO2      (44.)

#define M_CaOH2    (74.)
#define M_CaCO3    (100.)
#define M_CaCO3CSH (100.)
#define M_CSH      (342.) /* 3(CaO).2(SiO2).3(H2O) */
#define M_gel      (174.) /* 2(SiO2).3(H2O) */

#define M_K        (39.)  
#define M_KOH      (56.)

#define M_CaOH     (57.)
#define M_CaHCO3   (101.)
#define M_CaCO3aq  (100.)
#define M_CaOH2aq  (74.)


/* coefficients de diffusion moleculaire (dm2/s) */
#define d_oh       (5.273e-7)    /* (5.273e-7) d'apres TQN */
#define d_h        (9.310e-7)
#define d_ca       (7.92e-8)
#define d_h2co3    (7.2e-8)
#define d_hco3     (11.8e-8)
#define d_co3      (9.55e-8)
#define d_na       (1.33e-7)
#define d_naoh     (1.33e-7)
#define d_nahco3   (1.33e-7)
#define d_naco3    (1.33e-7)
#define d_k        (1.957e-7)
#define d_koh      (1.957e-7)

#define d_caoh     (7.92e-8) /*a modifier */
#define d_cahco3   (7.92e-8) /*a modifier */
#define d_caco3aq  (7.92e-8) /*a modifier */
#define d_caoh2aq  (7.92e-8) /*a modifier */

#define d_co2      (1.6e-3)

/* constantes d'equilibre (ref = 1 mole/L) */
#define k_e        (1.e-14)                /* autoprotolyse de l'eau */
#define k_h        (1.)                    /* cste de Henry */
#define k_co3      (4.570881896148751e3)   /* Equilibre de HCO3  <-> CO3 */
#define k_ca       (3.890451449942805e-9)  /* Equilibre de CaCO3 */
#define k_1        (2.187761623949552e-8)  /* Equilibre de H2CO3 <-> HCO3 */
#define k_2        (6.456542290346550e-6)  /* Equilibre de Ca(OH)2 */
#define k_naoh     (6.60693448e-15)        /* Equilibre de Na <-> NaOH */
#define k_nahco3   (0.5623413252)          /* Equilibre de Na <-> NaHCO3  */
#define k_naco3    (8.72971368e-10)        /* Equilibre de Na <-> NaCO3 */
#define k_caoh     (1.65958691e-13)        /* Equilibre de Ca <-> CaOH */
#define k_cahco3   (12.76438809)           /* Equilibre de Ca <-> CaHCO3 */
#define k_caco3    (7.852356346e-8)        /* Equilibre de Ca <-> CaCO3aqueux */
#define k_koh      (3.4673685e-15)         /* Equilibre de K + H20 <-> KOH + H */
#define k_caoh2    (1.)                    /* Equilibre de Ca(OH)2 en solution */

/* viscosite (Pa.s) */
#define mu_l       (1.e-3)

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (C/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (J/mole) */

/*puissance de la loi cinétique de dissolution de la portlandite*/
#define Xp    (0.5)

/* Fonctions */
static int    pm(const char *s) ;
static double concentration_oh(double,double,double) ;
static double concentration_oh_equilibre_approx(double,double) ;
static double concentration_oh_equilibre(double,double) ;
static double concentration_na_equilibre(double ,double ) ;
static double concentration_k_equilibre(double,double);
static double dn1_caoh2sdt(double,double) ;
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
/* Parametres */
static double phii,k_int,a_1,a_2,c_2,n_caoh20,n_csh0,T_csh ;
/*
static double mu_l ;
static double k_e,k_h,k_ca,k_co3,k_naoh,k_nahco3,k_naco3,k_1,k_2 ;
static double d_co2,d_ca,d_oh,d_h,d_h2co3,d_hco3,d_co3,d_na,d_naoh,d_nahco3,d_naco3 ;
static double v_ca,v_h2o,v_hco3,v_h,v_oh,v_h2co3,v_co3,v_na,v_na,v_naoh,v_nahco3,v_naco3 ;
static double v_csh,v_caoh2,v_caco3 ;
*/
static double p_g = 0. ;

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0)     return (0) ;
  else if(strcmp(s,"k_int") == 0)   return (1) ;
  else if(strcmp(s,"N_CaOH2") == 0) return (2) ;
  else if(strcmp(s,"T_csh") == 0)   return (3) ;
  else if(strcmp(s,"N_CSH") == 0)   return (4) ;
  else if(strcmp(s,"A_1") == 0)     return (5) ;
  else if(strcmp(s,"A_2") == 0)     return (6) ;
  else if(strcmp(s,"C_2") == 0)     return (7) ;
  else if(strcmp(s,"R_CaOH2") == 0) return (8) ;
  else if(strcmp(s,"courbes") == 0) return (9) ;
  else return(-1) ;
}

int dm47(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 10 ;
  
  mat->neq      = NEQ ;

  strcpy(mat->eqn[E_C],   "carbone") ;
  strcpy(mat->eqn[E_q],   "charge") ;
  strcpy(mat->eqn[E_mass],"masse") ;
  strcpy(mat->eqn[E_Ca],  "calcium") ;
  strcpy(mat->eqn[E_k],   "E_k") ;
  strcpy(mat->eqn[E_el],  "E_el") ;
  strcpy(mat->eqn[E_Na],  "sodium") ;
  strcpy(mat->eqn[E_K],  "potassium") ;

  strcpy(mat->inc[I_CO2],  "c_co2") ;
  strcpy(mat->inc[I_OH],   "c_oh") ;
  strcpy(mat->inc[I_P_l],  "p_l") ;
  strcpy(mat->inc[I_CaCO3],"c_caco3") ;
  strcpy(mat->inc[I_HCO3], "c_hco3") ;
  strcpy(mat->inc[I_psi],  "psi") ;
  strcpy(mat->inc[I_Na],   "c_na") ;
  strcpy(mat->inc[I_K],   "c_k") ;

  {
    /* Initialisation automatique a partir des donnees 
       tirees de la these de Mickael Thiery */
    double h     = 5.6e-6 ;  /* (moles/dm2/s) these MT p 223 */
    double D     = 6.5e-15 ; /* (moles/dm/s) these MT p 253 */
    double t_csh = 3000 ;    /* (s) these MT p 230 */
    double R_0   = 30.e-5 ; /* rayon du cristal de CaOH2 (dm) */
    double a_1   = 6000 ;     /* (dm/mole/s) these MT p 253 */
    double c_2   = h*R_0/D ;     /* (1/dm) these MT p 228 */

    mat->pr[pm("A_1")] = a_1 ;
    mat->pr[pm("C_2")] = c_2 ;
    mat->pr[pm("T_csh")] = k_h/t_csh ;

    dmat(mat,ficd,pm,n_donnees) ;
    
    R_0 = mat->pr[pm("R_CaOH2")] ;
    if(R_0 > 0.) { /* si R_0 est donne */
      double n_caoh20 = mat->pr[pm("N_CaOH2")] ; /* contenu en CaOH2 */
      double a_2 = 3*h/R_0*n_caoh20*v_caoh2 ; /* (dm/mole/s) these MT p 227 */
      double c_2 = h*R_0/D ;     /* (1/dm) these MT p 228 */
      mat->pr[pm("A_2")] = a_2 ;
      mat->pr[pm("C_2")] = c_2 ;
    }
  }
  
  return(mat->n) ;
}

int qm47(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 7 equations:\n\
\t- la conservation de la masse de C      (c_co2)\n\
\t- la conservation de la charge          (psi)\n\
\t- la conservation de la masse totale    (p_l)\n\
\t- la conservation de la masse de Ca     (c_caco3)\n\
\t- 1 equation de cinetique               (c_hco3)\n\
\t- Electroneutralite                     (c_oh)\n\
\t- la conservation de la masse de Na     (c_na)\n\ 
\t- la conservation de la masse de K      (c_k)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n\
\t pression : Pa !\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"k_int = 5e-19     # Permeabilite intrinseque (dm2)\n") ;
  fprintf(ficd,"N_CaOH2 = 6.1     # Contenu initial en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"N_CSH = 2.4       # contenu initial en CSH (moles/L)\n") ;
  fprintf(ficd,"T_csh = 1.8e-4    # k_h/T_csh = temps caracteristique de carbo des CSH (1/s)\n") ;
  fprintf(ficd,"A_1 = 150         # Coef de la cinetique 1 (dm/mole/s)\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Coef de la cinetique 2 (dm/mole/s)\n") ;
  fprintf(ficd,"C_2 = 0.14e6      # Coef de la cinetique 2 (1/dm)\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;  

  return(NEQ) ;
}


void tb47(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch47(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in47(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;
  double un = 1. ;
  
  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  phii    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  c_2     = el.mat->pr[pm("C_2")] ;
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ; 
  T_csh   = el.mat->pr[pm("T_csh")] ;
  n_csh0  = el.mat->pr[pm("N_CSH")] ;
  
  double x_natot = X_Na(1);
  double x_ktot = X_K(1);
  
  /* Contenus molaires */
  for(i=0;i<2;i++) {
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    double s_l     = courbe(p_c,el.mat->cb[0]) ;
    double s_g     = un - s_l ;
    /* molarites */
    double x_co2   = X_CO2(i) ;
    double x_oh    = X_OH(i) ;
    double x_hco3  = X_HCO3(i) ;
    double x_na    = X_Na(i) ;
    double x_k     = X_K(i) ;
    double x_h2co3,x_co3 ;
    double x_h,x_h2o ;
    double x_ca ;
    double x_naoh,x_nahco3,x_naco3 ;
    double x_koh;
    double x_caoh,x_cahco3,x_caco3aq,x_caoh2aq;
    
    
    double n_co2,n_h2co3,n_hco3,n_co3 ;
    double n_oh,n_h,n_h2o ;    
    double n_ca,n_caco3csh,n_gel,n_caco3,n_caoh2,n_csh ;
    double n_na,n_naoh,n_nahco3,n_naco3 ;
    double n_k,n_koh;
    double n_caoh,n_cahco3,n_caco3aq,n_caoh2aq ;
    double phi ;
    double rho_l ;

    if(x_hco3 < 0.) { /* Equilibre de la Portlandite */
      x_oh    = concentration_oh_equilibre(x_na,x_k) ;
      /*printf("x_ohint= %e mol.l-1\n", x_oh) ;*/
      x_hco3  = k_ca/(k_co3*k_2)*x_oh ;
      X_HCO3(i) = x_hco3 ;
    }
	
	    /*Processus Iteratif d'initialisation des variables x_na, x_ka et x_oh*/
    if(x_na == 0 && x_k != 0)
    {
    double delta = 1.;
    double tolerance = 1e-3;
    double x_oh_precedent ;
    do
    {
    x_oh_precedent = concentration_oh_equilibre(x_na,x_k) ;
    x_k = concentration_k_equilibre(x_ktot,x_oh);
    x_oh = concentration_oh_equilibre(x_na,x_k) ;
    delta = fabs((x_oh_precedent - x_oh)/x_oh);
    }while (delta>tolerance) ;  
    }
    
    else if(x_na !=0 && x_k ==0)
    {
    double delta = 1.;
    double tolerance = 1e-3;
    double x_oh_precedent ;
    do
    {
    x_oh_precedent = concentration_oh_equilibre(x_na,x_k) ;
    x_na = concentration_na_equilibre(x_natot,x_oh);
    x_oh = concentration_oh_equilibre(x_na,x_k) ;
    delta = fabs((x_oh_precedent - x_oh)/x_oh);
    }while (delta>tolerance) ;
    }
    
    else if(x_na !=0 && x_k != 0)
    {
    double delta = 1.;
    double tolerance = 1e-3;
    double x_oh_precedent ;
    do
    {
    x_oh_precedent = concentration_oh_equilibre(x_na,x_k) ;
    x_na = concentration_na_equilibre(x_natot,x_oh);
    x_k = concentration_k_equilibre(x_ktot,x_oh);
    x_oh = concentration_oh_equilibre(x_na,x_k) ;
    delta = fabs((x_oh_precedent - x_oh)/x_oh);
    }while (delta>tolerance) ;
    }   
    else
    {
    x_oh = concentration_oh_equilibre(x_na,x_k) ;
    }
        
    X_Na(i) = x_na;
    X_K(i) = x_k;    

    

    if(x_co2 < 0.) { /* Equilibre de H2CO3 <-> HCO3 */
      x_co2   = k_1/k_h*x_hco3/x_oh ;
      X_CO2(i) = x_co2 ;
    }

    

    x_oh    = concentration_oh(x_hco3,x_na,x_k) ;
    X_OH(i) = x_oh ;

    x_h2co3 = k_h*x_co2 ;
    x_co3   = k_co3*x_oh*x_hco3 ;
    x_h     = k_e/x_oh ;
    x_ca    = k_ca/x_co3 ;
    x_naoh  = k_naoh/k_e*x_na*x_oh ;
    x_nahco3= k_nahco3*x_na*x_hco3 ;
    x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    x_koh   = k_koh/k_e*x_k*x_oh ;
    x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    x_cahco3  = k_cahco3*x_hco3*x_ca ;
    x_caco3aq = k_caco3*k_ca/(k_e*k_co3) ;
    x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + x_koh*v_koh + x_caoh*v_caoh + v_cahco3*x_cahco3 + x_caco3aq*v_caco3aq + x_caoh2aq*v_caoh2aq))/v_h2o ;
    

    /* masse volumique liquide */
    rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_K*x_k + M_KOH*x_koh + M_CaOH*x_caoh + M_CaHCO3*x_cahco3 + M_CaCO3aq*x_caco3aq + M_CaOH2aq*x_caoh2aq ;

    /* solides */
    n_caco3    = N_CaCO3(i) ;
    n_csh      = n_csh0 ;
    n_caoh2    = n_caoh20 ;  
    n_caco3csh = 3*(n_csh0 - n_csh) ;
    n_gel      = (n_csh0 - n_csh) ;

    /* porosite */
    phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

    /* contenus molaires */
    n_co2   = phi*s_g*x_co2 ;
    n_oh    = phi*s_l*x_oh ;
    n_h2o   = phi*s_l*x_h2o ;
    n_hco3  = phi*s_l*x_hco3 ;
    n_h2co3 = phi*s_l*x_h2co3 ;
    n_co3   = phi*s_l*x_co3 ;
    n_h     = phi*s_l*x_h ;
    n_ca    = phi*s_l*x_ca ;
    n_na    = phi*s_l*x_na ;
    n_naoh  = phi*s_l*x_naoh ;
    n_nahco3= phi*s_l*x_nahco3 ;
    n_naco3 = phi*s_l*x_naco3 ;
    n_k     = phi*s_l*x_k ;
    n_koh   = phi*s_l*x_koh ;
    
    n_caoh  = phi*s_l*x_caoh ;
    n_cahco3= phi*s_l*x_cahco3 ;
    n_caco3aq=phi*s_l*x_caco3aq ;
    n_caoh2aq=phi*s_l*x_caoh2aq ;
    
    

    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_caco3 + n_caco3csh + n_nahco3 + n_naco3 + n_cahco3 + n_caco3aq ;
    N_Ca(i) = n_ca + n_caoh2 + n_caco3 + n_caco3csh + 3*n_csh + n_caoh + n_cahco3 + n_caco3aq + n_caoh2aq ;
    N_Na(i) = n_na + n_naoh + n_nahco3 + n_naco3 ; 
    N_K(i)  = n_k + n_koh ;

    /* masse totale */
    M(i)    = phi*s_g*M_CO2*x_co2 + phi*s_l*rho_l + M_CaOH2*n_caoh2 + M_CaCO3*n_caco3 + M_CaCO3CSH*n_caco3csh + M_CSH*n_csh + M_gel*n_gel ;

    /* cinetique */
    N_k(i)  = n_hco3 + n_co3 + n_nahco3 + n_naco3 + n_caco3 + n_cahco3 + n_caco3aq ;
    XI(i)   = phi*s_l*a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

    /* densite de charge */
    N_q(i)  = z_h*x_h + z_oh*x_oh + z_ca*x_ca + z_hco3*x_hco3 + z_co3*x_co3 + z_na*x_na + z_naco3*x_naco3 + z_k*x_k + z_caoh*x_caoh + z_cahco3*x_cahco3;

    /* contenus solides */
    N_CaOH2(i) = n_caoh2 ;
    N_CSH(i)   = n_csh ;
  }

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex47(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  int    i ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phii    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  c_2     = el.mat->pr[pm("C_2")] ;
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ;
  T_csh    = el.mat->pr[pm("T_csh")] ;
  n_csh0   = el.mat->pr[pm("N_CSH")] ;


  /* Contenus molaires */
  for(i=0;i<2;i++) {
    /* pressions */
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    /* molarites */
    double x_co2   = X_CO2(i) ;
    double x_oh    = X_OH(i) ;
    double x_hco3  = X_HCO3(i) ;
    
    /* solides */
    double n_caoh2 = N_CaOH2(i) ;
    double n_caco3 = N_CaCO3(i) ;
    double n_csh   = N_CSH(i) ;
    /* saturation */
    double s_l     = courbe(p_c,el.mat->cb[0]) ;
    /* solides */
    double n_caco3csh = 3*(n_csh0 - n_csh) ;
    /* porosite */
    double phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

    {
      double av      = n_caco3/n_caoh20 ;
      double dn1sdt  = a_2*dn1_caoh2sdt(av,c_2) ;
      /*cinétique de dissolution de la portlandite en loi log*/
      DN_CaOH2SDT(i) = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ;     
      /*cinétique de dissolution de la portlandite en loi puissance*/ 
    /*  DN_CaOH2SDT(i) = - dn1sdt*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),Xp)) ; */
    }
    DN_CSHSDT(i)   =  - phi*s_l*T_csh*x_co2 ;
  }
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct47(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phii    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  c_2     = el.mat->pr[pm("C_2")] ; 
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ;
  T_csh    = el.mat->pr[pm("T_csh")] ;
  n_csh0   = el.mat->pr[pm("N_CSH")] ;
  
  /* Contenus molaires */
  for(i=0;i<2;i++) {
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    double s_l     = courbe(p_c,el.mat->cb[0]) ;
    double s_g     = 1. - s_l ;

    /* molarites */
    double x_co2   = X_CO2(i) ;
    double x_oh    = X_OH(i) ;
    double x_hco3  = X_HCO3(i) ;
    double x_na    = X_Na(i) ;
    double x_k     = X_K(i) ;

    double x_h2co3 = k_h*x_co2 ;
    double x_co3   = k_co3*x_oh*x_hco3 ;
    double x_h     = k_e/x_oh ;
    double x_ca    = k_ca/x_co3 ;
    double x_naoh  = k_naoh/k_e*x_na*x_oh ;
    double x_nahco3= k_nahco3*x_na*x_hco3 ;
    double x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    double x_koh   = k_koh/k_e*x_k*x_oh ;
    double x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    double x_cahco3  = k_cahco3*x_hco3*x_ca ;
    double x_caco3aq = k_caco3*k_ca/(k_e*k_co3);
    double x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    double x_h2o   = (1. - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + x_koh*v_koh + x_caoh*v_caoh + v_cahco3*x_cahco3 + x_caco3aq*v_caco3aq + x_caoh2aq*v_caoh2aq))/v_h2o ;

    /* masse volumique liquide */
    double rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_K*x_k + M_KOH*x_koh +  M_CaOH*x_caoh + M_CaHCO3*x_cahco3 + M_CaCO3aq*x_caco3aq + M_CaOH2aq*x_caoh2aq;

    double n_co2,n_h2co3,n_hco3,n_co3 ;
    double n_oh,n_h,n_h2o ;
    double n_ca,n_caco3csh,n_gel,n_caoh2,n_csh ;
    double n_na,n_naoh,n_nahco3,n_naco3 ;
    double n_k,n_koh ;
    double n_caoh,n_cahco3,n_caco3aq,n_caoh2aq;
    double phi ;

    /* solides */
    double n_caco3    = N_CaCO3(i) ;
    {
      /* double av = n_caco3/n_caoh20 ; */
      double av = 1. - N_CaOH2n(i)/n_caoh20 ;
      double dn1sdt = a_2*dn1_caoh2sdt(av,c_2) ;
      /*cinétique de dissolution de la portlandite en loi log*/
      double dn_caoh2sdt = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ;    
      /*cinétique de dissolution de la portlandite en loi puissance*/ 
     /* double dn_caoh2sdt = - dn1sdt*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),Xp)) ; */
      
      n_caoh2    = N_CaOH2n(i) + dt*DN_CaOH2SDT(i) ;
      n_caoh2    = N_CaOH2n(i) + dt*dn_caoh2sdt ;
      if(n_caoh2 < 0.) n_caoh2 = 0. ;
    }
    n_csh      = N_CSHn(i)  + dt*DN_CSHSDT(i) ;
    if(n_csh < 0.) n_csh = 0. ;
    n_caco3csh = 3*(n_csh0 - n_csh) ;
    n_gel      = (n_csh0 - n_csh) ;
    
    /*       if(x_oh <= 0. || x_h2o <= 0. || x_hco3 <= 0. || x_na < 0. || x_k < 0. ){*/
      if(x_oh <= 0. || x_h2o <= 0. || x_hco3 <= 0.  ){
      printf("\n\
en x    = %e\n\
x_co2   = %e\n\
x_oh    = %e\n\
x_h2o   = %e\n\
x_hco3  = %e\n\
n_caco3 = %e\n\
x_na    = %e\n\
x_k     = %e\n\
x_naoh  = %e\n\
x_nahco3= %e\n\
x_naco3 = %e\n",x[i][0],x_co2,x_oh,x_h2o,x_hco3,n_caco3,x_na,x_k,x_naoh,x_nahco3,x_naco3) ;
      return(1) ;
    }
    
    /*if(x_na<=0.){
    x_na=0;
    x_naoh=0;
    x_nahco3=0;
    x_naco3=0;
    }*/

    /* porosite */
    phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

    /* contenus molaires */
    n_co2   = phi*s_g*x_co2 ;
    n_oh    = phi*s_l*x_oh ;
    n_h2o   = phi*s_l*x_h2o ;
    n_hco3  = phi*s_l*x_hco3 ;
    n_h2co3 = phi*s_l*x_h2co3 ;
    n_co3   = phi*s_l*x_co3 ;
    n_h     = phi*s_l*x_h ;
    n_ca    = phi*s_l*x_ca ;
    n_na    = phi*s_l*x_na ;
    n_naoh  = phi*s_l*x_naoh ;
    n_nahco3= phi*s_l*x_nahco3 ;
    n_naco3 = phi*s_l*x_naco3 ;
    n_k     = phi*s_l*x_k ;
    n_koh   = phi*s_l*x_koh ;
    n_caoh  = phi*s_l*x_caoh ;
    n_cahco3= phi*s_l*x_cahco3 ;
    n_caco3aq=phi*s_l*x_caco3aq ;
    n_caoh2aq=phi*s_l*x_caoh2aq ;

    /* contenus atomiques */
    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_caco3 + n_caco3csh + n_nahco3 + n_naco3 + n_cahco3 + n_caco3aq;
    N_Ca(i) = n_ca + n_caoh2 + n_caco3 + n_caco3csh + 3*n_csh + n_caoh + n_cahco3 + n_caco3aq + n_caoh2aq;
    N_Na(i) = n_na + n_naoh + n_nahco3 + n_naco3 ;
    N_K(i)  = n_k + n_koh ; 
   
    /* masse totale */
    M(i)    = phi*s_g*M_CO2*x_co2 + phi*s_l*rho_l + M_CaOH2*n_caoh2 + M_CaCO3*n_caco3 + M_CaCO3CSH*n_caco3csh + M_CSH*n_csh + M_gel*n_gel ;

    /* cinetique */
    N_k(i)  = n_hco3 + n_co3 + n_nahco3 + n_naco3 + n_caco3 + n_cahco3 + n_caco3aq ;
    XI(i)   = phi*s_l*a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

    /* densite de charge */
    N_q(i)  = z_h*x_h + z_oh*x_oh + z_ca*x_ca + z_hco3*x_hco3 + z_co3*x_co3 + z_na*x_na + z_naco3*x_naco3 + z_k*x_k + z_caoh*x_caoh + z_cahco3*x_cahco3 ;

    /* contenus solides */
    N_CaOH2(i) = n_caoh2 ;
    N_CSH(i)   = n_csh ;
  }

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int mx47(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])

  double trf_co2,trf_ca,trf_h,trf_oh,trf_h2co3,trf_hco3,trf_co3,trf_na,trf_naoh,trf_nahco3,trf_naco3,tr,trf_k,trf_koh, trf_caoh,trf_cahco3,trf_caco3aq,trf_caoh2aq ;
  double tre_ca,tre_h,tre_oh,tre_hco3,tre_co3,tre_na,tre_naco3,tre_q,tre_k,tre_caoh,tre_cahco3 ;
  double trd_h2co3,trd_hco3,trd_co3,trd_ca,trd_oh,trd_h,trd_na,trd_naoh,trd_nahco3,trd_naco3,trd_m,trd_k,trd_koh, trd_caoh,trd_cahco3,trd_caco3aq,trd_caoh2aq ;
  double dphisdx_oh,dphisdx_hco3,dphisdn_caco3 ;
  double ds_lsdp_l,ds_lsdp_c ;
  double dx_hsdx_oh[2],dx_co3sdx_oh[2],dx_co3sdx_hco3[2],dx_casdx_oh[2],dx_casdx_hco3[2] ;
  double dx_h2co3sdx_co2[2] ;
  double dx_h2osdx_co2,dx_h2osdx_hco3,dx_h2osdx_oh,dx_h2osdx_na,dx_h2osdx_k;
  double dx_naohsdx_oh[2],dx_naohsdx_na[2],dx_nahco3sdx_na[2],dx_nahco3sdx_hco3[2],dx_naco3sdx_hco3[2],dx_naco3sdx_na[2],dx_naco3sdx_oh[2] ;
  double dx_kohsdx_k[2],dx_kohsdx_oh[2];
  double dx_caohsdx_hco3[2],dx_cahco3sdx_oh[2],dx_caoh2aqsdx_oh[2],dx_caoh2aqsdx_hco3[2];
  double dn_caoh2sdx_oh,dn_caoh2sdx_hco3,dn_caoh2sdn_caco3 ;
  double dxisdx_oh,dxisdx_co2,dxisdx_hco3 ;
  double drho_lsdx_oh,drho_lsdx_hco3,drho_lsdx_co2,drho_lsdx_na,drho_lsdx_k ;
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double zero = 0.,un = 1.,deux = 2. ;
  double c[2] ;
  
  /*
    Initialisation 
  */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees 
  */
  phii    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  c_2     = el.mat->pr[pm("C_2")] ;
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ;
  T_csh   = el.mat->pr[pm("T_csh")] ;
  n_csh0  = el.mat->pr[pm("N_CSH")] ;

  /*
    CALCUL DE volume ET DE surf 
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  /*
    termes d'accumulation
  */
  for(i=0;i<2;i++) {
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    double s_l     = courbe(p_c,el.mat->cb[0]) ;
    double s_g     = un - s_l ;

    /* molarites */
    double x_co2   = X_CO2(i) ;
    double x_oh    = X_OH(i) ;
    double x_hco3  = X_HCO3(i) ;
    double x_na    = X_Na(i) ;
    double x_k     = X_K(i) ;

    double x_h2co3 = k_h*x_co2 ;
    double x_co3   = k_co3*x_oh*x_hco3 ;
    double x_h     = k_e/x_oh ;
    double x_ca    = k_ca/x_co3 ;
    double x_naoh  = k_naoh/k_e*x_na*x_oh ;
    double x_nahco3= k_nahco3*x_na*x_hco3 ;
    double x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    double x_koh   = k_koh/k_e*x_k*x_oh ;
    double x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    double x_cahco3  = k_cahco3*x_hco3*x_ca ;
    double x_caco3aq = k_caco3*k_ca/(k_e*k_co3) ;
    double x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    double x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + x_koh*v_koh + x_caoh*v_caoh + x_cahco3*v_cahco3 + x_caco3aq*v_caco3aq + x_caoh2aq*v_caoh2aq ))/v_h2o ;

    /* masse volumique liquide */
    double rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_K*x_k + M_KOH*x_koh + M_CaOH*x_caoh + M_CaHCO3*x_cahco3 + M_CaCO3aq*x_caco3aq + M_CaOH2aq*x_caoh2aq ;

    /* solides */
    double n_caco3    = N_CaCO3(i) ;
    double n_caoh2    = N_CaOH2(i) ;
    double n_csh      = N_CSH(i) ;
    double n_caco3csh = 3*(n_csh0 - n_csh) ;

    double xi      = a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;
     
    /* porosite */
    double phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;


    /* derivees */
    ds_lsdp_c = dcourbe(p_c,el.mat->cb[0]) ;
    ds_lsdp_l = -ds_lsdp_c ;

    dx_h2co3sdx_co2[i] = k_h ;

    dx_co3sdx_oh[i]    = k_co3*x_hco3 ;
    dx_co3sdx_hco3[i]  = k_co3*x_oh ;

    dx_hsdx_oh[i]      = - x_h/x_oh ;

    dx_casdx_oh[i]     = - x_ca/x_oh ;
    dx_casdx_hco3[i]   = - x_ca/x_hco3 ;

    dx_naohsdx_na[i]   =   k_naoh/k_e*x_oh ;
    dx_naohsdx_oh[i]   =   k_naoh/k_e*x_na ;
    
    dx_nahco3sdx_na[i]   = k_nahco3*x_hco3 ;
    dx_nahco3sdx_hco3[i] = k_nahco3*x_na  ;
    
    dx_naco3sdx_na[i]  =  k_naco3/k_e*x_oh*x_hco3 ;
    dx_naco3sdx_hco3[i]=  k_naco3/k_e*x_na*x_oh ;
    dx_naco3sdx_oh[i]  =  k_naco3/k_e*x_na*x_hco3 ;
    
    dx_kohsdx_k[i]     = x_oh*k_koh/k_e ;
    dx_kohsdx_oh[i]    = x_k*k_koh/k_e ;
    
    dx_caohsdx_hco3[i] = - x_caoh/x_hco3 ;
    
    dx_cahco3sdx_oh[i] = - x_cahco3/x_oh ;

    dx_caoh2aqsdx_oh[i]= x_caoh2aq/x_oh ;
    dx_caoh2aqsdx_hco3[i]= - x_caoh2aq/x_hco3 ;

    dx_h2osdx_co2      = -(v_h2co3*dx_h2co3sdx_co2[i])/v_h2o ;
    dx_h2osdx_hco3     = -(v_hco3 + v_co3*dx_co3sdx_hco3[i] + v_ca*dx_casdx_hco3[i] + v_nahco3*dx_nahco3sdx_hco3[i] + v_naco3*dx_naco3sdx_hco3[i] + v_caoh*dx_caohsdx_hco3[i] + v_caoh2aq*dx_caoh2aqsdx_hco3[i] )/v_h2o ;
    dx_h2osdx_oh       = -(v_co3*dx_co3sdx_oh[i] + v_oh + v_h*dx_hsdx_oh[i] + v_ca*dx_casdx_oh[i] + v_naoh*dx_naohsdx_oh[i] + v_naco3*dx_naco3sdx_oh[i] + v_koh*dx_kohsdx_oh[i] + v_cahco3*dx_cahco3sdx_oh[i] + v_caoh2aq*dx_caoh2aqsdx_oh[i] )/v_h2o ;
    dx_h2osdx_na       = -(v_na + v_naoh*dx_naohsdx_na[i] + v_nahco3*dx_nahco3sdx_na[i] + v_naco3*dx_naco3sdx_na[i])/v_h2o ;
    dx_h2osdx_k        = -(v_k + v_koh*dx_kohsdx_k[i])/v_h2o ;

    drho_lsdx_co2      = M_H2CO3*dx_h2co3sdx_co2[i] + M_H2O*dx_h2osdx_co2 ;
    drho_lsdx_oh       = M_H*dx_hsdx_oh[i] + M_OH + M_H2O*dx_h2osdx_oh + M_CO3*dx_co3sdx_oh[i] + M_Ca*dx_casdx_oh[i]  + M_NaOH*dx_naohsdx_oh[i] + M_NaCO3*dx_naco3sdx_oh[i] + M_KOH*dx_kohsdx_oh[i] + M_CaHCO3*dx_cahco3sdx_oh[i] + M_CaOH2aq*dx_caoh2aqsdx_oh[i] ;
    drho_lsdx_hco3     = M_H2O*dx_h2osdx_hco3 + M_HCO3 + M_CO3*dx_co3sdx_hco3[i] + M_Ca*dx_casdx_hco3[i] + M_NaHCO3*dx_nahco3sdx_hco3[i] + M_NaCO3*dx_naco3sdx_hco3[i] + M_CaOH*dx_caohsdx_hco3[i] + M_CaOH2aq*dx_caoh2aqsdx_hco3[i] ;
    drho_lsdx_na       = M_H2O*dx_h2osdx_na + M_Na + M_NaOH*dx_naohsdx_na[i] + M_NaHCO3*dx_nahco3sdx_na[i] + M_NaCO3*dx_naco3sdx_na[i] ;
    drho_lsdx_k        = M_H2O*dx_h2osdx_k + M_K + M_KOH*dx_kohsdx_k[i] ;

    if(n_caoh2 > 0.) {
      double av   = 1. - N_CaOH2n(i)/n_caoh20 ;
      double dn1  = dt*a_2*dn1_caoh2sdt(av,c_2) ;
      /*cinétique de dissolution de la portlandite en loi log*/
      dn_caoh2sdx_oh    =  dn1/x_oh ;
      dn_caoh2sdx_hco3  = -dn1/x_hco3 ;
      dn_caoh2sdn_caco3 =  0. ;     
      /*cinétique de dissolution de la portlandite en loi puissance*/
     /* dn_caoh2sdx_oh    =  Xp*dn1*k_ca/(k_co3*k_2*x_hco3)*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),(Xp-1))) ;
      dn_caoh2sdx_hco3  = -Xp*dn1*k_ca*x_oh/(k_co3*k_2*pow(x_hco3,2))*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),(Xp-1))) ;
      dn_caoh2sdn_caco3 =  0.  ;*/
      
    } else {
      dn_caoh2sdx_oh    = 0. ;
      dn_caoh2sdx_hco3  = 0. ;
      dn_caoh2sdn_caco3 = 0. ;
    }

    dphisdx_oh         = -v_caoh2*dn_caoh2sdx_oh ;
    dphisdx_hco3       = -v_caoh2*dn_caoh2sdx_hco3 ;
    dphisdn_caco3      = -v_caco3 - v_caoh2*dn_caoh2sdn_caco3 ;

    dxisdx_co2         =  a_1*x_oh*dx_h2co3sdx_co2[i] ;
    dxisdx_oh          =  a_1*x_h2co3 ;
    dxisdx_hco3        = -a_1*k_1 ;
    
    j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(-x_co2 + x_h2co3 + x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq ) ;
    K(E_C+j,I_CO2+j)   += volume[i]*phi*(s_g + s_l*dx_h2co3sdx_co2[i]) ;
    K(E_C+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_co3sdx_oh[i] + dx_naco3sdx_oh[i] + dx_cahco3sdx_oh[i])  + dphisdx_oh*(s_g*x_co2 + s_l*(x_h2co3 + x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq))) ;
    K(E_C+j,I_HCO3+j)  += volume[i]*(phi*s_l*(un + dx_co3sdx_hco3[i] + dx_naco3sdx_hco3[i] + dx_nahco3sdx_hco3[i]) + dphisdx_hco3*(s_g*x_co2 + s_l*(x_h2co3 + x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq ))) ;
    K(E_C+j,I_CaCO3+j) += volume[i]*(un + dphisdn_caco3*(s_g*x_co2 + s_l*(x_h2co3 + x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_caco3aq))) ;
    K(E_C+j,I_Na+j)    += volume[i]*(phi*s_l*(dx_nahco3sdx_na[i] + dx_naco3sdx_na[i])) ;

    /*
      Conservation de la charge  : div(w_q) = 0
    */

    /*
      Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
    */
    K(E_mass+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(-M_CO2*x_co2 + rho_l) ;
    K(E_mass+j,I_CO2+j)   += volume[i]*phi*(s_g*M_CO2 + s_l*drho_lsdx_co2) ;
    K(E_mass+j,I_OH+j)    += volume[i]*(phi*s_l*drho_lsdx_oh + dphisdx_oh*(s_l*rho_l + M_CO2*s_g*x_co2) + M_CaOH2*dn_caoh2sdx_oh) ;
    K(E_mass+j,I_HCO3+j)  += volume[i]*(phi*s_l*drho_lsdx_hco3 + dphisdx_hco3*(s_l*rho_l + M_CO2*s_g*x_co2) + M_CaOH2*dn_caoh2sdx_hco3) ;
    K(E_mass+j,I_CaCO3+j) += volume[i]*(dphisdn_caco3*(s_l*rho_l + M_CO2*s_g*x_co2 ) + M_CaOH2*dn_caoh2sdn_caco3 + M_CaCO3) ;
    K(E_mass+j,I_Na+j)    += volume[i]*phi*s_l*drho_lsdx_na ;
    K(E_mass+j,I_K+j)     += volume[i]*phi*s_l*drho_lsdx_k ;

    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(x_ca + x_caoh + x_cahco3 + x_caco3aq + x_caoh2aq) ;
    K(E_Ca+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_casdx_oh[i] + dx_cahco3sdx_oh[i] + dx_caoh2aqsdx_oh[i] ) + dn_caoh2sdx_oh + dphisdx_oh*s_l*(x_ca + x_caoh + x_cahco3 + x_caco3aq + x_caoh2aq)) ;
    K(E_Ca+j,I_HCO3+j)  += volume[i]*(phi*s_l*(dx_casdx_hco3[i] + dx_caohsdx_hco3[i] + dx_caoh2aqsdx_oh[i]) + dn_caoh2sdx_hco3 + dphisdx_hco3*s_l*(x_ca + x_caoh + x_cahco3 + x_caco3aq + x_caoh2aq)) ;
    K(E_Ca+j,I_CaCO3+j) += volume[i]*(un + dn_caoh2sdn_caco3 + dphisdn_caco3*s_l*(x_ca + x_caoh + x_cahco3 + x_caco3aq + x_caoh2aq)) ;

    /*
      Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
    */
    K(E_k+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 - dt*xi) ;
    K(E_k+j,I_CO2+j)   += volume[i]*phi*s_l*(-dt*dxisdx_co2) ;
    K(E_k+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_co3sdx_oh[i] + dx_naco3sdx_oh[i] + dx_cahco3sdx_oh[i] - dt*dxisdx_oh) + dphisdx_oh*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 - dt*xi)) ;
    K(E_k+j,I_HCO3+j)  += volume[i]*(phi*s_l*(un + dx_co3sdx_hco3[i] + dx_nahco3sdx_hco3[i] + dx_naco3sdx_hco3[i]  - dt*dxisdx_hco3) + dphisdx_hco3*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 - dt*xi)) ;
    K(E_k+j,I_CaCO3+j) += volume[i]*(un + dphisdn_caco3*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 - dt*xi)) ;
    K(E_k+j,I_Na+j)    += volume[i]*phi*s_l*(dx_nahco3sdx_na[i] + dx_naco3sdx_na[i]) ;
    /*
      Electroneutralite : q = 0
    */
    K(E_el+j,I_OH+j)    += volume[i]*(z_oh + z_h*dx_hsdx_oh[i] + z_ca*dx_casdx_oh[i] + z_co3*dx_co3sdx_oh[i] + z_naco3*dx_naco3sdx_oh[i] + z_cahco3*dx_cahco3sdx_oh[i]) ;
    K(E_el+j,I_HCO3+j)  += volume[i]*(z_hco3 + z_co3*dx_co3sdx_hco3[i] + z_ca*dx_casdx_hco3[i] + z_naco3*dx_naco3sdx_hco3[i] + z_caoh*dx_caohsdx_hco3[i]) ;
    K(E_el+j,I_Na+j)    += volume[i]*(z_na + z_naco3*dx_naco3sdx_na[i]) ;
    K(E_el+j,I_K+j)     += volume[i]*(z_k ) ;
    

    /*
      Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
    */
    K(E_Na+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(x_na + x_naoh + x_nahco3 + x_naco3) ;
    K(E_Na+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_naohsdx_oh[i] + dx_naco3sdx_oh[i]) + dphisdx_oh*s_l*(x_na + x_naoh + x_nahco3 + x_naco3)) ;
    K(E_Na+j,I_HCO3+j)  += volume[i]*(phi*s_l*(dx_nahco3sdx_hco3[i] + dx_naco3sdx_hco3[i]) + dphisdx_hco3*s_l*(x_na + x_naoh + x_nahco3 + x_naco3) );
    K(E_Na+j,I_CaCO3+j)  += volume[i]*dphisdn_caco3*s_l*(x_na + x_naoh + x_nahco3 + x_naco3 );
    K(E_Na+j,I_Na+j)    += volume[i]*phi*s_l*(1. + dx_naohsdx_na[i] + dx_nahco3sdx_na[i] + dx_naco3sdx_na[i]) ;
    
    /*
      Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
    K(E_K+j,I_P_l+j)    += volume[i]*phi*ds_lsdp_l*(x_k + x_koh) ;
    K(E_K+j,I_OH+j)     += volume[i]*(phi*s_l*(dx_kohsdx_oh[i]) + dphisdx_oh*s_l*(x_k + x_koh)) ;
    K(E_K+j,I_HCO3+j)   += volume[i]*dphisdx_hco3*s_l*(x_k + x_koh) ;
    K(E_K+j,I_CaCO3+j)  += volume[i]*dphisdn_caco3*s_l*(x_k + x_koh) ;
    K(E_K+j,I_K+j)      += volume[i]*phi*s_l*(1. + dx_kohsdx_k[i]) ;
    

    
    
    
  }

  /* termes d'ecoulement */
  tr        = dt*surf/dx ;

  trd_h2co3 = tr*KD_H2CO3 ;
  trd_hco3  = tr*KD_HCO3 ;
  trd_co3   = tr*KD_CO3 ;
  trd_oh    = tr*KD_OH ;
  trd_h     = tr*KD_H ;
  trd_ca    = tr*KD_Ca ;
  trd_na    = tr*KD_Na ;
  trd_naoh  = tr*KD_NaOH ;
  trd_nahco3= tr*KD_NaHCO3 ;
  trd_naco3 = tr*KD_NaCO3 ;
  trd_m     = tr*KD_m ;
  trd_k     = tr*KD_K ;
  trd_koh   = tr*KD_KOH ;
  trd_caoh  = tr*KD_CaOH ;
  trd_cahco3= tr*KD_CaHCO3 ;
  trd_caco3aq= tr*KD_CaCO3aq ;
  trd_caoh2aq= tr*KD_CaOH2aq ;

  trf_co2   = tr*KF_CO2 ;
  trf_h2co3 = tr*KF_H2CO3 ;
  trf_hco3  = tr*KF_HCO3 ;
  trf_co3   = tr*KF_CO3 ;
  trf_ca    = tr*KF_Ca ;
  trf_oh    = tr*KF_OH ;
  trf_h     = tr*KF_H ;
  trf_na    = tr*KF_Na ;
  trf_naoh  = tr*KF_NaOH ;
  trf_nahco3= tr*KF_NaHCO3 ;
  trf_naco3 = tr*KF_NaCO3 ;
  trf_k     = tr*KF_K ;
  trf_koh   = tr*KF_KOH ;
  trf_caoh  = tr*KF_CaOH ;
  trf_cahco3= tr*KF_CaHCO3 ;
  trf_caco3aq= tr*KF_CaCO3aq ;
  trf_caoh2aq= tr*KF_CaOH2aq ;  

  tre_hco3  = tr*Kpsi_HCO3 ;
  tre_co3   = tr*Kpsi_CO3 ;
  tre_ca    = tr*Kpsi_Ca ;
  tre_oh    = tr*Kpsi_OH ;
  tre_h     = tr*Kpsi_H ;
  tre_na    = tr*Kpsi_Na ;
  tre_naco3 = tr*Kpsi_NaCO3 ;
  tre_k     = tr*Kpsi_K ;
  tre_caoh  = tr*Kpsi_CaOH ;
  tre_cahco3= tr*Kpsi_CaHCO3 ;
  tre_q     = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_h2co3 + trd_hco3 + trd_co3 + trd_nahco3 + trd_naco3  + trd_cahco3 + trd_caco3aq ;
  }
  K(E_C,I_P_l)          += + c[0] ;
  K(E_C,I_P_l+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_P_l)      += - c[0] ;
  K(E_C+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co2 + trf_h2co3*dx_h2co3sdx_co2[i] ;
  }
  K(E_C,I_CO2)          += + c[0] ;
  K(E_C,I_CO2+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CO2)      += - c[0] ;
  K(E_C+NEQ,I_CO2+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3 + trf_co3*dx_co3sdx_hco3[i] + trf_nahco3*dx_nahco3sdx_hco3[i] + trf_naco3*dx_naco3sdx_hco3[i] ;
  }
  K(E_C,I_HCO3)         += + c[0] ;
  K(E_C,I_HCO3+NEQ)     += - c[1] ;
  K(E_C+NEQ,I_HCO3)     += - c[0] ;
  K(E_C+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co3*dx_co3sdx_oh[i] + trf_naco3*dx_naco3sdx_oh[i] + trf_cahco3*dx_cahco3sdx_oh[i];
  }
  K(E_C,I_OH)           += + c[0] ;
  K(E_C,I_OH+NEQ)       += - c[1] ;
  K(E_C+NEQ,I_OH)       += - c[0] ;
  K(E_C+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 + tre_naco3 + tre_cahco3;
  }
  K(E_C,I_psi)          += + c[0] ;
  K(E_C,I_psi+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_psi)      += - c[0] ;
  K(E_C+NEQ,I_psi+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
  c[i] = trf_nahco3*dx_nahco3sdx_na[i] + trf_naco3*dx_naco3sdx_na[i] ;
  }
  K(E_C,I_Na)           += + c[0] ;
  K(E_C,I_Na+NEQ)       += - c[1] ;
  K(E_C+NEQ,I_Na)       += - c[0] ;
  K(E_C+NEQ,I_Na+NEQ)   += + c[1] ;

  /*
    Conservation de la charge  : div(w_q) = 0
  */
  for(i=0;i<2;i++){
    c[i] = z_hco3*trf_hco3 + z_co3*trf_co3*dx_co3sdx_hco3[i] + z_ca*trf_ca*dx_casdx_hco3[i] + z_naco3*trf_naco3*dx_naco3sdx_hco3[i] + z_caoh*trf_caoh*dx_caohsdx_hco3[i] ;
  }
  K(E_q,I_HCO3)         += + c[0] ;
  K(E_q,I_HCO3+NEQ)     += - c[1] ;
  K(E_q+NEQ,I_HCO3)     += - c[0] ;
  K(E_q+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*dx_hsdx_oh[i] + z_oh*trf_oh + z_co3*trf_co3*dx_co3sdx_oh[i] + z_ca*trf_ca*dx_casdx_oh[i] + z_naco3*trf_naco3*dx_naco3sdx_oh[i] + z_cahco3*trf_cahco3*dx_cahco3sdx_oh[i] ;
  }
  K(E_q,I_OH)           += + c[0] ;
  K(E_q,I_OH+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_OH)       += - c[0] ;
  K(E_q+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_na*trf_na + z_naco3*trf_naco3*dx_naco3sdx_na[i] ;
  }
  K(E_q,I_Na)           += + c[0] ;
  K(E_q,I_Na+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Na)       += - c[0] ;
  K(E_q+NEQ,I_Na+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = z_k*trf_k  ;
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
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  K(E_mass,I_P_l)          += + trd_m ;
  K(E_mass,I_P_l+NEQ)      += - trd_m ;
  K(E_mass+NEQ,I_P_l)      += - trd_m ;
  K(E_mass+NEQ,I_P_l+NEQ)  += + trd_m ;

  for(i=0;i<2;i++){
    c[i] = M_CO2*trf_co2 ;
  }
  K(E_mass,I_CO2)          += + c[0] ;
  K(E_mass,I_CO2+NEQ)      += - c[1] ;
  K(E_mass+NEQ,I_CO2)      += - c[0] ;
  K(E_mass+NEQ,I_CO2+NEQ)  += + c[1] ;

  
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_ca + trd_caoh + trd_cahco3 + trd_caco3aq + trd_caoh2aq;
  }
  K(E_Ca,I_P_l)          += + c[0] ;
  K(E_Ca,I_P_l+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_P_l)      += - c[0] ;
  K(E_Ca+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*dx_casdx_hco3[i] + trf_caoh*dx_caohsdx_hco3[i] + trf_caoh2aq*dx_caoh2aqsdx_hco3[i] ;
  }
  K(E_Ca,I_HCO3)         += + c[0] ;
  K(E_Ca,I_HCO3+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_HCO3)     += - c[0] ;
  K(E_Ca+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*dx_casdx_oh[i] + trf_cahco3*dx_cahco3sdx_oh[i] + trf_caoh2aq*dx_caoh2aqsdx_oh[i] ;
  }
  K(E_Ca,I_OH)           += + c[0] ;
  K(E_Ca,I_OH+NEQ)       += - c[1] ;
  K(E_Ca+NEQ,I_OH)       += - c[0] ;
  K(E_Ca+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_caoh + tre_cahco3 ;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  
  /*
    Cinetique 1 : (n_k1 - n_kn) + dt * div(W_k) - dt * XI = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_hco3 + trd_co3 + trd_nahco3 + trd_naco3 + trd_caco3aq + trd_cahco3 ;
  }
  K(E_k,I_P_l)          += + c[0] ;
  K(E_k,I_P_l+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_P_l)      += - c[0] ;
  K(E_k+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3 + trf_co3*dx_co3sdx_hco3[i] + trd_nahco3*dx_nahco3sdx_hco3[i] + trd_naco3*dx_naco3sdx_hco3[i]  ;
  }
  K(E_k,I_HCO3)          += + c[0] ;
  K(E_k,I_HCO3+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_HCO3)      += - c[0] ;
  K(E_k+NEQ,I_HCO3+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co3*dx_co3sdx_oh[i] + trd_naco3*dx_naco3sdx_oh[i] + trd_cahco3*dx_cahco3sdx_oh[i] ;
  }
  K(E_k,I_OH)            += + c[0] ;
  K(E_k,I_OH+NEQ)        += - c[1] ;
  K(E_k+NEQ,I_OH)        += - c[0] ;
  K(E_k+NEQ,I_OH+NEQ)    += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 + tre_naco3 + tre_cahco3 ;
  }
  K(E_k,I_psi)           += + c[0] ;
  K(E_k,I_psi+NEQ)       += - c[1] ;
  K(E_k+NEQ,I_psi)       += - c[0] ;
  K(E_k+NEQ,I_psi+NEQ)   += + c[1] ;

  /*
    Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_na + trd_naoh + trd_nahco3 + trd_naco3 ;
  }
  K(E_Na,I_P_l)          += + c[0] ;
  K(E_Na,I_P_l+NEQ)      += - c[1] ;
  K(E_Na+NEQ,I_P_l)      += - c[0] ;
  K(E_Na+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_nahco3*dx_nahco3sdx_hco3[i] + trf_naco3*dx_naco3sdx_hco3[i] ;
  }
  K(E_Na,I_HCO3)         += + c[0] ;
  K(E_Na,I_HCO3+NEQ)     += - c[1] ;
  K(E_Na+NEQ,I_HCO3)     += - c[0] ;
  K(E_Na+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_naoh*dx_naohsdx_oh[i] + trf_naco3*dx_naco3sdx_oh[i];
  }
  K(E_Na,I_OH)           += + c[0] ;
  K(E_Na,I_OH+NEQ)       += - c[1] ;
  K(E_Na+NEQ,I_OH)       += - c[0] ;
  K(E_Na+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_na + tre_naco3 ;
  }
  K(E_Na,I_psi)          += + c[0] ;
  K(E_Na,I_psi+NEQ)      += - c[1] ;
  K(E_Na+NEQ,I_psi)      += - c[0] ;
  K(E_Na+NEQ,I_psi+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_na + trf_naoh*dx_naohsdx_na[i] + trf_nahco3*dx_nahco3sdx_na[i] + trf_naco3*dx_naco3sdx_na[i] ;
  }
  K(E_Na,I_Na)          += + c[0] ;
  K(E_Na,I_Na+NEQ)      += - c[1] ;
  K(E_Na+NEQ,I_Na)      += - c[0] ;
  K(E_Na+NEQ,I_Na+NEQ)  += + c[1] ;
  
  /*
    Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_k + trd_koh ;
  }
  K(E_K,I_P_l)          += + c[0] ;
  K(E_K,I_P_l+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_P_l)      += - c[0] ;
  K(E_K+NEQ,I_P_l+NEQ)  += + c[1] ;


  for(i=0;i<2;i++){
    c[i] = trf_koh*dx_kohsdx_oh[i] ;
  }
  K(E_K,I_OH)           += + c[0] ;
  K(E_K,I_OH+NEQ)       += - c[1] ;
  K(E_K+NEQ,I_OH)       += - c[0] ;
  K(E_K+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_k  ;
  }
  K(E_K,I_psi)          += + c[0] ;
  K(E_K,I_psi+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_psi)      += - c[0] ;
  K(E_K+NEQ,I_psi+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_k + trf_koh*dx_kohsdx_k[i]  ;
  }
  K(E_K,I_K)          += + c[0] ;
  K(E_K,I_K+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_K)      += - c[0] ;
  K(E_K+NEQ,I_K+NEQ)  += + c[1] ;

  return(0) ;

#undef K
}


void rs47(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])

  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i=0;i<NEQ*2;i++) r[i] = zero ;

  if(el.dim < dim) return ;
  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  R(0,E_C) -= volume[0]*(N_C(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C(1) - N_Cn(1)) - dt*surf*W_C ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;
  /*
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  R(0,E_mass) -= volume[0]*(M(0) - M_n(0)) + dt*surf*W_m ;
  R(1,E_mass) -= volume[1]*(M(1) - M_n(1)) - dt*surf*W_m ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
  */
  R(0,E_k) -= volume[0]*(N_k(0) - N_kn(0) - dt*XI(0)) + dt*surf*W_k ;
  R(1,E_k) -= volume[1]*(N_k(1) - N_kn(1) - dt*XI(1)) - dt*surf*W_k ;
  /*
    Electroneutralite : q = 0
  */
  R(0,E_el) -= volume[0]*N_q(0) ;
  R(1,E_el) -= volume[1]*N_q(1) ;
  /*
    Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  R(0,E_Na) -= volume[0]*(N_Na(0) - N_Nan(0)) + dt*surf*W_Na ;
  R(1,E_Na) -= volume[1]*(N_Na(1) - N_Nan(1)) - dt*surf*W_Na ;
  /*
    Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ;

#undef R
}


int so47(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double x_co2,x_h2co3,x_hco3,x_co3 ;
  double x_oh,x_h,x_h2o ;
  double x_ca,x_q ;
  double x_na,x_naoh,x_nahco3,x_naco3 ;
  double x_k,x_koh;
  double x_caoh,x_cahco3,x_caco3aq,x_caoh2aq ;
  double grd_ca,grd_oh,grd_h,grd_h2co3,grd_hco3,grd_co3,grd_psi ;
  double w_ca,w_h,w_oh,w_hco3,w_co3,w_h2co3,w_co2,w_m ;
  double n_caco3,n_caoh2,n_csh,n_caco3csh ;
  double grd_co2,grd_p_l,xi ;
  double s_l,p_c,p_l,phi ;
  double dx ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0.,un = 1.,deux = 2. ;
  double psi;
  double grd_k,w_k,grd_koh,w_koh,w_K;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phii     = el.mat->pr[pm("porosite")] ;
  k_int    = el.mat->pr[pm("k_int")] ;
  a_1      = el.mat->pr[pm("A_1")] ;
  a_2      = el.mat->pr[pm("A_2")] ;
  c_2      = el.mat->pr[pm("C_2")] ;
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ;
  T_csh    = el.mat->pr[pm("T_csh")] ;
  n_csh0   = el.mat->pr[pm("N_CSH")] ;

  /* initialisation */
  nso = 28 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pression */
  p_l    =  param(u,h_s,el.nn,I_P_l) ;
  /* saturation */
  p_c     = p_g - p_l ;
  s_l     = courbe(p_c,el.mat->cb[0]) ;
  /* molarites */
  x_co2  =  param(u,h_s,el.nn,I_CO2) ;
  x_oh   =  param(u,h_s,el.nn,I_OH) ;
  x_hco3 =  param(u,h_s,el.nn,I_HCO3) ;
  x_na   =  param(u,h_s,el.nn,I_Na) ;
  x_k    =  param(u,h_s,el.nn,I_K) ;

  x_h2co3 = k_h*x_co2 ;
  x_co3   = k_co3*x_oh*x_hco3 ;
  x_h     = k_e/x_oh ;
  x_ca    = k_ca/x_co3 ;
  x_naoh  = k_naoh/k_e*x_na*x_oh ;
  x_nahco3= k_nahco3*x_na*x_hco3 ;
  x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
  x_koh   = k_koh/k_e*x_k*x_oh ;
    x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    x_cahco3  = k_cahco3*x_hco3*x_ca ;
    x_caco3aq = k_caco3/k_e*x_hco3*x_ca*x_oh ;
    x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + x_koh*v_koh + x_caoh*v_caoh + v_cahco3*x_cahco3 + x_caco3aq*v_caco3aq + x_caoh2aq*v_caoh2aq))/v_h2o ;
    
  /* cinetique */
  xi      = a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

  /* densite de charge */
  x_q = 0.5*(N_q(0) + N_q(1)) ;

  /* contenus solides */
  n_caco3    = (N_CaCO3(0) + N_CaCO3(1))/deux ;
  n_caoh2    = (N_CaOH2(0) + N_CaOH2(1))/deux ;
  n_csh      = (N_CSH(0)   + N_CSH(1))/2. ;
  n_caco3csh = 3*(n_csh0 - n_csh) ;

  /* porosite */
  phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;
  

  psi = (PSI(1) + PSI(0))/2 ;


  /* Transferts */
  dx        = x[1][0] - x[0][0] ;
  grd_p_l   = (P_l(1) - P_l(0))/dx ;
  grd_co2   = (X_CO2(1) - X_CO2(0))/dx ;
  grd_ca    = (k_ca/k_co3)*(1/(X_OH(1)*X_HCO3(1)) - 1/(X_OH(0)*X_HCO3(0)))/dx ;
  grd_oh    = (X_OH(1) - X_OH(0))/dx ;
  grd_h     = (k_e/X_OH(1) - k_e/X_OH(0))/dx ;
  grd_h2co3 = k_h*(X_CO2(1) - X_CO2(0))/dx ;
  grd_hco3  = (X_HCO3(1) - X_HCO3(0))/dx ;
  grd_co3   = k_co3*(X_HCO3(1)*X_OH(1) - X_HCO3(0)*X_OH(0))/dx ;
  grd_psi   = (PSI(1) - PSI(0))/dx ;
  grd_k      = (X_K(1)      - X_K(0)     )/dx ;
  grd_koh    = ( k_koh/k_e*X_K(1)*X_OH(1)    -   k_koh/k_e*X_K(0)*X_OH(0)  )/dx ;
  
  /* flux */
  w_ca    = - KD_Ca*grd_p_l    - KF_Ca*grd_ca     - Kpsi_Ca*grd_psi ;
  w_h     = - KD_H*grd_p_l     - KF_H*grd_h       - Kpsi_H*grd_psi ;
  w_oh    = - KD_OH*grd_p_l    - KF_OH*grd_oh     - Kpsi_OH*grd_psi ;
  w_hco3  = - KD_HCO3*grd_p_l  - KF_HCO3*grd_hco3 - Kpsi_HCO3*grd_psi ;
  w_co3   = - KD_CO3*grd_p_l   - KF_CO3*grd_co3   - Kpsi_CO3*grd_psi ;
  w_h2co3 = - KD_H2CO3*grd_p_l - KF_H2CO3*grd_h2co3 ;
  w_co2   =                    - KF_CO2*grd_co2 ;
  w_m     = - KD_m*grd_p_l + M_CO2*w_co2 ;
  w_k     = - KD_K*grd_p_l     - KF_K*grd_k         - Kpsi_K*grd_psi     ;
  w_koh   = - KD_KOH*grd_p_l   - KF_KOH*grd_koh                          ;
  W_K     =   w_k + w_koh ;

  /* quantites exploitees */
  i = 0 ;
  strcpy(r[i].text,"x_co2") ; r[i].n = 1 ;
  r[i++].v[0] = x_co2 ;
  strcpy(r[i].text,"ph") ; r[i].n = 1 ;
  r[i++].v[0] = 14 + log(x_oh)/log(10.) ;
  strcpy(r[i].text,"n_csh") ; r[i].n = 1 ;
  r[i++].v[0] = n_csh ;
  strcpy(r[i].text,"porosite") ; r[i].n = 1 ;
  r[i++].v[0] = phi ;
  strcpy(r[i].text,"n_caoh2") ; r[i].n = 1 ;
  r[i++].v[0] = n_caoh2 ;
  strcpy(r[i].text,"x_ca") ; r[i].n = 1 ;
  r[i++].v[0] = x_ca ;
  strcpy(r[i].text,"x_co3") ; r[i].n = 1 ;
  r[i++].v[0] = x_co3 ;
  strcpy(r[i].text,"x_hco3") ; r[i].n = 1 ;
  r[i++].v[0] = x_hco3 ;
  strcpy(r[i].text,"n_caco3") ; r[i].n = 1 ;
  r[i++].v[0] = n_caco3 ;
  strcpy(r[i].text,"x_h") ; r[i].n = 1 ;
  r[i++].v[0] = x_h ;
  strcpy(r[i].text,"x_oh") ; r[i].n = 1 ;
  r[i++].v[0] = x_oh ;
  strcpy(r[i].text,"saturation") ; r[i].n = 1 ;
  r[i++].v[0] = s_l ;
  strcpy(r[i].text,"grad_psi") ; r[i].n = 1 ;
  r[i++].v[0] = grd_psi ;
  strcpy(r[i].text,"charge") ; r[i].n = 1 ;
  r[i++].v[0] = x_q ;
  strcpy(r[i].text,"x_na") ; r[i].n = 1 ;
  r[i++].v[0] = x_na ;
  strcpy(r[i].text,"x_naoh") ; r[i].n = 1 ;
  r[i++].v[0] = x_naoh ;
  strcpy(r[i].text,"x_nahco3") ; r[i].n = 1 ;
  r[i++].v[0] = x_nahco3 ;
  strcpy(r[i].text,"x_naco3") ; r[i].n = 1 ;
  r[i++].v[0] = x_naco3 ;
  strcpy(r[i].text,"x_k") ; r[i].n = 1 ;
  r[i++].v[0] = x_k ;
  strcpy(r[i].text,"x_koh") ; r[i].n = 1 ;
  r[i++].v[0] = x_koh ;
  strcpy(r[i].text,"x_caoh") ; r[i].n = 1 ;
  r[i++].v[0] = x_caoh ;
  strcpy(r[i].text,"x_cahco3") ; r[i].n = 1 ;
  r[i++].v[0] = x_cahco3 ;  
  strcpy(r[i].text,"x_caco3aq") ; r[i].n = 1 ;
  r[i++].v[0] = x_caco3aq ;
  strcpy(r[i].text,"x_caoh2aq") ; r[i].n = 1 ;
  r[i++].v[0] = x_caoh2aq ;
  strcpy(r[i].text,"p_l") ; r[i].n = 1 ;
  r[i++].v[0] = p_l ;
  strcpy(r[i].text,"n_caco3csh") ; r[i].n = 1 ;
  r[i++].v[0] = n_caco3csh ;  
  strcpy(r[i].text,"n_ktot") ; r[i].n = 1 ;
  r[i++].v[0] = (x_k +x_koh)*phi*s_l ;
  strcpy(r[i].text,"n_natot") ; r[i].n = 1 ;
  r[i++].v[0] = (x_na +x_naoh + x_nahco3 + x_naco3)*phi*s_l ; 
  return(nso) ;
}


void transfert(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom) 
/* Termes explicites (va)  */
{
  int    i ; 
  /*
    Donnees
  */
  phii    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ;
  n_csh0   = el.mat->pr[pm("N_CSH")] ;

  /* initialisation */
  for(i=0;i<NVE-4;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<2;i++) {
    /* pressions */
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    /* molarites */
    double x_co2   = X_CO2(i) ;
    double x_oh    = X_OH(i) ;
    double x_hco3  = X_HCO3(i) ;
    double x_na    = X_Na(i) ;
    double x_k     = X_K(i) ;
    /* solides */
    double n_caoh2 = N_CaOH2(i) ;
    double n_caco3 = N_CaCO3(i) ;
    double n_csh   = N_CSH(i) ;
    /* saturation */
    double s_l     = courbe(p_c,el.mat->cb[0]) ;
    double s_g     = 1. - s_l ;

    /* autres molarites */
    double x_h2co3 = k_h*x_co2 ;
    double x_co3   = k_co3*x_oh*x_hco3 ;
    double x_h     = k_e/x_oh ;
    double x_ca    = k_ca/x_co3 ;
    double x_naoh  = k_naoh/k_e*x_na*x_oh ;
    double x_nahco3= k_nahco3*x_na*x_hco3 ;
    double x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    double x_koh   = k_koh/k_e*x_k*x_oh ;
    double x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    double x_cahco3  = k_cahco3*x_hco3*x_ca ;
    double x_caco3aq = k_caco3/k_e*x_hco3*x_ca*x_oh ;
    double x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    double x_h2o   = (1. - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + x_koh*v_koh + x_caoh*v_caoh + x_cahco3*v_cahco3 + x_caco3aq*v_caco3aq + x_caoh2aq*v_caoh2aq ))/v_h2o ;

    /* masse volumique liquide */
    double rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_K*x_k + M_KOH*x_koh + M_CaOH*x_caoh + M_CaHCO3*x_cahco3 + M_CaCO3aq*x_caco3aq + M_CaOH2aq*x_caoh2aq ;
  
    /* autres solides */
    double n_caco3csh = 3*(n_csh0 - n_csh) ;
 
    /* porosite */
    double phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

    /* permeabilite */
    double k_l  = (k_int/mu_l)*courbe(p_c,el.mat->cb[1])*pow(phi/phii,3.)*pow(((1-phii)/(1-phi)),2.) ;
    /* tortuosite gaz */
    /*double tau  = pow(phi,1/3)*pow(s_g,7/3) ;*/
    double tau  = pow(phi,1.74)*pow(s_g,3.20) ;
    /* tortuosite liquide */
    double iff    = 0.00029*exp(9.95*phi)/(1+625*pow((1-s_l),4)) ;
 
    KD_Ca     += x_ca*k_l ;
    KD_H2CO3  += x_h2co3*k_l ;
    KD_HCO3   += x_hco3*k_l ;
    KD_CO3    += x_co3*k_l ;
    KD_OH     += x_oh*k_l ;
    KD_H      += x_h*k_l ;
    KD_Na     += x_na*k_l ;
    KD_NaOH   += x_naoh*k_l ;
    KD_NaHCO3 += x_nahco3*k_l ;
    KD_NaCO3  += x_naco3*k_l ;
    KD_K      += x_k*k_l ;
    KD_KOH    += x_koh*k_l ;
    KD_CaOH   += x_caoh*k_l ;
    KD_CaHCO3 += x_cahco3*k_l ;
    KD_CaCO3aq+= x_caco3aq*k_l ;
    KD_CaOH2aq+= x_caoh2aq*k_l ;
    KD_m      += rho_l*k_l ;
    
    KF_CO2    += phi*s_g*tau*d_co2 ;
    KF_Ca     += d_ca*iff ;
    KF_OH     += d_oh*iff ;
    KF_H      += d_h*iff ;
    KF_H2CO3  += d_h2co3*iff ;
    KF_HCO3   += d_hco3*iff ;
    KF_CO3    += d_co3*iff ;
    KF_Na     += d_na*iff ;
    KF_NaOH   += d_naoh*iff ;
    KF_NaHCO3 += d_nahco3*iff ;
    KF_NaCO3  += d_naco3*iff ;
    KF_K      += d_k*iff ;
    KF_KOH    += d_koh*iff ;
    KF_CaOH   += d_caoh*iff ;
    KF_CaHCO3 += d_cahco3*iff ;
    KF_CaCO3aq+= d_caco3aq*iff ;
    KF_CaOH2aq+= d_caoh2aq*iff ;
    
    Kpsi_Ca   += FARADAY/RT*KF_Ca*z_ca*x_ca ;
    Kpsi_HCO3 += FARADAY/RT*KF_HCO3*z_hco3*x_hco3 ;
    Kpsi_CO3  += FARADAY/RT*KF_CO3*z_co3*x_co3 ;
    Kpsi_OH   += FARADAY/RT*KF_OH*z_oh*x_oh ;
    Kpsi_H    += FARADAY/RT*KF_H*z_h*x_h ;
    Kpsi_Na   += FARADAY/RT*KF_Na*z_na*x_na ;
    Kpsi_NaCO3+= FARADAY/RT*KF_NaCO3*z_naco3*x_naco3 ;
    Kpsi_K    += FARADAY/RT*KF_K*z_k*x_k ;
    Kpsi_CaOH += FARADAY/RT*KF_CaOH*z_caoh*x_caoh ;
    Kpsi_CaHCO3 +=FARADAY/RT*KF_CaHCO3*z_cahco3*x_cahco3 ;
    Kpsi_q    += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca + z_na*Kpsi_Na + z_naco3*Kpsi_NaCO3 + z_k*Kpsi_K + z_caoh*Kpsi_CaOH + z_cahco3*Kpsi_CaHCO3 ;
  }
  
  /* moyenne */
  for(i=0;i<NVE-4;i++) va[i] *= 0.5 ;
}





void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les flux (f) */
{
  double x_h[2] ;
  double x_h2co3[2],x_co3[2] ;
  double x_ca[2] ;
  double x_naoh[2],x_nahco3[2],x_naco3[2] ;
  double x_koh[2] ;
  double x_caoh[2],x_cahco3[2],x_caco3aq[2],x_caoh2aq[2] ;

  double grd_p_l ;
  double grd_oh,grd_h ;
  double grd_co2,grd_h2co3,grd_hco3,grd_co3 ;
  double grd_ca ;
  double grd_na,grd_naoh,grd_nahco3,grd_naco3 ;
  double grd_k,grd_koh ;
  double grd_caoh,grd_cahco3,grd_caco3aq,grd_caoh2aq ;
  double grd_psi ;

  double w_h,w_oh ;
  double w_co2,w_hco3,w_co3,w_h2co3 ;
  double w_ca ;
  double w_na,w_naoh,w_nahco3,w_naco3 ;
  double w_k,w_koh ;
  double w_caoh,w_cahco3,w_caco3aq,w_caoh2aq ;
  double w_m ;
  double w_q ;

  double dx = x[1][0] - x[0][0] ;
  int    i ;

  for(i=0;i<2;i++) {
    double x_co2  = X_CO2(i) ;
    double x_oh   = X_OH(i) ;
    double x_hco3 = X_HCO3(i) ;
    double x_na   = X_Na(i) ;
    double x_k    = X_K(i) ;

    x_h2co3[i]  = k_h*x_co2 ;
    x_co3[i]    = k_co3*x_oh*x_hco3 ;
    x_h[i]      = k_e/x_oh ;
    x_ca[i]     = k_ca/x_co3[i] ;
    x_naoh[i]   = k_naoh/k_e*x_na*x_oh ;
    x_nahco3[i] = k_nahco3*x_na*x_hco3 ;
    x_naco3[i]  = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    x_koh[i]    = k_koh/k_e*x_k*x_oh ;
    x_caoh[i]   = k_caoh/k_e*x_oh*x_ca[i] ;
    x_cahco3[i] = k_cahco3*x_hco3*x_ca[i] ;
    x_caco3aq[i] = k_caco3*k_ca/(k_e*k_co3) ;
    x_caoh2aq[i] = k_caoh2*x_ca[i]*x_oh*x_oh ;
  }

  /* Gradients */
  grd_p_l    = (P_l(1)      - P_l(0)     )/dx ;
  grd_co2    = (X_CO2(1)    - X_CO2(0)   )/dx ;
  grd_oh     = (X_OH(1)     - X_OH(0)    )/dx ;
  grd_hco3   = (X_HCO3(1)   - X_HCO3(0)  )/dx ;
  grd_na     = (X_Na(1)     - X_Na(0)    )/dx ;
  grd_k      = (X_K(1)      - X_K(0)     )/dx ;
  grd_ca     = (x_ca[1]     - x_ca[0]    )/dx ;
  grd_h      = (x_h[1]      - x_h[0]     )/dx ;
  grd_h2co3  = (x_h2co3[1]  - x_h2co3[0] )/dx ;
  grd_co3    = (x_co3[1]    - x_co3[0]   )/dx ;
  grd_naoh   = (x_naoh[1]   - x_naoh[0]  )/dx ;
  grd_nahco3 = (x_nahco3[1] - x_nahco3[0])/dx ;
  grd_naco3  = (x_naco3[1]  - x_naco3[0] )/dx ;
  grd_koh    = (x_koh[1]    - x_koh[0]   )/dx ;
  grd_caoh   = (x_caoh[1]   - x_caoh[0]  )/dx ;
  grd_cahco3 = (x_cahco3[1] - x_cahco3[0])/dx ;
  grd_caco3aq= (x_caco3aq[1]- x_caco3aq[0])/dx ;
  grd_caoh2aq= (x_caoh2aq[1]- x_caoh2aq[0])/dx ;
  grd_psi    = (PSI(1)      - PSI(0)     )/dx ;

  /* Flux */
  w_ca    = - KD_Ca*grd_p_l    - KF_Ca*grd_ca       - Kpsi_Ca*grd_psi    ;
  w_h     = - KD_H*grd_p_l     - KF_H*grd_h         - Kpsi_H*grd_psi     ;
  w_oh    = - KD_OH*grd_p_l    - KF_OH*grd_oh       - Kpsi_OH*grd_psi    ;
  w_hco3  = - KD_HCO3*grd_p_l  - KF_HCO3*grd_hco3   - Kpsi_HCO3*grd_psi  ;
  w_co3   = - KD_CO3*grd_p_l   - KF_CO3*grd_co3     - Kpsi_CO3*grd_psi   ; 
  w_h2co3 = - KD_H2CO3*grd_p_l - KF_H2CO3*grd_h2co3                      ;
  w_na    = - KD_Na*grd_p_l    - KF_Na*grd_na       - Kpsi_Na*grd_psi    ;
  w_naoh  = - KD_NaOH*grd_p_l  - KF_NaOH*grd_naoh                        ;
  w_nahco3= - KD_NaHCO3*grd_p_l- KF_NaHCO3*grd_nahco3                    ;
  w_naco3 = - KD_NaCO3*grd_p_l - KF_NaCO3*grd_naco3 - Kpsi_NaCO3*grd_psi ;
  w_k     = - KD_K*grd_p_l     - KF_K*grd_k         - Kpsi_K*grd_psi     ;
  w_koh   = - KD_KOH*grd_p_l   - KF_KOH*grd_koh                          ;
  w_caoh  = - KD_CaOH*grd_p_l  - KF_CaOH*grd_caoh   - Kpsi_CaOH*grd_psi  ;
  w_cahco3= - KD_CaHCO3*grd_p_l- KF_CaHCO3*grd_cahco3-Kpsi_CaHCO3*grd_psi;
  w_caco3aq=- KD_CaCO3aq*grd_p_l-KF_CaCO3aq*grd_caco3aq                  ;
  w_caoh2aq=- KD_CaOH2aq*grd_p_l-KF_CaOH2aq*grd_caoh2aq                  ;
  w_co2   =                    - KF_CO2*grd_co2                          ;
  w_m     = - KD_m*grd_p_l + M_CO2*w_co2 ;
  w_q     =                    - z_h*KF_H*grd_h             \
                               - z_oh*KF_OH*grd_oh          \
                               - z_hco3*KF_HCO3*grd_hco3    \
                               - z_co3*KF_CO3*grd_co3       \
                               - z_ca*KF_Ca*grd_ca          \
                               - z_na*KF_Na*grd_na          \
                               - z_naco3*KF_NaCO3*grd_naco3 \
                               - z_k*KF_K*grd_k             \
                               - z_caoh*KF_CaOH*grd_caoh    \
                               - z_cahco3*KF_CaHCO3*grd_cahco3 \
                                                    - Kpsi_q*grd_psi     ;


  W_C     = w_co2 + w_h2co3 + w_hco3 + w_co3 + w_nahco3 + w_naco3 + w_cahco3 + w_caco3aq ;
  W_Ca    = w_ca + w_caoh + w_cahco3 + w_caco3aq + w_caoh2aq ;
  W_Na    = w_na + w_naoh + w_nahco3 + w_naco3 ;
  W_m     = w_m ;
  W_k     = w_hco3 + w_co3 + w_nahco3 + w_naco3 + w_caco3aq + w_cahco3 ;
  W_q     = w_q ;
  W_K     = w_k + w_koh ;
}

double concentration_oh_equilibre_approx(double x_na, double x_k)
/* racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double a = z_co3*k_ca/k_2 + z_naco3*k_naco3/k_e*k_ca/(k_co3*k_2)*x_na ;
  double b = z_oh + z_hco3*k_ca/(k_co3*k_2) ;
  double c = z_na*x_na + z_k*x_k,cc = c/b ;
  double d = z_h*k_e ;
  double e = z_ca*k_2,ee = e/b ;
  /* 
     on neglige a et d
     solution de bx^3 + cx^2 + e = 0 
     que l'on met sous la forme x^3 + px + q = 0
  */
  double p = -cc*cc/3,q = 2*cc*cc*cc/27 + ee ;
  double r = sqrt(q*q/4 + p*p*p/27) ;
  double uu = pow(-q/2. + r,1/3.), vv = pow(-q/2. - r,1/3.) ;
  double x = uu + vv ; /* formule de Cardan */
  double x_oh  = x - cc/3 ;
  return(x_oh) ;
}

double concentration_oh_equilibre(double x_na,double x_k)
/*amélioration : méthode de newton*/
{
  double a = z_co3*k_ca/k_2 + z_naco3*k_naco3/k_e*k_ca/(k_co3*k_2)*x_na ;
  double b = z_oh + z_hco3*k_ca/(k_co3*k_2) ;
  double c = z_na*x_na + z_k*x_k ;
  double d = z_h*k_e + z_caoh*k_caoh*k_2/k_e + z_cahco3*k_cahco3*k_ca/k_co3 ;
  double e = z_ca*k_2 ;
  double x_oh = concentration_oh_equilibre_approx(x_na,x_k);
  double tolerance = 1e-3 ; 
  double delta =1;
  double fonction ;
  double derivee ;
  
  do
  {
  derivee = 4*a*pow(x_oh,3) + 3*b*pow(x_oh,2) + 2*c*x_oh + d ;
  fonction = a*pow(x_oh,4) + b*pow(x_oh,3) + c*pow(x_oh,2) + d*x_oh + e ;
  x_oh = x_oh -fonction/derivee ;
  fonction = a*pow(x_oh,4) + b*pow(x_oh,3) + c*pow(x_oh,2) + d*x_oh + e ;
  delta = fabs((fonction/derivee)/x_oh) ;
  printf("x_oh= %e mol.l-1\n", x_oh) ;
  printf("delta= %e\n", delta);
  }while( delta > tolerance);
  return(x_oh) ;

}

double concentration_na_equilibre(double x_na,double x_oh)
{
    double x_natotvise = x_na ;
    double x_naoh ;
    double x_natot ;
    double delta=1.;
    double tolerance=1e-3;
    
    do
    {
    x_naoh  = k_naoh/k_e*x_na*x_oh ;
    x_natot = x_na + x_naoh ;
    delta = (x_natotvise - x_natot)/x_natot;
    x_na = x_na + delta*x_na ;
    printf("x_natot= %e mol.l-1\n", x_natot) ;
    printf("delta= %e\n", delta); 
     
    }while(fabs(delta) > tolerance );
return(x_na) ;
}

double concentration_k_equilibre(double x_k,double x_oh)
{
    double x_ktotvise = x_k ;
    double x_koh ;
    double x_ktot ;
    double delta=1.;
    double tolerance=1e-3;
    
    do
    {
    x_koh  = k_koh/k_e*x_k*x_oh ; ;
    x_ktot = x_k + x_koh ;
    delta = (x_ktotvise - x_ktot)/x_ktot;
    x_k = x_k + delta*x_k ;
    printf("x_ktot= %e mol.l-1\n", x_ktot) ;
    printf("delta= %e\n", delta); 
     
    }while(fabs(delta) > tolerance );
return(x_k) ;
}



double concentration_oh(double x_hco3,double x_na, double x_k)
{
  double a = z_oh + z_co3*k_co3*x_hco3 + z_naco3*k_naco3/k_e*x_na*x_hco3 ;
  double b = (z_hco3*x_hco3 + z_na*x_na + z_k*x_k)/a ;
  double c = (z_h*k_e + z_caoh*k_caoh*k_2/k_e + z_cahco3*k_cahco3*k_ca/k_co3 + z_ca*k_ca/(k_co3*x_hco3))/a ;
  double d = b*b - 4*c ;
  return(0.5*(-b + sqrt(d))) ;
}

double dn1_caoh2sdt(double av,double c_2)
{
#define AV         ((av < 1.) ? av : 1.)
  double rp = (av < 1.) ? pow(1 - av,1./3) : 0. ;
  double rc = pow(1 - AV + v_caco3/v_caoh2*AV,1./3) ;
  /*
  double alpha2 = -1./3*av*av - 2./3*av + 1 ;
  double alpha  = -5.29478*av*av*av*av + 8.6069*av*av*av - 4.2444*av*av + 0.9325*av ;
  */

  return((rc > 0.) ? rp*rp/(1 + c_2*rp*(1 - rp/rc)) : 0.) ;
  /* return(alpha2/(1 + c_2*alpha)) ; */
#undef AV
}

