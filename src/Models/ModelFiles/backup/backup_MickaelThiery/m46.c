/*
  Trois reactions avec cinetique :
   1. acido-basique : HCO3- + H20 <-> H2CO3 + OH-    (loi de Dankwerts)
   2. dissolution de portlandite : Ca(OH)2 <-> Ca2+ + 2OH-    (loi en log)
   3. Carbo CSH  : CSH + 2H2CO3 <-> 3CaCO3 + 2SiO2.3H2O + 3H2O (loi macro)
  Cristaux de portlandite spherique, 
  Diffusion gaz MT,
  Transferts avec diffusion ionique
  Electroneutralite
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  46
#define TITLE "Carbonatation du beton (avril 2008)"
#define AUTHORS "Thiery"

#include "OldMethods.h"

/* Macros */
#define NEQ     (6)
#define NVE     (24)
#define NVI     (21)

#define E_C     (0)
#define E_q     (1)
#define E_mass  (2)
#define E_Ca    (3)
#define E_k     (4)
#define E_el    (5)

#define I_CO2   (0)
#define I_OH    (5)
#define I_P_l   (2)
#define I_CaCO3 (3)
#define I_HCO3  (4)
#define I_psi   (1)


#define X_CO2(n)   (u[(n)][I_CO2])
#define X_OH(n)    (u[(n)][I_OH])
#define N_CaCO3(n) (u[(n)][I_CaCO3])
#define P_l(n)     (u[(n)][I_P_l])
#define X_HCO3(n)  (u[(n)][I_HCO3])
#define PSI(n)     (u[(n)][I_psi])

#define N_C(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define M(n)       (f[(4+n)])
#define N_Ca(n)    (f[(6+n)])
#define N_k(n)     (f[(8+n)])
#define XI(n)      (f[(10+n)])
#define W_C        (f[12])
#define W_q        (f[13])
#define W_m        (f[14])
#define W_Ca       (f[15])
#define W_k        (f[16])
#define N_CaOH2(n) (f[(17+n)])
#define N_CSH(n)   (f[(19+n)])

#define KD_Ca      (va[(0)])
#define KD_OH      (va[(1)])
#define KD_H       (va[(2)])
#define KD_H2CO3   (va[(3)])
#define KD_HCO3    (va[(4)])
#define KD_CO3     (va[(5)])
#define KD_m       (va[(6)])

#define KF_CO2     (va[(7)])
#define KF_Ca      (va[(8)])
#define KF_OH      (va[(9)])
#define KF_H       (va[(10)])
#define KF_H2CO3   (va[(11)])
#define KF_HCO3    (va[(12)])
#define KF_CO3     (va[(13)])

#define Kpsi_Ca    (va[(14)])
#define Kpsi_OH    (va[(15)])
#define Kpsi_H     (va[(16)])
#define Kpsi_HCO3  (va[(17)])
#define Kpsi_CO3   (va[(18)])
#define Kpsi_q     (va[(19)])

#define DN_CSHSDT(n)   (va[(20+n)])
#define DN_CaOH2SDT(n) (va[(22+n)])


/* les valences */
#define z_ca    (2.)
#define z_h     (1.)
#define z_oh    (-1.)
#define z_hco3  (-1.)
#define z_co3   (-2.)

/* volumes molaires partiels des ions (dm3/mole) */
#define v_h        (-5.50e-3)
#define v_oh       (23.50e-3)
#define v_h2o      (18.e-3)
#define v_h2co3    (50.e-3)
#define v_hco3     (50.e-3)
#define v_co3      (-2.3e-3)
#define v_ca       (-18.7e-3)

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

#define M_CO2      (44.)

#define M_CaOH2    (74.)
#define M_CaCO3    (100.)
#define M_CaCO3CSH (100.)
#define M_CSH      (342.) /* 3(CaO).2(SiO2).3(H2O) */
#define M_gel      (174.) /* 2(SiO2).3(H2O) */

/* coefficients de diffusion moleculaire (dm2/s) */
#define d_oh       (5.273e-7)    /* (5.273e-7) d'apres TQN */
#define d_h        (9.310e-7)
#define d_ca       (7.92e-8)
#define d_h2co3    (7.2e-8)
#define d_hco3     (11.8e-8)
#define d_co3      (9.55e-8)
#define d_co2      (1.6e-3)

/* constantes d'equilibre (ref = 1 mole/L) */
#define k_e        (1.e-14)                /* autoprotolyse de l'eau */
#define k_h        (1.)                    /* cste de Henry */
#define k_co3      (4.570881896148751e3)   /* Equilibre de HCO3  <-> CO3 */
#define k_ca       (3.890451449942805e-9)  /* Equilibre de CaCO3 */
#define k_1        (2.187761623949552e-8)  /* Equilibre de H2CO3 <-> HCO3 */
#define k_2        (6.456542290346550e-6)  /* Equilibre de Ca(OH)2 */

/* viscosite (Pa.s) */
#define mu_l       (1.e-3)

/* constantes physiques */
#define FARADAY   (9.64846e4) /* Faraday (C/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143e3 et T=293 (Pa.dm3/mole) */


/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double) ;
static double concentration_oh_equilibre(void) ;
static double dn1_caoh2sdt(double,double) ;
/* static double ddn1_caoh2sdt(double,double) ; */
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
/* Parametres */
static double phii,k_int,a_1,a_2,c_2,n_caoh20,n_csh0,T_csh ;
/*
static double mu_l ;
static double k_e,k_h,k_ca,k_co3,k_1,k_2 ;
static double d_co2,d_ca,d_oh,d_h,d_h2co3,d_hco3,d_co3 ;
static double v_ca,v_h2o,v_hco3,v_h,v_oh,v_h2co3,v_co3 ;
static double v_csh,v_caoh2,v_caco3 ;
*/
static double p_g = 0. ;

int pm(char *s)
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

int dm46(int dim,mate_t *mat,FILE *ficd)
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

  strcpy(mat->inc[I_CO2],  "c_co2") ;
  strcpy(mat->inc[I_OH],   "c_oh") ;
  strcpy(mat->inc[I_P_l],  "p_l") ;
  strcpy(mat->inc[I_CaCO3],"c_caco3") ;
  strcpy(mat->inc[I_HCO3], "c_hco3") ;
  strcpy(mat->inc[I_psi],  "psi") ;

  {
    /* Initialisation automatique a partir des donnees 
       tirees de la these de Mickael Thiery */
    double h     = 5.6e-6 ;  /* (moles/dm2/s) these MT p 223 */
    double D     = 1.5e-15 ; /* (moles/dm/s) these MT p 253 */
    double t_csh = 3000 ;    /* (s) these MT p 230 */
    double R_0   = 110.e-5 ; /* rayon du cristal de CaOH2 (dm) */
    double a_1   = 175 ;     /* (dm/mole/s) these MT p 253 */
    double c_2   = h/D ;     /* (1/dm) these MT p 228 */

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

int qm46(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 6 equations:\n\
\t- la conservation de la masse de C   (c_co2)\n\
\t- la conservation de la charge       (psi)\n\
\t- la conservation de la masse totale (p_l)\n\
\t- la conservation de la masse de Ca  (c_caco3)\n\
\t- 1 equation de cinetique            (c_hco3)\n\
\t- Electroneutralite                  (c_oh)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n\
\t pression : Pa !\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"k_int = 5e-19     # Permeabilite intrinseque\n") ;
  fprintf(ficd,"N_CaOH2 = 6.1     # Contenu initial en moles de Ca(OH)2\n") ;
  fprintf(ficd,"N_CSH = 2.4       # contenu initial en moles de CSH\n") ;
  fprintf(ficd,"A_1 = 150         # Coef de la cinetique 1\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Coef de la cinetique 2\n") ;
  fprintf(ficd,"C_2 = 0.14e6      # Coef de la cinetique 2\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;  

  return(NEQ) ;
}


void tb46(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch46(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in46(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{

  double n_co2,n_h2co3,n_hco3,n_co3,n_oh,n_h,n_h2o,n_ca ;
  double n_caco3csh,n_gel,n_caco3,n_caoh2,n_csh ;
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double phi ;
  double s_l,s_g,p_c,p_l ;
  double rho_l ;
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

  
  /* Contenus molaires */
  for(i=0;i<2;i++) {
    p_l     = P_l(i) ;
    p_c     = p_g - p_l ;
    s_l     = courbe(p_c,el.mat->cb[0]) ;
    s_g     = un - s_l ;
    
    /* molarites */
    x_hco3  = X_HCO3(i) ;
    x_oh    = X_OH(i) ;
    x_co2   = X_CO2(i) ;

    if(x_hco3 < 0.) { /* Equilibre de la Portlandite */
      x_oh    = concentration_oh_equilibre() ;
      x_hco3  = X_HCO3(i) = k_ca/(k_co3*k_2)*x_oh ;
    }

    if(x_co2 < 0.) { /* Equilibre de H2CO3 <-> HCO3 */
      x_co2   = X_CO2(i) = k_1/k_h*x_hco3/x_oh ;
    }

    x_oh    = X_OH(i) = concentration_oh(x_hco3) ;

    x_h2co3 = k_h*x_co2 ;
    x_co3   = k_co3*x_oh*x_hco3 ;
    x_h     = k_e/x_oh ;
    x_ca    = k_ca/x_co3 ;
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;

    /* masse volumique liquide */
    rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca ;

    /* solides */
    n_caco3    = N_CaCO3(i) ;
    n_csh      = n_csh0 ;
    n_caoh2    = n_caoh20 ;  
    n_caco3csh = 3*(n_csh0 - n_csh) ;
    n_gel      = (n_csh0 - n_csh) ;

    /* porosite */
    phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;
    /* phi     = phii - dv_caoh2*n_caco3 - dv_csh*(n_csh0 - n_csh) ; */

    /* contenus molaires */
    n_co2   = phi*s_g*x_co2 ;
    n_oh    = phi*s_l*x_oh ;
    n_h2o   = phi*s_l*x_h2o ;
    n_hco3  = phi*s_l*x_hco3 ;
    n_h2co3 = phi*s_l*x_h2co3 ;
    n_co3   = phi*s_l*x_co3 ;
    n_h     = phi*s_l*x_h ;
    n_ca    = phi*s_l*x_ca ;

    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_caco3 + n_caco3csh ;
    N_Ca(i) = n_ca + n_caoh2 + n_caco3 + n_caco3csh + 3*n_csh ;

    /* masse totale */
    M(i)    = phi*s_g*M_CO2*x_co2 + phi*s_l*rho_l + M_CaOH2*n_caoh2 + M_CaCO3*n_caco3 + M_CaCO3CSH*n_caco3csh + M_CSH*n_csh + M_gel*n_gel ;

    /* cinetique */
    N_k(i)  = n_hco3 + n_co3 + n_caco3 ;
    XI(i)   = phi*s_l*a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

    /* densite de charge */
    N_q(i)  = z_h*x_h + z_oh*x_oh + z_ca*x_ca + z_hco3*x_hco3 + z_co3*x_co3 ;

    /* contenus solides */
    N_CaOH2(i) = n_caoh2 ;
    N_CSH(i)   = n_csh ;
  }

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex46(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  double x_co2,x_hco3,x_oh ;
  double n_caoh2,n_caco3,n_csh,n_caco3csh ;
  double s_l,p_c,p_l,phi ;
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
    p_l     = P_l(i) ;
    p_c     = p_g - p_l ;
    s_l     = courbe(p_c,el.mat->cb[0]) ;

    /* molarites */
    x_co2   = X_CO2(i) ;
    x_oh    = X_OH(i) ;
    x_hco3  = X_HCO3(i) ;

    /* solides */
    n_caoh2    = N_CaOH2(i) ;
    n_caco3    = N_CaCO3(i) ;
    n_csh      = N_CSH(i) ;
    n_caco3csh = 3*(n_csh0 - n_csh) ;

    /* porosite */
    phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

    {
      double av      = n_caco3/n_caoh20 ;
      double dn1sdt  = a_2*dn1_caoh2sdt(av,c_2) ;
      DN_CaOH2SDT(i) = dn1sdt*log(k_ca*x_oh/(k_co3*k_2*x_hco3)) ;
    }
    DN_CSHSDT(i)   =  - phi*s_l*T_csh*x_co2 ;
  }
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct46(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define N_CaOH2n(n) (f_n[(17+n)])
#define N_CSHn(n)   (f_n[(19+n)])

  double n_co2,n_h2co3,n_hco3,n_co3,n_oh,n_h,n_h2o,n_ca ;
  double n_caco3,n_gel,n_caco3csh,n_caoh2,n_csh ;
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double s_l,s_g,p_c,p_l,phi ;
  double rho_l ;
  int    i ;
  double un = 1. ;
  
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
    p_l     = P_l(i) ;
    p_c     = p_g - p_l ;
    s_l     = courbe(p_c,el.mat->cb[0]) ;
    s_g     = un - s_l ;

    /* molarites */
    x_co2   = X_CO2(i) ;
    x_oh    = X_OH(i) ;
    x_hco3  = X_HCO3(i) ;

    x_h2co3 = k_h*x_co2 ;
    x_co3   = k_co3*x_oh*x_hco3 ;
    x_h     = k_e/x_oh ;
    x_ca    = k_ca/x_co3 ;
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;

    /* masse volumique liquide */
    rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca ;

    /* solides */
    n_caco3    = N_CaCO3(i) ;
    {
      double av  = 1. - N_CaOH2n(i)/n_caoh20 ;
      double dn1 = dt*a_2*dn1_caoh2sdt(av,c_2) ;

      n_caoh2    = N_CaOH2n(i) + dt*DN_CaOH2SDT(i) ;
      n_caoh2    = N_CaOH2n(i) + dn1*log(k_ca*x_oh/(k_co3*k_2*x_hco3)) ;
      if(n_caoh2 < 0.) n_caoh2 = 0. ;
    }
    n_csh      = N_CSHn(i)  + dt*DN_CSHSDT(i) ;
    if(n_csh < 0.) n_csh = 0. ;
    n_caco3csh = 3*(n_csh0 - n_csh) ;
    n_gel      = (n_csh0 - n_csh) ;
    
    /* if(x_co2 <= 0. || x_oh <= 0. || x_h2o <= 0. || x_hco3 <= 0.) { */
    if(x_oh <= 0. || x_h2o <= 0. || x_hco3 <= 0.) {
      printf("\n\
en x    = %e\n\
x_co2   = %e\n\
x_oh    = %e\n\
x_h2o   = %e\n\
x_hco3  = %e\n\
n_caco3 = %e\n",x[i][0],x_co2,x_oh,x_h2o,x_hco3,n_caco3) ;
      return(1) ;
    }

    /* porosite */
    phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

    /* contenus molaires */
    n_co2   = phi*s_g*x_co2 ;
    n_h     = phi*s_l*x_h ;
    n_oh    = phi*s_l*x_oh ;
    n_h2o   = phi*s_l*x_h2o ;
    n_hco3  = phi*s_l*x_hco3 ;
    n_h2co3 = phi*s_l*x_h2co3 ;
    n_co3   = phi*s_l*x_co3 ;
    n_ca    = phi*s_l*x_ca ;

    /* contenus atomiques */
    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_caco3 + n_caco3csh ;
    N_Ca(i) = n_ca + n_caoh2 + n_caco3 + n_caco3csh + 3*n_csh ;

    /* masse totale */
    M(i)    = phi*s_g*M_CO2*x_co2 + phi*s_l*rho_l + M_CaOH2*n_caoh2 + M_CaCO3*n_caco3 + M_CaCO3CSH*n_caco3csh + M_CSH*n_csh + M_gel*n_gel ;

    /* cinetique */
    N_k(i)  = n_hco3 + n_co3 + n_caco3 ;
    XI(i)   = phi*s_l*a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

    /* densite de charge */
    N_q(i)     = z_h*x_h + z_oh*x_oh + z_ca*x_ca + z_hco3*x_hco3 + z_co3*x_co3 ;
    /* contenus solides */
    N_CaOH2(i) = n_caoh2 ;
    N_CSH(i)   = n_csh ;
  }

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;

  return(0) ;

#undef N_CaOH2n
#undef N_CSHn
}


int mx46(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define N_CaOH2n(n) (f_n[(17+n)])

#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double n_caco3,n_caoh2,n_csh,n_caco3csh ;
  double s_l,s_g,p_c,p_l,xi ;
  double trf_co2,trf_ca,trf_h,trf_oh,trf_h2co3,trf_hco3,trf_co3,tr ;
  double tre_ca,tre_h,tre_oh,tre_hco3,tre_co3,tre_q ;
  double trd_h2co3,trd_hco3,trd_co3,trd_ca,trd_oh,trd_h,trd_m ;
  double phi ;
  double dphisdx_oh,dphisdx_hco3,dphisdn_caco3 ;
  double ds_lsdp_l,ds_lsdp_c ;
  double dx_hsdx_oh[2],dx_co3sdx_oh[2],dx_co3sdx_hco3[2],dx_casdx_oh[2],dx_casdx_hco3[2] ;
  double dx_h2co3sdx_co2[2] ;
  double dx_h2osdx_co2,dx_h2osdx_hco3,dx_h2osdx_oh ;
  double dn_caoh2sdx_oh,dn_caoh2sdx_hco3,dn_caoh2sdn_caco3 ;
  double dxisdx_oh,dxisdx_co2,dxisdx_hco3 ;
  double rho_l ;
  double drho_lsdx_oh,drho_lsdx_hco3,drho_lsdx_co2 ;
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
    j = i*NEQ ;

    p_l     = P_l(i) ;
    p_c     = p_g - p_l ;
    s_l     = courbe(p_c,el.mat->cb[0]) ;
    s_g     = un - s_l ;

    /* molarites */
    x_co2   = X_CO2(i) ;
    x_oh    = X_OH(i) ;
    x_hco3  = X_HCO3(i) ;

    x_h2co3 = k_h*x_co2 ;
    x_co3   = k_co3*x_oh*x_hco3 ;
    x_h     = k_e/x_oh ;
    x_ca    = k_ca/x_co3 ;
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;

    /* masse volumique liquide */
    rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca ;

    /* solides */
    n_caco3    = N_CaCO3(i) ;
    n_caoh2    = N_CaOH2(i) ;
    n_csh      = N_CSH(i) ;
    n_caco3csh = 3*(n_csh0 - n_csh) ;

    xi      = a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;
     
    /* porosite */
    phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;


    /* derivees */
    ds_lsdp_c = dcourbe(p_c,el.mat->cb[0]) ;
    ds_lsdp_l = -ds_lsdp_c ;

    dx_h2co3sdx_co2[i] = k_h ;

    dx_co3sdx_oh[i]    = k_co3*x_hco3 ;
    dx_co3sdx_hco3[i]  = k_co3*x_oh ;

    dx_hsdx_oh[i]      = - x_h/x_oh ;

    dx_casdx_oh[i]     = - x_ca/x_oh ;
    dx_casdx_hco3[i]   = - x_ca/x_hco3 ;

    dx_h2osdx_co2      = -(v_h2co3*dx_h2co3sdx_co2[i])/v_h2o ;
    dx_h2osdx_hco3     = -(v_hco3 + v_co3*dx_co3sdx_hco3[i] + v_ca*dx_casdx_hco3[i])/v_h2o ;
    dx_h2osdx_oh       = -(v_co3*dx_co3sdx_oh[i] + v_oh + v_h*dx_hsdx_oh[i] + v_ca*dx_casdx_oh[i])/v_h2o ;

    drho_lsdx_co2      = M_H2CO3*dx_h2co3sdx_co2[i] + M_H2O*dx_h2osdx_co2 ;
    drho_lsdx_oh       = M_H*dx_hsdx_oh[i] + M_OH + M_H2O*dx_h2osdx_oh + M_CO3*dx_co3sdx_oh[i] + M_Ca*dx_casdx_oh[i] ;
    drho_lsdx_hco3     = M_H2O*dx_h2osdx_hco3 + M_HCO3 + M_CO3*dx_co3sdx_hco3[i] + M_Ca*dx_casdx_hco3[i] ;

    if(n_caoh2 > 0.) {
      double av   = 1. - N_CaOH2n(i)/n_caoh20 ;
      double dn1  = dt*a_2*dn1_caoh2sdt(av,c_2) ;
      /* double ddn1 = dt*a_2*ddn1_caoh2sdt(av,c_2)/n_caoh20 ; */

      dn_caoh2sdx_oh    =  dn1/x_oh ;
      dn_caoh2sdx_hco3  = -dn1/x_hco3 ;
      /* dn_caoh2sdn_caco3 =  ddn1*log(k_ca*x_oh/(k_co3*k_2*x_hco3)) ; */
      dn_caoh2sdn_caco3 =  0. ;
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
    
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(-x_co2 + x_h2co3 + x_hco3 + x_co3) ;
    K(E_C+j,I_CO2+j)   += volume[i]*phi*(s_g + s_l*dx_h2co3sdx_co2[i]) ;
    K(E_C+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_co3sdx_oh[i]) + dphisdx_oh*(s_g*x_co2 + s_l*(x_h2co3 + x_hco3 + x_co3))) ;
    K(E_C+j,I_HCO3+j)  += volume[i]*(phi*s_l*(un + dx_co3sdx_hco3[i]) + dphisdx_hco3*(s_g*x_co2 + s_l*(x_h2co3 + x_hco3 + x_co3))) ;
    K(E_C+j,I_CaCO3+j) += volume[i]*(un + dphisdn_caco3*(s_g*x_co2 + s_l*(x_h2co3 + x_hco3 + x_co3))) ;

    /*
      Conservation de la charge  : div(w_q) = 0
    */

    /*
      Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
    */
    K(E_mass+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(-M_CO2*x_co2 + rho_l) ;
    K(E_mass+j,I_CO2+j)   += volume[i]*phi*(s_g*M_CO2 + s_l*drho_lsdx_co2) ;
    K(E_mass+j,I_OH+j)    += volume[i]*(phi*s_l*drho_lsdx_oh + dphisdx_oh*s_l*rho_l + M_CaOH2*dn_caoh2sdx_oh) ;
    K(E_mass+j,I_HCO3+j)  += volume[i]*(phi*s_l*drho_lsdx_hco3 + dphisdx_hco3*s_l*rho_l + M_CaOH2*dn_caoh2sdx_hco3) ;
    K(E_mass+j,I_CaCO3+j) += volume[i]*(dphisdn_caco3*s_l*rho_l + M_CaOH2*dn_caoh2sdn_caco3 + M_CaCO3) ;

    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(x_ca) ;
    K(E_Ca+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_casdx_oh[i]) + dn_caoh2sdx_oh + dphisdx_oh*s_l*(x_ca)) ;
    K(E_Ca+j,I_HCO3+j)  += volume[i]*(phi*s_l*(dx_casdx_hco3[i]) + dn_caoh2sdx_hco3 + dphisdx_hco3*s_l*(x_ca)) ;
    K(E_Ca+j,I_CaCO3+j) += volume[i]*(un + dn_caoh2sdn_caco3 + dphisdn_caco3*s_l*(x_ca)) ;

    /*
      Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
    */
    K(E_k+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(x_hco3 + x_co3 - dt*xi) ;
    K(E_k+j,I_CO2+j)   += volume[i]*phi*s_l*(-dt*dxisdx_co2) ;
    K(E_k+j,I_OH+j)    += volume[i]*(phi*s_l*(dx_co3sdx_oh[i] - dt*dxisdx_oh) + dphisdx_oh*s_l*(x_co3 + x_hco3 - dt*xi)) ;
    K(E_k+j,I_HCO3+j)  += volume[i]*(phi*s_l*(un + dx_co3sdx_hco3[i] - dt*dxisdx_hco3) + dphisdx_hco3*s_l*(x_co3 + x_hco3 - dt*xi)) ;
    K(E_k+j,I_CaCO3+j) += volume[i]*(un + dphisdn_caco3*s_l*(x_co3 + x_hco3 - dt*xi)) ;

    /*
      Electroneutralite : q = 0
    */
    K(E_el+j,I_OH+j)    += volume[i]*(z_oh + z_h*dx_hsdx_oh[i] + z_ca*dx_casdx_oh[i] + z_co3*dx_co3sdx_oh[i]) ;
    K(E_el+j,I_HCO3+j)  += volume[i]*(z_hco3 + z_co3*dx_co3sdx_hco3[i] + z_ca*dx_casdx_hco3[i]) ;
  }

  /* termes d'ecoulement */
  tr        = dt*surf/dx ;

  trd_h2co3 = tr*KD_H2CO3 ;
  trd_hco3  = tr*KD_HCO3 ;
  trd_co3   = tr*KD_CO3 ;
  trd_oh    = tr*KD_OH ;
  trd_h     = tr*KD_H ;
  trd_ca    = tr*KD_Ca ;
  trd_m     = tr*KD_m ;

  trf_co2   = tr*KF_CO2 ;
  trf_h2co3 = tr*KF_H2CO3 ;
  trf_hco3  = tr*KF_HCO3 ;
  trf_co3   = tr*KF_CO3 ;
  trf_ca    = tr*KF_Ca ;
  trf_oh    = tr*KF_OH ;
  trf_h     = tr*KF_H ;

  tre_hco3  = tr*Kpsi_HCO3 ;
  tre_co3   = tr*Kpsi_CO3 ;
  tre_ca    = tr*Kpsi_Ca ;
  tre_oh    = tr*Kpsi_OH ;
  tre_h     = tr*Kpsi_H ;
  tre_q     = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_h2co3 + trd_hco3 + trd_co3 ;
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
    c[i] = trf_hco3 + trf_co3*dx_co3sdx_hco3[i] ;
  }
  K(E_C,I_HCO3)         += + c[0] ;
  K(E_C,I_HCO3+NEQ)     += - c[1] ;
  K(E_C+NEQ,I_HCO3)     += - c[0] ;
  K(E_C+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co3*dx_co3sdx_oh[i] ;
  }
  K(E_C,I_OH)           += + c[0] ;
  K(E_C,I_OH+NEQ)       += - c[1] ;
  K(E_C+NEQ,I_OH)       += - c[0] ;
  K(E_C+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 ;
  }
  K(E_C,I_psi)          += + c[0] ;
  K(E_C,I_psi+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_psi)      += - c[0] ;
  K(E_C+NEQ,I_psi+NEQ)  += + c[1] ;

  /*
    Conservation de la charge  : div(w_q) = 0
  */
  for(i=0;i<2;i++){
    c[i] = z_hco3*trf_hco3 + z_co3*trf_co3*dx_co3sdx_hco3[i] + z_ca*trf_ca*dx_casdx_hco3[i] ;
  }
  K(E_q,I_HCO3)         += + c[0] ;
  K(E_q,I_HCO3+NEQ)     += - c[1] ;
  K(E_q+NEQ,I_HCO3)     += - c[0] ;
  K(E_q+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*dx_hsdx_oh[i] + z_oh*trf_oh + z_co3*trf_co3*dx_co3sdx_oh[i] + z_ca*trf_ca*dx_casdx_oh[i] ;
  }
  K(E_q,I_OH)           += + c[0] ;
  K(E_q,I_OH+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_OH)       += - c[0] ;
  K(E_q+NEQ,I_OH+NEQ)   += + c[1] ;

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
    c[i] = trd_ca ;
  }
  K(E_Ca,I_P_l)          += + c[0] ;
  K(E_Ca,I_P_l+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_P_l)      += - c[0] ;
  K(E_Ca+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*dx_casdx_hco3[i] ;
  }
  K(E_Ca,I_HCO3)         += + c[0] ;
  K(E_Ca,I_HCO3+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_HCO3)     += - c[0] ;
  K(E_Ca+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*dx_casdx_oh[i] ;
  }
  K(E_Ca,I_OH)           += + c[0] ;
  K(E_Ca,I_OH+NEQ)       += - c[1] ;
  K(E_Ca+NEQ,I_OH)       += - c[0] ;
  K(E_Ca+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_ca ;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  
  /*
    Cinetique 1 : (n_k1 - n_kn) + dt * div(W_k) - dt * XI = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_hco3 + trd_co3 ;
  }
  K(E_k,I_P_l)          += + c[0] ;
  K(E_k,I_P_l+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_P_l)      += - c[0] ;
  K(E_k+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3 + trf_co3*dx_co3sdx_hco3[i] ;
  }
  K(E_k,I_HCO3)          += + c[0] ;
  K(E_k,I_HCO3+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_HCO3)      += - c[0] ;
  K(E_k+NEQ,I_HCO3+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co3*dx_co3sdx_oh[i] ;
  }
  K(E_k,I_OH)            += + c[0] ;
  K(E_k,I_OH+NEQ)        += - c[1] ;
  K(E_k+NEQ,I_OH)        += - c[0] ;
  K(E_k+NEQ,I_OH+NEQ)    += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 ;
  }
  K(E_k,I_psi)           += + c[0] ;
  K(E_k,I_psi+NEQ)       += - c[1] ;
  K(E_k+NEQ,I_psi)       += - c[0] ;
  K(E_k+NEQ,I_psi+NEQ)   += + c[1] ;

  return(0) ;

#undef N_CaOH2n

#undef K
}


void rs46(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define M_n(n)     (f_n[(4+n)])
#define N_Can(n)   (f_n[(6+n)])
#define N_kn(n)    (f_n[(8+n)])

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

#undef N_Cn
#undef N_qn
#undef M_n
#undef N_Can
#undef N_kn

#undef R
}


int so46(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double x_co2,x_h2co3,x_hco3,x_co3 ;
  double x_oh,x_h,x_h2o ;
  double x_ca,x_q ;
  double grd_ca,grd_oh,grd_h,grd_h2co3,grd_hco3,grd_co3,grd_psi ;
  double w_ca,w_h,w_oh,w_hco3,w_co3,w_h2co3,w_co2,w_m ;
  double n_caco3,n_caoh2,n_csh,n_caco3csh ;
  double grd_co2,grd_p_l,xi ;
  double s_l,p_c,p_l,phi ;
  double dx ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0.,un = 1.,deux = 2. ;

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
  nso = 14 ;
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

  x_h2co3 = k_h*x_co2 ;
  x_co3   = k_co3*x_oh*x_hco3 ;
  x_h     = k_e/x_oh ;
  x_ca    = k_ca/x_co3 ;
  x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;

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
  
  /* flux */
  w_ca    = - KD_Ca*grd_p_l    - KF_Ca*grd_ca     - Kpsi_Ca*grd_psi ;
  w_h     = - KD_H*grd_p_l     - KF_H*grd_h       - Kpsi_H*grd_psi ;
  w_oh    = - KD_OH*grd_p_l    - KF_OH*grd_oh     - Kpsi_OH*grd_psi ;
  w_hco3  = - KD_HCO3*grd_p_l  - KF_HCO3*grd_hco3 - Kpsi_HCO3*grd_psi ;
  w_co3   = - KD_CO3*grd_p_l   - KF_CO3*grd_co3   - Kpsi_CO3*grd_psi ;
  w_h2co3 = - KD_H2CO3*grd_p_l - KF_H2CO3*grd_h2co3 ;
  w_co2   =                    - KF_CO2*grd_co2 ;
  w_m     = - KD_m*grd_p_l ;

  /* quantites exploitees */
  strcpy(r[0].text,"x_co2") ; r[0].n = 1 ;
  r[0].v[0] = x_co2 ;
  strcpy(r[1].text,"ph") ; r[1].n = 1 ;
  r[1].v[0] = 14 + log(x_oh)/log(10.) ;
  strcpy(r[2].text,"n_csh") ; r[2].n = 1 ;
  r[2].v[0] = n_csh ;
  strcpy(r[3].text,"porosite") ; r[3].n = 1 ;
  r[3].v[0] = phi ;
  strcpy(r[4].text,"n_caoh2") ; r[4].n = 1 ;
  r[4].v[0] = n_caoh2 ;
  strcpy(r[5].text,"x_ca") ; r[5].n = 1 ;
  r[5].v[0] = x_ca ;
  strcpy(r[6].text,"x_co3") ; r[6].n = 1 ;
  r[6].v[0] = x_co3 ;
  strcpy(r[7].text,"x_hco3") ; r[7].n = 1 ;
  r[7].v[0] = x_hco3 ;
  strcpy(r[8].text,"n_caco3") ; r[8].n = 1 ;
  r[8].v[0] = n_caco3 ;
  strcpy(r[9].text,"x_h") ; r[9].n = 1 ;
  r[9].v[0] = x_h ;
  strcpy(r[10].text,"x_oh") ; r[10].n = 1 ;
  r[10].v[0] = x_oh ;
  strcpy(r[11].text,"saturation") ; r[11].n = 1 ;
  r[11].v[0] = s_l ;
  strcpy(r[12].text,"grad_psi") ; r[12].n = 1 ;
  r[12].v[0] = grd_psi ;
  strcpy(r[13].text,"charge") ; r[13].n = 1 ;
  r[13].v[0] = x_q ;
  return(nso) ;
}


void transfert(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom) 
/* Termes explicites (va)  */
{
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double n_caoh2,n_caco3,n_csh,n_caco3csh ;
  double s_l,s_g,p_c,p_l,k_l,tau,phi,iff ;
  double rho_l ;
  double un = 1.,deux = 2. ;
  
  /*
    Donnees
  */
  phii     = el.mat->pr[pm("porosite")] ;
  k_int    = el.mat->pr[pm("k_int")] ;
  n_caoh20 = el.mat->pr[pm("N_CaOH2")] ;
  n_csh0   = el.mat->pr[pm("N_CSH")] ;


  /*
    Coefficients de transfert
  */
  p_l  = (P_l(0) + P_l(1))/deux ;
  p_c  = p_g - p_l ;
  s_l  = courbe(p_c,el.mat->cb[0]) ;
  s_g  = un - s_l ;

  /* molarites */
  x_co2   = (X_CO2(0)   + X_CO2(1))/2. ;
  x_oh    = (X_OH(0)    + X_OH(1))/2. ;
  x_hco3  = (X_HCO3(0)  + X_HCO3(1))/2. ;

  x_h2co3 = k_h*x_co2 ;
  x_co3   = k_co3*x_oh*x_hco3 ;
  x_h     = k_e/x_oh ;
  x_ca    = k_ca/x_co3 ;
  x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;

  /* masse volumique liquide */
  rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 + M_Ca*x_ca ;
  
  /* solides */
  n_caoh2    = (N_CaOH2(0) + N_CaOH2(1))/2. ;
  n_caco3    = (N_CaCO3(0) + N_CaCO3(1))/2. ;
  n_csh      = (N_CSH(0)   + N_CSH(1))/2. ;
  n_caco3csh = 3*(n_csh0 - n_csh) ;
 
  /* porosite */
  phi = phii + v_caoh2*(n_caoh20 - n_caoh2) + v_csh*(n_csh0 - n_csh) + v_caco3*(-n_caco3) + v_caco3csh*(-n_caco3csh) ;

  /* permeabilite */
  k_l  = (k_int/mu_l)*courbe(p_c,el.mat->cb[1])*pow(phi/phii,3.)*pow(((1-phii)/(1-phi)),2.) ;
  /* tortuosite gaz */
  tau  = pow(phi,1.74)*pow(s_g,3.20) ;
  /* tortuosite liquide */
  iff    = 0.00029*exp(9.95*phi)/(1+625*pow((1-s_l),4)) ;
 
  KD_Ca     = x_ca*k_l ;
  KD_H2CO3  = x_h2co3*k_l ;
  KD_HCO3   = x_hco3*k_l ;
  KD_CO3    = x_co3*k_l ;
  KD_OH     = x_oh*k_l ;
  KD_H      = x_h*k_l ;
  KD_m      = rho_l*k_l ;

  KF_CO2    = phi*s_g*tau*d_co2 ;
  KF_Ca     = d_ca*iff ;
  KF_OH     = d_oh*iff ;
  KF_H      = d_h*iff ;
  KF_H2CO3  = d_h2co3*iff ;
  KF_HCO3   = d_hco3*iff ;
  KF_CO3    = d_co3*iff ;

  Kpsi_Ca   = FARADAY/RT*KF_Ca*z_ca*x_ca ;
  Kpsi_HCO3 = FARADAY/RT*KF_HCO3*z_hco3*x_hco3 ;
  Kpsi_CO3  = FARADAY/RT*KF_CO3*z_co3*x_co3 ;
  Kpsi_OH   = FARADAY/RT*KF_OH*z_oh*x_oh ;
  Kpsi_H    = FARADAY/RT*KF_H*z_h*x_h ;
  Kpsi_q    = z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca ;
}


void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double grd_co2,grd_p_l,grd_ca,grd_oh,grd_h,grd_h2co3,grd_hco3,grd_co3,grd_psi ;
  double w_ca,w_h,w_oh,w_hco3,w_co3,w_h2co3,w_co2,w_m,w_q ;
  double dx ;

  /* Gradients */
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

  /* Flux */
  w_ca    = - KD_Ca*grd_p_l    - KF_Ca*grd_ca     - Kpsi_Ca*grd_psi ;
  w_h     = - KD_H*grd_p_l     - KF_H*grd_h       - Kpsi_H*grd_psi ;
  w_oh    = - KD_OH*grd_p_l    - KF_OH*grd_oh     - Kpsi_OH*grd_psi ;
  w_hco3  = - KD_HCO3*grd_p_l  - KF_HCO3*grd_hco3 - Kpsi_HCO3*grd_psi ;
  w_co3   = - KD_CO3*grd_p_l   - KF_CO3*grd_co3   - Kpsi_CO3*grd_psi ;
  w_h2co3 = - KD_H2CO3*grd_p_l - KF_H2CO3*grd_h2co3 ;
  w_co2   =                    - KF_CO2*grd_co2 ;
  w_m     = - KD_m*grd_p_l + M_CO2*w_co2 ;
  w_q     =                    - z_h*KF_H*grd_h          \
                               - z_oh*KF_OH*grd_oh       \
                               - z_hco3*KF_HCO3*grd_hco3 \
                               - z_co3*KF_CO3*grd_co3    \
                               - z_ca*KF_Ca*grd_ca - Kpsi_q*grd_psi ;

  W_C     = w_co2 + w_h2co3 + w_hco3 + w_co3 ;
  W_Ca    = w_ca ;
  W_m     = w_m ;
  W_k     = w_hco3 + w_co3 ;
  W_q     = w_q ;
}

#undef NEQ

#undef NVE
#undef NVI

#undef E_C
#undef E_q
#undef E_mass
#undef E_Ca
#undef E_k
#undef E_el

#undef I_CO2
#undef I_OH
#undef I_P_l
#undef I_CaCO3
#undef I_HCO3
#undef I_psi
  
#undef X_CO2
#undef X_OH
#undef X_HCO3
#undef N_CaCO3
#undef P_l
#undef PSI

#undef N_C
#undef N_q
#undef M
#undef N_Ca
#undef N_k
#undef N_CaOH2
#undef W_C
#undef W_q
#undef W_m
#undef W_Ca
#undef W_k
#undef XI
#undef N_CSH

#undef KD_Ca
#undef KD_OH
#undef KD_H
#undef KD_H2CO3
#undef KD_HCO3
#undef KD_CO3
#undef KD_m
#undef KF_CO2
#undef KF_Ca
#undef KF_OH
#undef KF_H
#undef KF_H2CO3
#undef KF_HCO3
#undef KF_CO3
#undef Kpsi_Ca
#undef Kpsi_OH
#undef Kpsi_H
#undef Kpsi_HCO3
#undef Kpsi_CO3
#undef Kpsi_q

#undef DN_CSHSDT
#undef DN_CaOH2SDT


double concentration_oh_equilibre(void)
/* racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double a = z_co3*k_ca/k_2 ;
  double b = z_oh + z_hco3*k_ca/(k_co3*k_2) ;
  double d = z_h*k_e ;
  double e = z_ca*k_2 ;
  /* on neglige a, d */
  double x_oh  = pow(-e/b,1./3) ;
  return(x_oh) ;
}

double concentration_oh(double x_hco3)
{
  double a = z_oh + z_co3*k_co3*x_hco3 ;
  double b = z_hco3*x_hco3/a ;
  double c = (z_h*k_e + z_ca*k_ca/(k_co3*x_hco3))/a ;
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

/*
double ddn1_caoh2sdt(double av,double c_2)
{
  double dav = 1.e-3 ;
  return((dn1_caoh2sdt(av + dav,c_2) - dn1_caoh2sdt(av,c_2))/dav) ;
}
*/


#undef z_ca
#undef z_h
#undef z_oh
#undef z_hco3
#undef z_co3

#undef v_h
#undef v_oh
#undef v_h2o
#undef v_h2co3
#undef v_hco3
#undef v_co3
#undef v_ca

#undef v_caoh2
#undef v_caco3
#undef v_csh
#undef vcaco3csh

#undef M_Ca
#undef M_H2CO3
#undef M_HCO3
#undef M_CO3
#undef M_OH
#undef M_H
#undef M_H2O

#undef M_CO2

#undef M_CaOH2
#undef M_CaCO3
#undef M_CaCO3CSH
#undef M_CSH
#undef M_gel

#undef d_oh
#undef d_h
#undef d_ca
#undef d_h2co3
#undef d_hco3
#undef d_co3
#undef d_co2

#undef k_e
#undef k_h
#undef k_co3
#undef k_ca
#undef k_1
#undef k_2

#undef mu_l

#undef FARADAY
#undef RT
