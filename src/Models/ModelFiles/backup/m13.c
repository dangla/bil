#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  13
#define TITLE "Carbonatation"
#define AUTHORS "Thiery"

#include "OldMethods.h"

/* Macros */
#define NEQ     (5)
#define E_C     (0)
#define E_O     (1)
#define E_H     (2)
#define E_Ca    (3)
#define E_k     (4)
#define I_CO2   (0)
#define I_OH    (1)
#define I_P_l   (2)
#define I_CaCO3 (3)
#define I_HCO3  (4)
/* Fonctions */
static int    pm(char *s) ;
/* Parametres */
static double k_h,k_ca,k_co3,k_e,phi,d_co2,k_int,mu_l,v_ca,v_h2o,v_hco3,v_h,v_oh,v_h2co3,v_co3,a_1,k_1,a_2,k_2,n_caoh2,p_g = 0. ;

int pm(char *s)
{
  if(strcmp(s,"K_henry") == 0) return (0) ;
  else if(strcmp(s,"K_ca") == 0) return (1) ;
  else if(strcmp(s,"K_co3") == 0) return (2) ;
  else if(strcmp(s,"K_eau") == 0) return (3) ;
  else if(strcmp(s,"porosite") == 0) return (4) ;
  else if(strcmp(s,"D_co2") == 0) return (5) ;
  else if(strcmp(s,"k_int") == 0) return (6) ;
  else if(strcmp(s,"mu_l") == 0) return (7) ;
  else if(strcmp(s,"V_ca") == 0) return (8) ;
  else if(strcmp(s,"V_h2o") == 0) return (9) ;
  else if(strcmp(s,"V_hco3") == 0) return (10) ;
  else if(strcmp(s,"V_h") == 0) return (11) ;
  else if(strcmp(s,"V_oh") == 0) return (12) ;
  else if(strcmp(s,"V_h2co3") == 0) return (13) ;
  else if(strcmp(s,"V_co3") == 0) return (14) ;
  else if(strcmp(s,"A_1") == 0) return (15) ;
  else if(strcmp(s,"K_1") == 0) return (16) ;
  else if(strcmp(s,"A_2") == 0) return (17) ;
  else if(strcmp(s,"K_2") == 0) return (18) ;
  else if(strcmp(s,"N_CaOH2") == 0) return (19) ;
  else if(strcmp(s,"courbes") == 0) return (20) ;
  else { 
    printf("donnee \"%s\" non connue (pm13)\n",s) ; exit(0) ;
  }
}

int dm13(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 21 ;

  mat->neq      = NEQ ;
  strcpy(mat->eqn[E_C],"carbone") ;
  strcpy(mat->eqn[E_O],"oxygene") ;
  strcpy(mat->eqn[E_H],"hydrogene") ;
  strcpy(mat->eqn[E_Ca],"calcium") ;
  strcpy(mat->eqn[E_k],"k") ;
  strcpy(mat->inc[I_CO2],"c_co2") ;
  strcpy(mat->inc[I_OH],"c_oh") ;
  strcpy(mat->inc[I_P_l],"p_l") ;
  strcpy(mat->inc[I_CaCO3],"c_caco3") ;
  strcpy(mat->inc[I_HCO3],"c_hco3") ;
  
  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm13(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 5 equations:\n\
\t 1. l\'equation de conservation de C (c_co2)\n\
\t 2. l\'equation de conservation de O (c_oh)\n\
\t 3. l\'equation de conservation de H (p_l)\n\
\t 4. l\'equation de conservation de Ca (c_caco3)\n\
\t 5. 1 equation de cinetique (c_hco3)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"K_henry = 1.      # La constante de Henry\n") ;
  fprintf(ficd,"K_ca = 3.89e-9    # Constante d\'equilibre de CaCO3\n") ;
  fprintf(ficd,"K_co3 = 4.571e3   # Constante d\'equilibre (CO3(2-))/(HCO3-)(OH-)\n") ;
  fprintf(ficd,"K_eau = 1.e-14    # Produit ionique de l\'eau\n") ;
  fprintf(ficd,"k_int = 5e-19     # Permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 1.e-3      # Viscosite du l\'eau\n") ;
  fprintf(ficd,"D_co2 = 0.5e-3    # Coefficient de diffusion de CO2\n") ;
  fprintf(ficd,"V_oh = 18.e-3     # Volume molaire partiel de OH-\n") ;
  fprintf(ficd,"V_h2o = 18.e-3    # Volume molaire partiel de H2O\n") ;
  fprintf(ficd,"V_hco3 = 50.e-3   # Volume molaire partiel de HCO3-\n") ;
  fprintf(ficd,"V_h = 0           # Volume molaire partiel de H+\n") ;
  fprintf(ficd,"V_ca = 0          # Volume molaire partiel de Ca2+\n") ;
  fprintf(ficd,"V_h2co3 = 50.e-3  # Volume molaire partiel de H2CO3\n") ;
  fprintf(ficd,"V_co3 = 50.e-3    # Volume molaire partiel de CO3(2-\n") ;
  fprintf(ficd,"A_1 = 150         # Coef de la cinetique 1\n") ;
  fprintf(ficd,"K_1 = 2.18776e-8  # Cste d\'equilibre 1\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Coef de la cinetique \n") ;
  fprintf(ficd,"K_2 = 6.45654e-6  # Cste d\'equilibre 2\n") ;
  fprintf(ficd,"N_CaOH2 = 6.1     # Contenu initial en moles de portlandite\n") 
;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;  
  
  return(NEQ) ;
}

void tb13(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 19 ;
  el->n_ve = 6 ;
}

void ch13(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in13(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define X_CO2(n)   (u[(n)][I_CO2])
#define X_OH(n)    (u[(n)][I_OH])
#define N_CaCO3(n) (u[(n)][I_CaCO3])
#define P_l(n)     (u[(n)][I_P_l])
#define X_HCO3(n)  (u[(n)][I_HCO3])

#define N_C(n)     (f[(n)])
#define N_O(n)     (f[(2+n)])
#define N_H(n)     (f[(4+n)])
#define N_Ca(n)    (f[(6+n)])
#define N_k(n)     (f[(8+n)])
#define XI(n)      (f[(10+n)])
#define W_C        (f[12])
#define W_O        (f[13])
#define W_H        (f[14])
#define W_Ca       (f[15])
#define W_k        (f[16])
#define N_CaOH2(n) (f[(17+n)])

#define KD_C       (va[(0)])
#define KD_O       (va[(1)])
#define KD_H       (va[(2)])
#define KD_Ca      (va[(3)])
#define KD_k       (va[(4)])
#define KF_C       (va[(5)])
  double n_co2,n_h2co3,n_hco3,n_co3,n_oh,n_h,n_h2o,n_ca ;
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double grd_co2,grd_p_l,k_l,tau ;
  double s_l,s_g,p_c ;
  double dx ;
  int    i ;
  double un = 1.,deux = 2. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  k_h     = el.mat->pr[pm("K_henry")] ;
  k_ca    = el.mat->pr[pm("K_ca")] ;
  k_co3   = el.mat->pr[pm("K_co3")] ;
  k_e     = el.mat->pr[pm("K_eau")] ;
  phi     = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  d_co2   = el.mat->pr[pm("D_co2")] ;
  v_h2co3 = el.mat->pr[pm("V_h2co3")] ;
  v_hco3  = el.mat->pr[pm("V_hco3")] ;
  v_co3   = el.mat->pr[pm("V_co3")] ;
  v_h     = el.mat->pr[pm("V_h")] ;
  v_oh    = el.mat->pr[pm("V_oh")] ;
  v_h2o   = el.mat->pr[pm("V_h2o")] ;
  v_ca    = el.mat->pr[pm("V_ca")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  k_1     = el.mat->pr[pm("K_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  k_2     = el.mat->pr[pm("K_2")] ;
  n_caoh2 = el.mat->pr[pm("N_CaOH2")] ;
  
  /* Contenus molaires */
  for(i=0;i<2;i++) {
    p_c     = p_g - P_l(i) ;
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
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3\
		     + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;
    /* contenus molaires */
    n_co2   = phi*s_g*x_co2 ;
    n_oh    = phi*s_l*x_oh ;
    n_h2o   = phi*s_l*x_h2o ;
    n_hco3  = phi*s_l*x_hco3 ;
    n_h2co3 = phi*s_l*x_h2co3 ;
    n_co3   = phi*s_l*x_co3 ;
    n_h     = phi*s_l*x_h ;
    n_ca    = phi*s_l*x_ca ;
    
    N_CaOH2(i) = n_caoh2 ;
    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + N_CaCO3(i) ;
    N_O(i)  = n_h2co3 + n_hco3 + n_co3 + n_oh + n_h2o + N_CaCO3(i) \
            + 2*N_CaOH2(i) ;
    N_H(i)  = 2*n_h2co3 + n_hco3 + n_h + n_oh + 2*n_h2o + 2*N_CaOH2(i) ;
    N_Ca(i) = n_ca + N_CaOH2(i) + N_CaCO3(i) ;
    N_k(i)  = n_hco3 + n_co3 + N_CaCO3(i) ;
    XI(i)   = phi*s_l*a_1*log(x_h2co3*x_oh/(x_hco3*k_1)) ;
  }
  /* Coefficient de transfert */
  p_c  = p_g - (P_l(0)+P_l(1))/deux ;
  s_l  = courbe(p_c,el.mat->cb[0]) ;
  s_g  = un - s_l ;
  k_l  = k_int/mu_l*courbe(p_c,el.mat->cb[1]) ;
  tau  = pow(phi,1./3.)*pow(s_g,7./3.) ;

  x_co2   = (X_CO2(0)+X_CO2(1))/deux ;
  x_oh    = (X_OH(0)+X_OH(1))/deux ;
  x_hco3  = (X_HCO3(0)+X_HCO3(1))/deux ;
  x_h2co3 = k_h*x_co2 ;
  x_co3   = k_co3*x_oh*x_hco3 ;
  x_h     = k_e/x_oh ;
  x_ca    = k_ca/x_co3 ;
  x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3\
		   + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;
  KD_C  = (x_h2co3 + x_hco3 + x_co3)*k_l ;
  KD_O  = (x_h2co3 + x_hco3 + x_co3 + x_oh + x_h2o)*k_l ;
  KD_H  = (2*x_h2co3 + x_hco3 + x_h + x_oh + 2*x_h2o)*k_l ;
  KD_Ca = x_ca*k_l ;
  KD_k  = (x_hco3 + x_co3)*k_l ;
  KF_C  = phi*s_g*tau*d_co2 ;
  
  /* Flux */
  dx      = x[1][0] - x[0][0] ;
  grd_co2 = (X_CO2(1) - X_CO2(0))/dx ;
  grd_p_l = (P_l(1) - P_l(0))/dx ;
  
  W_C   = - KF_C*grd_co2 - KD_C*grd_p_l ;
  W_O   = - KD_O*grd_p_l ;
  W_H   = - KD_H*grd_p_l ;
  W_Ca  = - KD_Ca*grd_p_l ;
  W_k   = - KD_k*grd_p_l ;
  
#undef X_CO2
#undef X_OH
#undef X_HCO3
#undef N_CaCO3
#undef P_l

#undef N_C
#undef N_O
#undef N_H
#undef N_Ca
#undef N_k
#undef N_CaOH2
#undef W_C
#undef W_O
#undef W_H
#undef W_Ca
#undef W_k
#undef XI

#undef KD_C
#undef KD_O
#undef KD_H
#undef KD_Ca
#undef KD_k
#undef KF_C
}

int ex13(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define X_CO2(n)   (u[(n)][I_CO2])
#define X_OH(n)    (u[(n)][I_OH])
#define N_CaCO3(n) (u[(n)][I_CaCO3])
#define P_l(n)     (u[(n)][I_P_l])
#define X_HCO3(n)  (u[(n)][I_HCO3])

#define KD_C       (va[(0)])
#define KD_O       (va[(1)])
#define KD_H       (va[(2)])
#define KD_Ca      (va[(3)])
#define KD_k       (va[(4)])
#define KF_C       (va[(5)])
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double s_l,s_g,p_c,k_l,tau ;
  double un = 1.,deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  d_co2   = el.mat->pr[pm("D_co2")] ;
  /*
    Coefficients de transfert
  */
  p_c  = p_g - (P_l(0)+P_l(1))/deux ;
  s_l  = courbe(p_c,el.mat->cb[0]) ;
  s_g  = un - s_l ;
  k_l  = k_int/mu_l*courbe(p_c,el.mat->cb[1]) ;
  tau  = pow(phi,1./3.)*pow(s_g,7./3.) ;
  
  x_co2   = (X_CO2(0)+X_CO2(1))/deux ;
  x_oh    = (X_OH(0)+X_OH(1))/deux ;
  x_hco3  = (X_HCO3(0) + X_HCO3(1))/deux ;
  x_h2co3 = k_h*x_co2 ;
  x_co3   = k_co3*x_oh*x_hco3 ;
  x_h     = k_e/x_oh ;
  x_ca    = k_ca/x_co3 ;
  x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3\
		   + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;
  
  KD_C  = (x_h2co3 + x_hco3 + x_co3)*k_l ;
  KD_O  = (x_h2co3 + x_hco3 + x_co3 + x_oh + x_h2o)*k_l ;
  KD_H  = (2*x_h2co3 + x_hco3 + x_h + x_oh + 2*x_h2o)*k_l ;
  KD_Ca = x_ca*k_l ;
  KD_k  = (x_hco3 + x_co3)*k_l ;
  KF_C  = phi*s_g*tau*d_co2 ;
  return(0) ;
  
#undef X_CO2
#undef X_OH
#undef X_HCO3
#undef N_CaCO3
#undef P_l

#undef KD_C
#undef KD_O
#undef KD_H
#undef KD_Ca
#undef KD_k
#undef KF_C
}

int ct13(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define X_CO2(n)   (u_1[(n)][I_CO2])
#define X_OH(n)    (u_1[(n)][I_OH])
#define N_CaCO3(n) (u_1[(n)][I_CaCO3])
#define P_l(n)     (u_1[(n)][I_P_l])
#define X_HCO3(n)  (u_1[(n)][I_HCO3])

#define N_C(n)     (f_1[(n)])
#define N_O(n)     (f_1[(2+n)])
#define N_H(n)     (f_1[(4+n)])
#define N_Ca(n)    (f_1[(6+n)])
#define N_k(n)     (f_1[(8+n)])
#define XI(n)      (f_1[(10+n)])
#define W_C        (f_1[12])
#define W_O        (f_1[13])
#define W_H        (f_1[14])
#define W_Ca       (f_1[15])
#define W_k        (f_1[16])
#define N_CaOH2(n) (f_1[(17+n)])

#define N_CaOH2n(n) (f_n[(17+n)])

#define KD_C       (va[(0)])
#define KD_O       (va[(1)])
#define KD_H       (va[(2)])
#define KD_Ca      (va[(3)])
#define KD_k       (va[(4)])
#define KF_C       (va[(5)])
  double n_co2,n_h2co3,n_hco3,n_co3,n_oh,n_h,n_h2o,n_ca ;
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double grd_co2,grd_p_l ;
  double s_l,s_g,p_c ;
  double dx ;
  int    i ;
  double un = 1. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  k_h     = el.mat->pr[pm("K_henry")] ;
  k_ca    = el.mat->pr[pm("K_ca")] ;
  k_co3   = el.mat->pr[pm("K_co3")] ;
  k_e     = el.mat->pr[pm("K_eau")] ;
  phi     = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  d_co2   = el.mat->pr[pm("D_co2")] ;
  v_h2co3 = el.mat->pr[pm("V_h2co3")] ;
  v_hco3  = el.mat->pr[pm("V_hco3")] ;
  v_co3   = el.mat->pr[pm("V_co3")] ;
  v_h     = el.mat->pr[pm("V_h")] ;
  v_oh    = el.mat->pr[pm("V_oh")] ;
  v_h2o   = el.mat->pr[pm("V_h2o")] ;
  v_ca    = el.mat->pr[pm("V_ca")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  k_1     = el.mat->pr[pm("K_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  k_2     = el.mat->pr[pm("K_2")] ;

  /* Contenus molaires */
  for(i=0;i<2;i++) {
    p_c     = p_g - P_l(i) ;
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
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3\
		     + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;
    if(x_co2 <= 0. || x_oh <= 0. || x_h2o <= 0. || x_hco3 <= 0.) {
      printf("\n\
en x    = %e\n\
x_co2   = %e\n\
x_oh    = %e\n\
x_h2o   = %e\n\
x_hco3  = %e\n",x[i][0],x_co2,x_oh,x_h2o,x_hco3) ;
      return(1) ;
    }
    /* contenus molaires */
    n_co2   = phi*s_g*x_co2 ;
    n_oh    = phi*s_l*x_oh ;
    n_h2o   = phi*s_l*x_h2o ;
    n_hco3  = phi*s_l*x_hco3 ;
    n_h2co3 = phi*s_l*x_h2co3 ;
    n_co3   = phi*s_l*x_co3 ;
    n_h     = phi*s_l*x_h ;
    n_ca    = phi*s_l*x_ca ;
    
    N_CaOH2(i) = N_CaOH2n(i) + dt*a_2*log(k_ca*x_oh/(k_co3*k_2*x_hco3)) ;
    if(N_CaOH2(i) < 0.) N_CaOH2(i) = 0. ;
    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + N_CaCO3(i);
    N_O(i)  = n_h2co3 + n_hco3 + n_co3 + n_oh + n_h2o + N_CaCO3(i) \
            + 2*N_CaOH2(i) ;
    N_H(i)  = 2*n_h2co3 + n_hco3 + n_h + n_oh + 2*n_h2o + 2*N_CaOH2(i) ;
    N_Ca(i) = n_ca + N_CaOH2(i) + N_CaCO3(i) ;
    N_k(i)  = n_hco3 + n_co3 + N_CaCO3(i) ;
    XI(i)   = phi*s_l*a_1*log(x_h2co3*x_oh/(x_hco3*k_1)) ;
  }
  /* flux */
  dx      = x[1][0] - x[0][0] ;
  grd_co2 = (X_CO2(1) - X_CO2(0))/dx ;
  grd_p_l = (P_l(1) - P_l(0))/dx ;
  
  W_C   = - KF_C*grd_co2 - KD_C*grd_p_l ;
  W_O   = - KD_O*grd_p_l ;
  W_H   = - KD_H*grd_p_l ;
  W_Ca  = - KD_Ca*grd_p_l ;
  W_k   = - KD_k*grd_p_l ;

  return 0 ;

#undef X_CO2
#undef X_OH
#undef N_CaCO3
#undef P_l
#undef X_HCO3

#undef N_C
#undef N_O
#undef N_H
#undef N_Ca
#undef N_k
#undef N_CaOH2
#undef W_C
#undef W_O
#undef W_H
#undef W_Ca
#undef W_k
#undef XI

#undef N_CaOH2n

#undef KD_C
#undef KD_O
#undef KD_H
#undef KD_Ca
#undef KD_k
#undef KF_C
}

int mx13(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define X_CO2(n)   (u_1[(n)][I_CO2])
#define X_OH(n)    (u_1[(n)][I_OH])
#define N_CaCO3(n) (u_1[(n)][I_CaCO3])
#define P_l(n)     (u_1[(n)][I_P_l])
#define X_HCO3(n)  (u_1[(n)][I_HCO3])

#define N_CaOH2(n) (f_1[(17+n)])

#define KD_C       (va[(0)])
#define KD_O       (va[(1)])
#define KD_H       (va[(2)])
#define KD_Ca      (va[(3)])
#define KD_k       (va[(4)])
#define KF_C       (va[(5)])

#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double s_l,s_g,p_c,ds_lsdp_c,xi ;
  double tr_co2,tr_p_l,tr ;
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double zero = 0.,un = 1.,deux = 2. ;
  /*
    Initialisation 
  */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees 
  */
  k_h     = el.mat->pr[pm("K_henry")] ;
  k_ca    = el.mat->pr[pm("K_ca")] ;
  k_co3   = el.mat->pr[pm("K_co3")] ;
  k_e     = el.mat->pr[pm("K_eau")] ;
  phi     = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  d_co2   = el.mat->pr[pm("D_co2")] ;
  v_h2co3 = el.mat->pr[pm("V_h2co3")] ;
  v_hco3  = el.mat->pr[pm("V_hco3")] ;
  v_co3   = el.mat->pr[pm("V_co3")] ;
  v_h     = el.mat->pr[pm("V_h")] ;
  v_oh    = el.mat->pr[pm("V_oh")] ;
  v_h2o   = el.mat->pr[pm("V_h2o")] ;
  v_ca    = el.mat->pr[pm("V_ca")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  k_1     = el.mat->pr[pm("K_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  k_2     = el.mat->pr[pm("K_2")] ;
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
    p_c       = p_g - P_l(i) ;
    s_l       = courbe(p_c,el.mat->cb[0]) ;
    ds_lsdp_c = dcourbe(p_c,el.mat->cb[0]) ;
    s_g       = un - s_l ;
    /* molarites */
    x_co2   = X_CO2(i) ;
    x_oh    = X_OH(i) ;
    x_hco3  = X_HCO3(i) ;
    x_h2co3 = k_h*x_co2 ;
    xi      = a_1*log(x_h2co3*x_oh/(x_hco3*k_1)) ;
    x_co3   = k_co3*x_oh*x_hco3 ;
    x_h     = k_e/x_oh ;
    x_ca    = k_ca/x_co3 ;
    x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3\
		     + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;
    if(N_CaOH2(i) > 0.) a_2 = el.mat->pr[pm("A_2")] ; else a_2 = 0. ;
    
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_CO2+j)   += volume[i]*phi*(s_g + s_l*k_h) ;
    K(E_C+j,I_OH+j)    += volume[i]*phi*s_l*(k_co3*x_hco3) ;
    K(E_C+j,I_HCO3+j)  += volume[i]*phi*s_l*(un + k_co3*x_oh) ;
    K(E_C+j,I_CaCO3+j) += volume[i] ;
    K(E_C+j,I_P_l+j)   += volume[i]*phi*(-ds_lsdp_c)*(-x_co2 \
			+ x_h2co3 + x_hco3 + x_co3) ;
    /*
      Conservation de O (oxygene)  : (n_O1 - n_On) + dt * div(w_O) = 0
    */
    K(E_O+j,I_CO2+j)   += volume[i]*phi*s_l*(k_h - k_h*v_h2co3/v_h2o) ;
    K(E_O+j,I_OH+j)    += volume[i]*(phi*s_l*(k_co3*x_hco3 + un \
					      - (k_co3*x_hco3*v_co3 \
						 - x_h/x_oh*v_h + v_oh \
						 - x_ca/x_oh*v_ca)/v_h2o) \
				     + 2*dt*a_2/x_oh) ;
    K(E_O+j,I_HCO3+j)  += volume[i]*(phi*s_l*(un + k_co3*x_oh \
					      - (v_hco3 + k_co3*x_oh*v_co3 \
						 - x_ca/x_hco3*v_ca)/v_h2o) \
				     - 2*dt*a_2/x_hco3) ;
    K(E_O+j,I_CaCO3+j) += volume[i] ;
    K(E_O+j,I_P_l+j)   += volume[i]*phi*(-ds_lsdp_c)*(x_h2co3 + x_hco3 \
                        + x_co3 + x_oh + x_h2o) ;
    /*
      Conservation de H (hydrogene) : (n_H1 - n_Hn) + dt * div(w_H) = 0
    */
    K(E_H+j,I_CO2+j)   += volume[i]*phi*s_l*(2*k_h - 2*k_h*v_h2co3/v_h2o) ;
    K(E_H+j,I_OH+j)    += volume[i]*(phi*s_l*(-x_h/x_oh + un \
					      - 2*(k_co3*x_hco3*v_co3 \
						   - x_h/x_oh*v_h + v_oh \
						   - x_ca/x_oh*v_ca)/v_h2o) \
				     + 2*dt*a_2/x_oh) ;
    K(E_H+j,I_HCO3+j)  += volume[i]*(phi*s_l*(un \
					      - 2*(v_hco3 + k_co3*x_oh*v_co3 \
						   - x_ca/x_hco3*v_ca)/v_h2o) \
				     - 2*dt*a_2/x_hco3) ;
    K(E_H+j,I_P_l+j)   += volume[i]*phi*(-ds_lsdp_c)*(2*x_h2co3 + x_hco3 \
			+ x_h + x_oh + 2*x_h2o) ;
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_OH+j)    += volume[i]*(phi*s_l*(-x_ca/x_oh) + dt*a_2/x_oh) ;
    K(E_Ca+j,I_HCO3+j)  += volume[i]*(phi*s_l*(-x_ca/x_hco3) - dt*a_2/x_hco3) ;
    K(E_Ca+j,I_CaCO3+j) += volume[i] ;
    K(E_Ca+j,I_P_l+j)   += volume[i]*phi*(-ds_lsdp_c)*x_ca ;
    /*
      Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
    */
    K(E_k+j,I_CO2+j)   += volume[i]*phi*s_l*(-dt*a_1/x_co2) ;
    K(E_k+j,I_OH+j)    += volume[i]*phi*s_l*(k_co3*x_hco3 - dt*a_1/x_oh) ;
    K(E_k+j,I_CaCO3+j) += volume[i] ;
    K(E_k+j,I_HCO3+j)  += volume[i]*phi*s_l*(un + k_co3*x_oh + dt*a_1/x_hco3) ;
    K(E_k+j,I_P_l+j)   += volume[i]*phi*(-ds_lsdp_c)*(x_hco3 + x_co3 - dt*xi) ;
  }
  /* termes d'ecoulement */
  tr     = dt*surf/dx ;
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  tr_co2 = tr*KF_C ;
  tr_p_l = tr*KD_C ;
  K(E_C,I_CO2)         += + tr_co2 ;
  K(E_C,I_CO2+NEQ)     += - tr_co2 ;
  K(E_C,I_P_l)         += + tr_p_l ;
  K(E_C,I_P_l+NEQ)     += - tr_p_l ;

  K(E_C+NEQ,I_CO2)     += - tr_co2 ;
  K(E_C+NEQ,I_CO2+NEQ) += + tr_co2 ;
  K(E_C+NEQ,I_P_l)     += - tr_p_l ;
  K(E_C+NEQ,I_P_l+NEQ) += + tr_p_l ;
  /*
    Conservation de O (oxygene)  : (n_O1 - n_On) + dt * div(w_O) = 0
  */
  tr_p_l = tr*KD_O ;
  K(E_O,I_P_l)          += + tr_p_l ;
  K(E_O,I_P_l+NEQ)      += - tr_p_l ;

  K(E_O+NEQ,I_P_l)      += - tr_p_l ;
  K(E_O+NEQ,I_P_l+NEQ)  += + tr_p_l ;
  /*
    Conservation de H (hydrogene) : (n_H1 - n_Hn) + dt * div(w_H) = 0
  */
  tr_p_l = tr*KD_H ;
  K(E_H,I_P_l)          += + tr_p_l ;
  K(E_H,I_P_l+NEQ)      += - tr_p_l ;

  K(E_H+NEQ,I_P_l)      += - tr_p_l ;
  K(E_H+NEQ,I_P_l+NEQ)  += + tr_p_l ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  tr_p_l = tr*KD_Ca ;
  K(E_Ca,I_P_l)         += + tr_p_l ;
  K(E_Ca,I_P_l+NEQ)     += - tr_p_l ;

  K(E_Ca+NEQ,I_P_l)     += - tr_p_l ;
  K(E_Ca+NEQ,I_P_l+NEQ) += + tr_p_l ;
  /*
    Cinetique 1 : (n_k1 - n_kn) + dt * div(W_k) - dt * XI = 0
  */
  tr_p_l = tr*KD_k ;
  K(E_k,I_P_l)          += + tr_p_l ;
  K(E_k,I_P_l+NEQ)      += - tr_p_l ;

  K(E_k+NEQ,I_P_l)      += - tr_p_l ;
  K(E_k+NEQ,I_P_l+NEQ)  += + tr_p_l ;

  /*
  printf("\nmatrice\n") ;
  for(i=0;i<NEQ;i++) {
    printf("\n") ;
    for(j=0;j<NEQ;j++) printf("%e ",K(i,j)) ;
  }
  */

  return(0) ;

#undef X_CO2
#undef X_OH
#undef X_HCO3
#undef N_CaCO3
#undef P_l

#undef N_CaOH2

#undef KD_C
#undef KD_O
#undef KD_H
#undef KD_Ca
#undef KD_k
#undef KF_C

#undef K
}

void rs13(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define X_CO2(n)   (u_1[(n)][I_CO2])
#define X_OH(n)    (u_1[(n)][I_OH])
#define N_CaCO3(n) (u_1[(n)][I_CaCO3])
#define P_l(n)     (u_1[(n)][I_P_l])
#define X_HCO3(n)  (u_1[(n)][I_HCO3])

#define N_C1(n)    (f_1[(n)])
#define N_O1(n)    (f_1[(2+n)])
#define N_H1(n)    (f_1[(4+n)])
#define N_Ca1(n)   (f_1[(6+n)])
#define N_k1(n)    (f_1[(8+n)])
#define XI(n)      (f_1[(10+n)])
#define W_C        (f_1[12])
#define W_O        (f_1[13])
#define W_H        (f_1[14])
#define W_Ca       (f_1[15])
#define W_k        (f_1[16])

#define N_Cn(n)    (f_n[(n)])
#define N_On(n)    (f_n[(2+n)])
#define N_Hn(n)    (f_n[(4+n)])
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
  R(0,E_C) -= volume[0]*(N_C1(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C1(1) - N_Cn(1)) - dt*surf*W_C ;
  /*
    Conservation de O (oxygene)  : (n_O1 - n_On) + dt * div(w_O) = 0
  */
  R(0,E_O) -= volume[0]*(N_O1(0) - N_On(0)) + dt*surf*W_O ;
  R(1,E_O) -= volume[1]*(N_O1(1) - N_On(1)) - dt*surf*W_O ;
  /*
    Conservation de H (hydrogene) : (n_H1 - n_Hn) + dt * div(w_H) = 0
  */
  R(0,E_H) -= volume[0]*(N_H1(0) - N_Hn(0)) + dt*surf*W_H ;
  R(1,E_H) -= volume[1]*(N_H1(1) - N_Hn(1)) - dt*surf*W_H ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca1(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca1(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
  */
  R(0,E_k) -= volume[0]*(N_k1(0) - N_kn(0) - dt*XI(0)) + dt*surf*W_k ;
  R(1,E_k) -= volume[1]*(N_k1(1) - N_kn(1) - dt*XI(1)) - dt*surf*W_k ;
  
#undef X_CO2
#undef X_OH
#undef N_CaCO3
#undef P_l
#undef X_HCO3

#undef N_C1
#undef N_O1
#undef N_H1
#undef N_Ca1
#undef N_k1
#undef W_C
#undef W_O
#undef W_H
#undef W_Ca
#undef W_k
#undef XI

#undef N_Cn
#undef N_On
#undef N_Hn
#undef N_Can
#undef N_kn

#undef R
}

int so13(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define X_CO2(n)   (u[(n)][I_CO2])
#define X_OH(n)    (u[(n)][I_OH])
#define N_CaCO3(n) (u[(n)][I_CaCO3])
#define P_l(n)     (u[(n)][I_P_l])
#define X_HCO3(n)  (u[(n)][I_HCO3])

#define N_C(n)     (f[(n)])
#define N_O(n)     (f[(2+n)])
#define N_H(n)     (f[(4+n)])
#define N_Ca(n)    (f[(6+n)])
#define N_k(n)     (f[(8+n)])
#define XI(n)      (f[(10+n)])
#define W_C        (f[12])
#define W_O        (f[13])
#define W_H        (f[14])
#define W_Ca       (f[15])
#define W_k        (f[16])
#define N_CaOH2(n) (f[(17+n)])
  int    k ;
  double x_co2,x_h2co3,x_hco3,x_co3,x_oh,x_h,x_h2o,x_ca ;
  double n_c,n_o,n_h,n_ca,n_caco3,n_caoh2 ;
  double xi ;
  double s_l,p_c,p_l ;
  double dx ;
  double un = 1.,deux = 2. ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  k_h     = el.mat->pr[pm("K_henry")] ;
  k_ca    = el.mat->pr[pm("K_ca")] ;
  k_co3   = el.mat->pr[pm("K_co3")] ;
  k_e     = el.mat->pr[pm("K_eau")] ;
  phi     = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  d_co2   = el.mat->pr[pm("D_co2")] ;
  v_h2co3 = el.mat->pr[pm("V_h2co3")] ;
  v_hco3  = el.mat->pr[pm("V_hco3")] ;
  v_co3   = el.mat->pr[pm("V_co3")] ;
  v_h     = el.mat->pr[pm("V_h")] ;
  v_oh    = el.mat->pr[pm("V_oh")] ;
  v_h2o   = el.mat->pr[pm("V_h2o")] ;
  v_ca    = el.mat->pr[pm("V_ca")] ;
  a_1     = el.mat->pr[pm("A_1")] ;
  k_1     = el.mat->pr[pm("K_1")] ;
  a_2     = el.mat->pr[pm("A_2")] ;
  k_2     = el.mat->pr[pm("K_2")] ;

  dx      = x[1][0] - x[0][0] ;
  p_l     = (P_l(0) + P_l(1))/deux ;
  p_c     = p_g - p_l ;
  s_l     = courbe(p_c,el.mat->cb[0]) ;
  /* molarites */
  x_co2   = (X_CO2(0) + X_CO2(1))/deux ;
  x_oh    = (X_OH(0) + X_OH(1))/deux ;
  x_hco3  = (X_HCO3(0) + X_HCO3(1))/deux ;
  xi      = (XI(0) + XI(1))/deux ;
  x_h2co3 = k_h*x_co2 ;
  x_co3   = k_co3*x_oh*x_hco3 ;
  x_h     = k_e/x_oh ;
  x_ca    = k_ca/x_co3 ;
  x_h2o   = (un - (x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3\
		   + x_h*v_h + x_oh*v_oh + x_ca*v_ca))/v_h2o ;
  /* contenus molaires */
  n_c     = (N_C(0) + N_C(1))/deux ;
  n_o     = (N_O(0) + N_O(1))/deux ;
  n_h     = (N_H(0) + N_H(1))/deux ;
  n_ca    = (N_Ca(0) + N_Ca(1))/deux ;
  n_caco3 = (N_CaCO3(0) + N_CaCO3(1))/deux ;
  n_caoh2 = (N_CaOH2(0) + N_CaOH2(1))/deux ;
  /* quantites exploitees */
  k = 0 ;
  strcpy(r[k].text,"x_co2") ; r[k].n = 1 ;
  r[k++].v[0] = x_co2 ;
  strcpy(r[k].text,"x_oh") ; r[k].n = 1 ;
  r[k++].v[0] = x_oh ;
  strcpy(r[k].text,"x_h2o") ; r[k].n = 1 ;
  r[k++].v[0]= x_h2o ;
  strcpy(r[k].text,"n_caco3") ; r[k].n = 1 ;
  r[k++].v[0] = n_caco3 ;
  strcpy(r[k].text,"s_l") ; r[k].n = 1 ;
  r[k++].v[0]= s_l ;
  strcpy(r[k].text,"xi") ; r[k].n = 1 ;
  r[k++].v[0]= xi ;
  strcpy(r[k].text,"x_hco3") ; r[k].n = 1 ;
  r[k++].v[0] = x_hco3 ;
  strcpy(r[k].text,"x_caoh2") ; r[k].n = 1 ;
  r[k++].v[0] = n_caoh2 ;
  strcpy(r[k].text,"x_ca") ; r[k].n = 1 ;
  r[k++].v[0] = x_ca ;
  strcpy(r[k].text,"n_c") ; r[k].n = 1 ;
  r[k++].v[0]= n_c ;
  strcpy(r[k].text,"n_o+2n_c") ; r[k].n = 1 ;
  r[k++].v[0]= n_o + 2*n_c ;
  strcpy(r[k].text,"n_h") ; r[k].n = 1 ;
  r[k++].v[0] = n_h ;
  strcpy(r[k].text,"n_ca") ; r[k].n = 1 ;
  r[k++].v[0] = n_ca ;
  strcpy(r[k].text,"masse") ; r[k].n = 1 ;
  r[k++].v[0] = n_c*12 + (n_o+2*n_c)*16 + n_h + n_ca*40 ;
  return(k) ;

#undef X_CO2
#undef X_OH
#undef N_CaCO3
#undef P_l
#undef X_HCO3

#undef N_C
#undef N_O
#undef N_H
#undef N_Ca
#undef N_k
#undef W_C
#undef W_O
#undef W_H
#undef W_Ca
#undef W_k
#undef XI
#undef N_CaOH2
}
