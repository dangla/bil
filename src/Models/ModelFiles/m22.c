#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  22
#define TITLE "Transport ionique dans le beton (1D)"
#define AUTHORS "Dridi"

#include "OldMethods.h"

/* Macros */
#define NEQ       (7)

#define NVI       (18)
#define NVE       (21)

#define E_mass    (0)
#define E_fe      (1)
#define E_q       (2)
#define E_o2      (3)
#define E_cat     (4)
#define E_ani     (5)
#define E_feoh2_l (6)

#define I_p_l     (0)
#define I_c_fe    (6)
#define I_psi     (2)
#define I_c_o2    (3)
#define I_c_cat   (4)
#define I_c_ani   (5)
#define I_c_feoh2 (1)


#define P_l(n)      (u[(n)][I_p_l])
#define C_Fe(n)     (u[(n)][I_c_fe])
#define PSI(n)      (u[(n)][I_psi])
#define C_O2(n)     (u[(n)][I_c_o2])
#define C_cat(n)    (u[(n)][I_c_cat])
#define C_ani(n)    (u[(n)][I_c_ani])
#define C_FeOH2(n)  (u[(n)][I_c_feoh2])

#define M(n)        (f[(n)])
#define N_Fe(n)     (f[(4+n)])
#define N_O2(n)     (f[(6+n)])
#define N_cat(n)    (f[(8+n)])
#define N_ani(n)    (f[(10+n)])
#define W_m         (f[12])
#define W_Fe        (f[13])
#define W_q         (f[14])
#define W_O2        (f[15])
#define W_cat       (f[16])
#define W_ani       (f[17])

#define M_n(n)      (f_n[(n)])
#define N_Fen(n)    (f_n[(4+n)])
#define N_O2n(n)    (f_n[(6+n)])
#define N_catn(n)   (f_n[(8+n)])
#define N_anin(n)   (f_n[(10+n)])

#define KD_m        (va[(0)])
#define KD_H        (va[(1)])
#define KD_OH       (va[(2)])
#define KD_Fe       (va[(3)])
#define KD_FeOH2_l  (va[(4)])
#define KD_O2       (va[(5)])
#define KD_cat      (va[(6)])
#define KD_ani      (va[(7)])

#define KF_OH       (va[(8)])
#define KF_H        (va[(9)])
#define KF_Fe       (va[(10)])
#define KF_FeOH2_l  (va[(11)])
#define KF_O2_l     (va[(12)])
#define KF_O2_g     (va[(13)])
#define KF_cat      (va[(14)])
#define KF_ani      (va[(15)])

#define KE_H        (va[(16)])
#define KE_OH       (va[(17)])
#define KE_Fe       (va[(18)])
#define KE_cat      (va[(19)])
#define KE_ani      (va[(20)])

/* valences */
#define z_fe       (2.)
#define z_oh       (-1.)
#define z_h        (1.)
#define z_cat      (1.)
#define z_ani      (-1.)

/* Masses molaires (unite arbitraire = M_h) */
#define M_fe       (55.85)
#define M_oh       (17.)
#define M_h        (1.)
#define M_h2o      (18.)
#define M_o2       (32.)
#define M_cat      (40.)
#define M_ani      (40.)
#define M_feoh2    (89.85)

/* coefficients de diffusion moleculaire (m2/s) */
#define d_oh       (5.273e-9)
#define d_h        (9.310e-9)
#define d_o2_l     (2.51e-9)
#define d_o2_g     (1.e-7)
#define d_fe       (0.72e-9)
#define d_feoh2_l  (0.72e-9)
#define d_cat      (1.33e-9)
#define d_ani      (1.18e-9)

/* constante d equilibre */
#define k_eau      (1.e-8)
#define k_hen      (0.03045)    /* (1.25e-5)*(RT) */
#define k_feoh2_l  (3.1429e-2)
#define s_feoh2    (1.995e-5)

/* constantes physiques */
#define FARADAY    (9.64846e4) /* cste de Faraday (C/mole) */
#define RT         (2436.)     /* produit de R et T */
                               /* R = 8.3143 J/K/mol cste des gaz parfait */
                               /* T = 293 K */

/* viscosite */
#define mu_l       (1.e-3)
/* proprietes physiques de l'eau */
#define c_h2o      (55.5e3)    /* molarite */

/* Fonctions */
static int    pm(const char *s) ;
static double saturation(double,crbe_t) ;
static double dsaturation(double,crbe_t) ;
static double tortuosite_l(double,mate_t*) ;
static double tortuosite_g(double,mate_t*) ;
static double concentration_oh(double,double) ;
static double concentration_feoh2_l(double,double) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
/* Parametres */
static double phi,k_int,p_g=0 ;
/*
static double d_oh,d_h,d_o2_l,d_o2_g,d_fe,d_feoh2_l,d_cat,d_ani ;
static double k_eau,k_feoh2_l,k_hen,s_feoh2 ;
static double c_h2o,mu_l ;
static double FARADAY,RT ;
static const double M_h = 1,M_o = 16,M_fe = 55.85,M_cat = 40,M_ani = 40 ;
static const double M_h2o = 18,M_o2 = 32,M_oh = 17,M_feoh2 = 89.85 ;
*/

/* Macros pour les fonctions saturation */
#define SATURATION(a,b)   courbe(a,b)
#define DSATURATION(a,b)  dcourbe(a,b)
#define SATURATION(a,b)   saturation(a,b)
#define DSATURATION(a,b)  dsaturation(a,b)


int pm(const char *s)
{
if(strcmp(s,"phi") == 0) return (0) ;
else if(strcmp(s,"k_int") == 0) return (1) ;
else if(strcmp(s,"c_h2o") == 0) return (2) ;
else if(strcmp(s,"mu_l") == 0) return (3) ;
else if(strcmp(s,"RT_0") == 0) return (4) ;
else if(strcmp(s,"K_eau") == 0) return (5) ;
else if(strcmp(s,"K_Fe(OH)2") == 0) return (6) ;
else if(strcmp(s,"K_Far") == 0) return (7) ;
else if(strcmp(s,"D_oh") == 0) return (8) ;
else if(strcmp(s,"D_h") == 0) return (9) ;
else if(strcmp(s,"D_o2_l") == 0) return (10) ;
else if(strcmp(s,"D_o2_g") == 0) return (11) ;
else if(strcmp(s,"D_fe") == 0) return (12) ;
else if(strcmp(s,"D_feoh2") == 0) return (13) ;
else if(strcmp(s,"D_cat") == 0) return (14) ;
else if(strcmp(s,"D_ani") == 0) return (15) ;
else if(strcmp(s,"K_Hen") == 0) return (16) ;
else if(strcmp(s,"S_Fe(OH)2") == 0) return (17) ;
else
{ printf("donnee \"%s\" non connue (pm22)\n",s) ; exit(0) ; }
}

int dm22(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 19 ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_mass],   "E_mass") ;
  strcpy(mat->eqn[E_fe],     "E_fe") ;
  strcpy(mat->eqn[E_q],      "E_q") ;
  strcpy(mat->eqn[E_o2],     "E_o2") ;
  strcpy(mat->eqn[E_cat],    "E_cat") ;
  strcpy(mat->eqn[E_ani],    "E_ani") ;
  strcpy(mat->eqn[E_feoh2_l],"E_feoh2_l") ;

  strcpy(mat->inc[I_p_l],    "p_l") ;
  strcpy(mat->inc[I_c_fe],   "c_fe") ;
  strcpy(mat->inc[I_psi],    "psi") ;
  strcpy(mat->inc[I_c_o2],   "c_o2") ;
  strcpy(mat->inc[I_c_cat],  "c_cat") ;
  strcpy(mat->inc[I_c_ani],  "c_ani") ;
  strcpy(mat->inc[I_c_feoh2],"c_feoh2") ;

  dmat(mat,ficd,pm,n_donnees) ;
  
  return(mat->n) ;
}

int qm22(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de :\n\
\t 1. Conservation de O2                         (c_o2)    \n\
\t 2. Conservation de la masse liquide et O2 gaz (p_l)     \n\
\t 3. Conservation de Fe                         (c_feoh2) \n\
\t 4. Conservation de la charge                  (psi)     \n\
\t 5. Conservation de K+                         (c_cat)   \n\
\t 6. Conservation de A-                         (c_ani)   \n\
\t 7. Equilibre de Fe(OH)2(l)                    (c_fe)    \n\
") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"phi = 0.28            # La porosite\n") ;
  fprintf(ficd,"k_int = 3.e-21        # Permeabilite intrinseque\n") ;
  fprintf(ficd,"c_h2o = 5.55e4        # Concentration molaire de l\'eau\n") ;
  fprintf(ficd,"mu_l = 1.e-3          # Viscosite du liquide\n") ;
  fprintf(ficd,"RT = 2479             # Constante des gaz parfaits fois la temperature\n") ;
  fprintf(ficd,"K_eau = 1.e-8         # Constante d\'equilibre pour l\'eau\n") ;
  fprintf(ficd,"K_Fe(OH)2 = 3.1429e-2 # Constante d\'equilibre pour Fe(OH)2\n") ;
  fprintf(ficd,"S_Fe(OH)2 = 1.995e-5  # Solubilite de Fe(OH)2 solide\n") ;
  fprintf(ficd,"K_Far = 96485         # Constante de Faraday\n") ;
  fprintf(ficd,"K_Hen = 1.25e-5       # Constante de Henry pour O2\n") ;
  fprintf(ficd,"D_oh = 5.24e-9        # Diffusion moleculaire de OH-\n") ;
  fprintf(ficd,"D_h = 9.32e-9         # Diffusion moleculaire de H+\n") ;
  fprintf(ficd,"D_o2_l = 2.51e-9      # Diffusion moleculaire de O2 liquide\n") ;
  fprintf(ficd,"D_o2_g = 1.e-7        # Diffusion moleculaire de O2 gazeux\n") ;
  fprintf(ficd,"D_fe = 0.72e-9        # Diffusion moleculaire de Fe2+\n") ;
  fprintf(ficd,"D_feoh2 = 0.72e-9     # Diffusion moleculaire de Fe(OH)2\n") ;
  fprintf(ficd,"D_cat = 1.33e-9       # Diffusion moleculaire des cations\n") ;
  fprintf(ficd,"D_ani = 1.18e-9       # Diffusion moleculaire des anions\n") ;
  fprintf(ficd,"courbes = my_file     # Nom du fichier : p_c S_l k_rl tau_l/tau_l_sat tau_g/tau_g_sat\n") ;
  
  return(NEQ) ;
}

void tb22(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}

void ch22(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in22(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double p_l,p_c,s_l,s_g ;
  double c_h,c_oh,c_fe,c_feoh2_l,c_o2_l,c_cat,c_ani,c_o2_g,c_feoh2 ;
  double rho_l ;
  int    i ;
  double un = 1. ;

  if(el.dim < dim) return ;

  /* Donnees */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  /*
  c_h2o   = el.mat->pr[pm("c_h2o")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  k_eau   = el.mat->pr[pm("K_eau")] ;
  k_feoh2_l = el.mat->pr[pm("K_Fe(OH)2")] ;
  FARADAY   = el.mat->pr[pm("K_Far")] ;
  d_oh   = el.mat->pr[pm("D_oh")] ;
  d_h    = el.mat->pr[pm("D_h")] ;
  d_o2_l = el.mat->pr[pm("D_o2_l")] ;
  d_o2_g = el.mat->pr[pm("D_o2_g")] ;
  d_fe   = el.mat->pr[pm("D_fe")] ;
  d_feoh2_l = el.mat->pr[pm("D_feoh2")] ;
  d_cat  = el.mat->pr[pm("D_cat")] ;
  d_ani  = el.mat->pr[pm("D_ani")] ;
  RT    = el.mat->pr[pm("RT")] ;
  k_hen   = el.mat->pr[pm("K_Hen")]*RT ;
  s_feoh2 = el.mat->pr[pm("S_Fe(OH)2")] ;
  */

  
  /* Contenus massiques */
  for(i=0;i<2;i++) {
    p_l       = P_l(i) ;
    p_c       = p_g - p_l ;
    s_l       = SATURATION(p_c,el.mat->cb[0]) ;
    s_g       = un - s_l ; 
    
    /* concentrations */
    c_fe      = C_Fe(i) ;
    c_o2_l    = C_O2(i) ;
    c_cat     = C_cat(i) ;
    c_ani     = C_ani(i) ;
    c_feoh2   = C_FeOH2(i) ;
    c_o2_g    = c_o2_l/k_hen ;
    c_oh      = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
    c_h       = k_eau/c_oh ;
    c_feoh2_l = concentration_feoh2_l(c_feoh2,s_feoh2) ;

    /* masse volumique de la solution + le solide feoh2 */
    rho_l = M_h2o*c_h2o + M_h*c_h + M_oh*c_oh + M_o2*c_o2_l + M_fe*c_fe + M_feoh2*c_feoh2 + M_cat*c_cat + M_ani*c_ani ;

    /* contenus molaires */
    N_O2(i)   = phi*(s_l*c_o2_l + s_g*c_o2_g) ;
    N_Fe(i)   = phi*s_l*(c_fe + c_feoh2) ;
    N_cat(i)  = phi*s_l*c_cat ;
    N_ani(i)  = phi*s_l*c_ani ;

    /* masse liquide et O2 gaz */
    M(i)      = phi*s_l*rho_l + phi*s_g*M_o2*c_o2_g ;
  }

  {
    ex_t ex22 ;
    ex22(x,u,f,va,el,dim,geom,0.) ;
  }

  flux(x,u,f,va,el,dim,geom) ;
} 

int ex22(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double p_l,p_c,s_l,k_l,tau_l,tau_g,s_g ;
  double c_h,c_oh,c_fe,c_feoh2_l,c_o2_l,c_cat,c_ani,c_feoh2 ;
  double rho_l ;
  double un = 1. ;
  int    i ;

  if(el.dim < dim) return(0) ;

  /* Donnees */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;

  for(i=0;i<NVE;i++) va[i] = 0. ;

  for(i=0;i<2;i++) {
    p_l       = P_l(i) ;
    c_fe      = C_Fe(i) ;
    c_o2_l    = C_O2(i) ;
    c_cat     = C_cat(i) ;
    c_ani     = C_ani(i) ;
    c_feoh2   = C_FeOH2(i) ;
  
    c_oh      = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
    c_h       = k_eau/c_oh ;
    c_feoh2_l = concentration_feoh2_l(c_feoh2,s_feoh2) ;

    /* masse volumique de la solution */
    rho_l = M_h2o*c_h2o + M_h*c_h + M_oh*c_oh + M_o2*c_o2_l + M_fe*c_fe + M_feoh2*c_feoh2_l + M_cat*c_cat + M_ani*c_ani ;
    
    if(c_oh < 0. || c_fe < 0. || c_cat < 0. || c_ani < 0. || c_o2_l < 0. || c_feoh2 < 0.) {
      double y = x[i][0] ;
      printf("\n\
    En x    = %e\n\
    c_oh    = %e\n\
    c_fe    = %e\n\
    c_cat   = %e\n\
    c_ani   = %e\n\
    c_o2    = %e\n\
    c_feoh2 = %e\n"\
	     ,y,c_oh,c_fe,c_cat,c_ani,c_o2_l,c_feoh2) ;
      return(1) ;
    }

    p_c    = p_g - p_l ;
    s_l    = SATURATION(p_c,el.mat->cb[0]) ;
    s_g    = un - s_l ;
    k_l    = k_int/mu_l*courbe(p_c,el.mat->cb[1]) ;
    
    tau_l  = tortuosite_l(p_c,el.mat) ;
    tau_g  = tortuosite_g(p_c,el.mat) ;
    
    
    /* coefficients de transfert */
    KD_H       += c_h*k_l ;
    KD_OH      += c_oh*k_l ;
    KD_Fe      += c_fe*k_l ;
    KD_FeOH2_l += c_feoh2_l*k_l ;
    KD_O2      += c_o2_l*k_l ;
    KD_cat     += c_cat*k_l ;
    KD_ani     += c_ani*k_l ;
    KD_m       += rho_l*k_l ;
    
    KF_OH      += phi*s_l*tau_l*d_oh ;
    KF_H       += phi*s_l*tau_l*d_h ;
    KF_Fe      += phi*s_l*tau_l*d_fe ;
    KF_FeOH2_l += phi*s_l*tau_l*d_feoh2_l ;
    KF_O2_l    += phi*s_l*tau_l*d_o2_l ;
    KF_O2_g    += phi*s_g*tau_g*d_o2_g ;
    KF_cat     += phi*s_l*tau_l*d_cat ;
    KF_ani     += phi*s_l*tau_l*d_ani ;
    
    KE_H       += z_h*c_h*FARADAY/RT*KF_H ;
    KE_OH      += z_oh*c_oh*FARADAY/RT*KF_OH ;
    KE_Fe      += z_fe*c_fe*FARADAY/RT*KF_Fe ;
    KE_cat     += z_cat*c_cat*FARADAY/RT*KF_cat ;
    KE_ani     += z_ani*c_ani*FARADAY/RT*KF_ani ;
  }

  /* moyenne */
  for(i=0;i<NVE;i++) va[i] *= 0.5 ;

  return(0) ;
}

int ct22(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double p_c,s_l,s_g,p_l ;
  double c_h,c_oh,c_fe,c_feoh2_l,c_o2_l,c_cat,c_ani,c_o2_g,c_feoh2 ;
  double rho_l ;
  int    i ;
  
  if(el.dim < dim) return(0) ;

  /* Donnees */
  phi     = el.mat->pr[pm("phi")] ;

  /* Contenus massiques */
  for(i=0;i<2;i++) {
    /* pression */
    p_l       = P_l(i) ;
    p_c       = p_g - p_l ;
    /* saturation */
    s_l       = SATURATION(p_c,el.mat->cb[0]) ;
    s_g       = 1. - s_l ;
    /* concentrations */

    /* Essai */
    /*
    if(C_Fe(i) < 0.) C_Fe(i) = 0. ;
    if(C_FeOH2(i) < 0.) C_FeOH2(i) = 0. ;
    */
    /* Fin de l'essai */

    c_fe      = C_Fe(i) ;
    c_o2_l    = C_O2(i) ;
    c_cat     = C_cat(i) ;
    c_ani     = C_ani(i) ;
    c_feoh2   = C_FeOH2(i) ;

    c_o2_g    = c_o2_l/k_hen ;
    c_oh      = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
    c_h       = k_eau/c_oh ;
    c_feoh2_l = concentration_feoh2_l(c_feoh2,s_feoh2) ;

    /* masse volumique liquide et le solide feoh2_s */
    rho_l = M_h2o*c_h2o + M_h*c_h + M_oh*c_oh + M_o2*c_o2_l + M_fe*c_fe + M_feoh2*c_feoh2 + M_cat*c_cat + M_ani*c_ani ;

    /* contenus molaires */
    N_O2(i)   = phi*(s_l*c_o2_l + s_g*c_o2_g) ;
    N_Fe(i)   = phi*s_l*(c_fe + c_feoh2) ;
    N_cat(i)  = phi*s_l*c_cat ; 
    N_ani(i)  = phi*s_l*c_ani ;

    /* masse liquide et O2 gaz */
    M(i)      = phi*s_l*rho_l + phi*s_g*M_o2*c_o2_g ;

    if(c_oh < 0. || c_fe < 0. || c_cat < 0. || c_ani < 0. || c_o2_l < 0. || c_feoh2 < 0.) {
      printf("\n\
      En x    = %e\n\
      c_oh    = %e\n\
      c_fe    = %e\n\
      c_cat   = %e\n\
      c_ani   = %e\n\
      c_o2    = %e\n\
      c_feoh2 = %e\n"\
	     ,x[i][0],c_oh,c_fe,c_cat,c_ani,c_o2_l,c_feoh2) ;
      return(1) ;
    }
  }

  flux(x,u,f,va,el,dim,geom) ;
  return(0) ;
} 


int mx22(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double p_c,s_l,ds_lsdp_c,s_g;
  double c_h,c_oh,c_fe,c_feoh2_l,c_o2_l,c_cat,c_ani,c_o2_g,c_feoh2 ;
  double rho_l ;
  double dc_hsdc_fe[2],dc_hsdc_cat[2],dc_hsdc_ani[2] ;
  double dc_ohsdc_fe[2],dc_ohsdc_cat[2],dc_ohsdc_ani[2] ;
  double dc_feoh2_lsdc_feoh2[2] ;
  double dc_o2_gsdc_o2_l ;
  double drho_lsdc_fe,drho_lsdc_o2_l,drho_lsdc_feoh2,drho_lsdc_cat,drho_lsdc_ani ;
  double tr ;
  double trd_h,trd_oh,trd_fe,trd_feoh2_l,trd_o2,trd_cat,trd_ani,trd_m ;
  double trf_oh,trf_h,trf_fe,trf_feoh2_l,trf_o2,trf_o2_l,trf_o2_g,trf_cat,trf_ani ;
  double tre_h,tre_oh,tre_fe,tre_cat,tre_ani ;
  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  double c[2] ;

  /* initialisation */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = zero ;
  
  if(el.dim < dim) return(0) ;

  /* Donnees */
  phi     = el.mat->pr[pm("phi")] ;

  /* CALCUL DE volume ET DE surf */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++) {
  volume[i] = fabs(dx)*0.5 ;
  if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ;
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;

 /* termes d'accumulation */
  for(i=0;i<2;i++) {
    p_c       = p_g - P_l(i) ;
    s_l       = SATURATION(p_c,el.mat->cb[0]) ;
    s_g       = un - s_l ;
    ds_lsdp_c = DSATURATION(p_c,el.mat->cb[0]) ;
    /* concentrations */
    c_fe      = C_Fe(i) ;
    c_o2_l    = C_O2(i) ;
    c_cat     = C_cat(i) ;
    c_ani     = C_ani(i) ;
    c_feoh2   = C_FeOH2(i) ;

    c_o2_g    = c_o2_l/k_hen ;
    c_oh      = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
    c_h       = k_eau/c_oh ;
    c_feoh2_l = concentration_feoh2_l(c_feoh2,s_feoh2) ;

    /* masse volumique liquide et le solide feoh2_s */
    rho_l = M_h2o*c_h2o + M_h*c_h + M_oh*c_oh + M_o2*c_o2_l + M_fe*c_fe + M_feoh2*c_feoh2 + M_cat*c_cat + M_ani*c_ani ;


    /* derivees */
    dc_ohsdc_fe[i]     = z_fe*c_oh/(c_h + c_oh) ;
    dc_ohsdc_cat[i]    = z_cat*c_oh/(c_h + c_oh) ;
    dc_ohsdc_ani[i]    = z_ani*c_oh/(c_h + c_oh) ;
    
    dc_hsdc_fe[i]      = -k_eau/(c_oh*c_oh)*dc_ohsdc_fe[i] ;
    dc_hsdc_cat[i]     = -k_eau/(c_oh*c_oh)*dc_ohsdc_cat[i] ;
    dc_hsdc_ani[i]     = -k_eau/(c_oh*c_oh)*dc_ohsdc_ani[i] ;

    dc_o2_gsdc_o2_l = 1./k_hen ;
    
    /*
      dc_feoh2_lsdc_fe[i]  = (c_oh*c_oh + 2*c_fe*c_oh*dc_ohsdc_fe[i])/k_feoh2_l ;
      dc_feoh2_lsdc_cat[i] = 2*c_fe*c_oh*dc_ohsdc_cat[i]/k_feoh2_l ;
      dc_feoh2_lsdc_ani[i] = 2*c_fe*c_oh*dc_ohsdc_ani[i]/k_feoh2_l ;
    */

    dc_feoh2_lsdc_feoh2[i] = (c_feoh2 < s_feoh2) ? 1 : 0 ;

    drho_lsdc_fe    = M_h*dc_hsdc_fe[i] + M_oh*dc_ohsdc_fe[i] + M_fe ;
    drho_lsdc_o2_l  = M_o2 ;
    drho_lsdc_feoh2 = M_feoh2 ;
    drho_lsdc_cat   = M_h*dc_hsdc_cat[i] + M_oh*dc_ohsdc_cat[i] + M_cat ;
    drho_lsdc_ani   = M_h*dc_hsdc_ani[i] + M_oh*dc_ohsdc_ani[i] + M_ani ;

    /*
    CONSERVATION DE LA MASSE LIQUIDE + O2 GAZ : (m_1 - m_n) + dt * div(w_m1) = 0
    */
    K(E_mass+i*NEQ,I_p_l+i*NEQ)       += volume[i]*phi*(-ds_lsdp_c)*(rho_l - M_o2*c_o2_g) ;
    K(E_mass+i*NEQ,I_c_o2+i*NEQ)      += volume[i]*phi*(s_l*drho_lsdc_o2_l + s_g*M_o2*dc_o2_gsdc_o2_l) ;
    K(E_mass+i*NEQ,I_c_fe+i*NEQ)      += volume[i]*phi*s_l*(drho_lsdc_fe) ;
    K(E_mass+i*NEQ,I_c_cat+i*NEQ)     += volume[i]*phi*s_l*(drho_lsdc_cat) ;
    K(E_mass+i*NEQ,I_c_ani+i*NEQ)     += volume[i]*phi*s_l*(drho_lsdc_ani) ;
    K(E_mass+i*NEQ,I_c_feoh2+i*NEQ)   += volume[i]*phi*s_l*(drho_lsdc_feoh2) ; 
    /*
    CONSERVATION DE L'ELEMENT Fe = Fe2+ + Fe(OH)2  :  (m_fe1 - m_fen) + dt * div(w_fe1) = 0
    */
    K(E_fe+i*NEQ,I_p_l+i*NEQ)       += volume[i]*phi*(-ds_lsdp_c)*(c_fe + c_feoh2) ;
    K(E_fe+i*NEQ,I_c_fe+i*NEQ)      += volume[i]*phi*s_l ;
    K(E_fe+i*NEQ,I_c_feoh2+i*NEQ)   += volume[i]*phi*s_l ;
    /*
    CONSERVATION DE LA CHARGE Q = (H+) - (OH-) + 2(Fe2+) + (K+) - (A-)   :  0 + dt * div(i) = 0
    */

    /*
    CONSERVATION DE L'OXYGENE    :  (m_o21 - m_o2n) + dt * div(w_o21) = 0
    */
    K(E_o2+i*NEQ,I_p_l+i*NEQ)      += volume[i]*phi*(-ds_lsdp_c)*(c_o2_l - c_o2_g) ;
    K(E_o2+i*NEQ,I_c_o2+i*NEQ)     += volume[i]*phi*(s_l + s_g*dc_o2_gsdc_o2_l) ;
    /*
    CONSERVATION DES CATIONS    :  (m_cat1 - m_catn) + dt * div(w_cat1) = 0
    */
    K(E_cat+i*NEQ,I_p_l+i*NEQ)     += volume[i]*phi*(-ds_lsdp_c)*c_cat ;
    K(E_cat+i*NEQ,I_c_cat+i*NEQ)   += volume[i]*phi*s_l ;
    /*
    CONSERVATION DES ANIONS    :  (m_ani1 - m_anin) + dt * div(w_ani1) = 0
    */
    K(E_ani+i*NEQ,I_p_l+i*NEQ)     += volume[i]*phi*(-ds_lsdp_c)*c_ani ;
    K(E_ani+i*NEQ,I_c_ani+i*NEQ)   += volume[i]*phi*s_l ;
    /* 
    EQUILIBRE CHIMIQUE DE LA DISSOCIATION DE Fe(OH)2_l EN Fe2+ et 2OH- 
    */
    K(E_feoh2_l+i*NEQ,I_c_fe+i*NEQ)     += volume[i]*(-c_oh*c_oh/k_feoh2_l - 2*c_fe*c_oh*dc_ohsdc_fe[i]/k_feoh2_l) ;
    K(E_feoh2_l+i*NEQ,I_c_cat+i*NEQ)    += volume[i]*(-2*c_fe*c_oh*dc_ohsdc_cat[i]/k_feoh2_l) ;
    K(E_feoh2_l+i*NEQ,I_c_ani+i*NEQ)    += volume[i]*(-2*c_fe*c_oh*dc_ohsdc_ani[i]/k_feoh2_l) ;
    K(E_feoh2_l+i*NEQ,I_c_feoh2+i*NEQ)  += volume[i]*(dc_feoh2_lsdc_feoh2[i]) ;
  }
 /*
 termes d'ecoulement
 */
  tr        = dt*surf/dx ;

  trd_h     = tr*KD_H ;
  trd_oh    = tr*KD_OH ;
  trd_fe    = tr*KD_Fe ;
  trd_feoh2_l = tr*KD_FeOH2_l ;
  trd_o2    = tr*KD_O2 ;
  trd_cat   = tr*KD_cat ;
  trd_ani   = tr*KD_ani ;
  trd_m     = tr*KD_m ;

  trf_oh    = tr*KF_OH ;
  trf_h     = tr*KF_H ;
  trf_fe    = tr*KF_Fe ;
  trf_feoh2_l = tr*KF_FeOH2_l ;
  trf_o2_l  = tr*KF_O2_l ;
  trf_o2_g  = tr*KF_O2_g/k_hen ;
  trf_cat   = tr*KF_cat ;
  trf_ani   = tr*KF_ani ;
  trf_o2    = trf_o2_l + trf_o2_g ;

  tre_h     = tr*KE_H ;
  tre_oh    = tr*KE_OH ;
  tre_fe    = tr*KE_Fe ;
  tre_cat   = tr*KE_cat ;
  tre_ani   = tr*KE_ani ;

  /*
    CONSERVATION DE LA MASSE LIQUIDE + O2 GAZ
  */
  for(i=0;i<2;i++) {
    c[i] = trd_m ;
  }
  K(E_mass,I_p_l)            += + c[0] ;
  K(E_mass,I_p_l+NEQ)        += - c[1] ;
  K(E_mass+NEQ,I_p_l)        += - c[0] ;
  K(E_mass+NEQ,I_p_l+NEQ)    += + c[1] ;
  
  for(i=0;i<2;i++) {
    c[i] = M_o2*trf_o2_g ;
  }
  K(E_mass,I_c_o2)           += + c[0] ;
  K(E_mass,I_c_o2+NEQ)       += - c[1] ;
  K(E_mass+NEQ,I_c_o2)       += - c[0] ;
  K(E_mass+NEQ,I_c_o2+NEQ)   += + c[1] ;

  /*
    CONSERVATION DE L'ELEMENT Fe
  */
  for(i=0;i<2;i++) {
    c[i] = trd_fe + trd_feoh2_l ;
  }
  K(E_fe,I_p_l)              += + c[0] ;
  K(E_fe,I_p_l+NEQ)          += - c[1] ;
  K(E_fe+NEQ,I_p_l)          += - c[0] ;
  K(E_fe+NEQ,I_p_l+NEQ)      += + c[1] ;
  
  for(i=0;i<2;i++) {
    c[i] = trf_fe ;
  }
  K(E_fe,I_c_fe)             += + c[0] ;
  K(E_fe,I_c_fe+NEQ)         += - c[1] ;
  K(E_fe+NEQ,I_c_fe)         += - c[0] ;
  K(E_fe+NEQ,I_c_fe+NEQ)     += + c[1] ;
  
  K(E_fe,I_psi)              += + tre_fe ;
  K(E_fe,I_psi+NEQ)          += - tre_fe ;
  K(E_fe+NEQ,I_psi)          += - tre_fe ;
  K(E_fe+NEQ,I_psi+NEQ)      += + tre_fe ;
  
  for(i=0;i<2;i++) {
    c[i] = trf_feoh2_l*dc_feoh2_lsdc_feoh2[i] ;
  }
  K(E_fe,I_c_feoh2)          += + c[0] ;
  K(E_fe,I_c_feoh2+NEQ)      += - c[1] ;
  K(E_fe+NEQ,I_c_feoh2)      += - c[0] ;
  K(E_fe+NEQ,I_c_feoh2+NEQ)  += + c[1] ;

  /*
    CONSERVATION DE LA CHARGE    :  0 + dt * div(i) = 0
    i = j_h - j_oh + 2j_fe + j_cat - j_ani
  */
  for(i=0;i<2;i++) {
    c[i] = z_h*trf_h*dc_hsdc_fe[i] + z_oh*trf_oh*dc_ohsdc_fe[i] + z_fe*trf_fe ;
  }
  K(E_q,I_c_fe)             += + c[0] ;
  K(E_q,I_c_fe+NEQ)         += - c[1] ;
  K(E_q+NEQ,I_c_fe)         += - c[0] ;
  K(E_q+NEQ,I_c_fe+NEQ)     += + c[1] ;
  
  for(i=0;i<2;i++) {
    c[i] = z_h*trf_h*dc_hsdc_cat[i] + z_oh*trf_oh*dc_ohsdc_cat[i] + z_cat*trf_cat ;
  }
  K(E_q,I_c_cat)            += + c[0] ;
  K(E_q,I_c_cat+NEQ)        += - c[1] ;
  K(E_q+NEQ,I_c_cat)        += - c[0] ;
  K(E_q+NEQ,I_c_cat+NEQ)    += + c[1] ;
  
  for(i=0;i<2;i++) {
    c[i] = z_h*trf_h*dc_hsdc_ani[i] + z_oh*trf_oh*dc_ohsdc_ani[i] + z_ani*trf_ani ;
  }
  K(E_q,I_c_ani)            += + c[0] ;
  K(E_q,I_c_ani+NEQ)        += - c[1] ;
  K(E_q+NEQ,I_c_ani)        += - c[0] ;
  K(E_q+NEQ,I_c_ani+NEQ)    += + c[1] ;
  
  for(i=0;i<2;i++) {
    c[i] = z_h*tre_h + z_oh*tre_oh + z_fe*tre_fe + z_cat*tre_cat + z_ani*tre_ani ;
  }
  K(E_q,I_psi)              += + c[0] ;
  K(E_q,I_psi+NEQ)          += - c[1] ;
  K(E_q+NEQ,I_psi)          += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)      += + c[1] ;

  /*
    CONSERVATION DE L'OXYGENE O2
  */
  K(E_o2,I_p_l)             += + trd_o2 ;
  K(E_o2,I_p_l+NEQ)         += - trd_o2 ;
  K(E_o2+NEQ,I_p_l)         += - trd_o2 ;
  K(E_o2+NEQ,I_p_l+NEQ)     += + trd_o2 ;
  
  K(E_o2,I_c_o2)            += + trf_o2 ;
  K(E_o2,I_c_o2+NEQ)        += - trf_o2 ;
  K(E_o2+NEQ,I_c_o2)        += - trf_o2 ;
  K(E_o2+NEQ,I_c_o2+NEQ)    += + trf_o2 ;

  /*
    CONSERVATION DES CATIONS K+ 
  */
  K(E_cat,I_p_l)            += + trd_cat ;
  K(E_cat,I_p_l+NEQ)        += - trd_cat ;
  K(E_cat+NEQ,I_p_l)        += - trd_cat ;
  K(E_cat+NEQ,I_p_l+NEQ)    += + trd_cat ;
  
  K(E_cat,I_psi)            += + tre_cat ;
  K(E_cat,I_psi+NEQ)        += - tre_cat ;
  K(E_cat+NEQ,I_psi)        += - tre_cat ;
  K(E_cat+NEQ,I_psi+NEQ)    += + tre_cat ;
  
  K(E_cat,I_c_cat)          += + trf_cat ;
  K(E_cat,I_c_cat+NEQ)      += - trf_cat ;
  K(E_cat+NEQ,I_c_cat)      += - trf_cat ;
  K(E_cat+NEQ,I_c_cat+NEQ)  += + trf_cat ;

  /*
    CONSERVATION DES ANIONS A-
  */
  K(E_ani,I_p_l)            += + trd_ani ;
  K(E_ani,I_p_l+NEQ)        += - trd_ani ;
  K(E_ani+NEQ,I_p_l)        += - trd_ani ;
  K(E_ani+NEQ,I_p_l+NEQ)    += + trd_ani ;
  
  K(E_ani,I_psi)            += + tre_ani ;
  K(E_ani,I_psi+NEQ)        += - tre_ani ;
  K(E_ani+NEQ,I_psi)        += - tre_ani ;
  K(E_ani+NEQ,I_psi+NEQ)    += + tre_ani ;
  
  K(E_ani,I_c_ani)          += + trf_ani ;
  K(E_ani,I_c_ani+NEQ)      += - trf_ani ;
  K(E_ani+NEQ,I_c_ani)      += - trf_ani ;
  K(E_ani+NEQ,I_c_ani+NEQ)  += + trf_ani ;
  
  return(0) ;

#undef K
}

void rs22(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double c_fe,c_cat,c_ani,c_feoh2,c_oh,c_feoh2_l ;
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;
  
  if(el.dim < dim) return ;

  /* CALCUL DE volume ET DE surf */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)*0.5 ;
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ;
    }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;

  /* CONSERVATION DE LA MASSE LIQUIDE + O2 GAZ */
  R(0,E_mass) -= volume[0]*(M(0) - M_n(0)) + dt*surf*W_m ;
  R(1,E_mass) -= volume[1]*(M(1) - M_n(1)) - dt*surf*W_m ;
  /* CONSERVATION DE L'ELEMENT Fe = Fe2+ + Fe(OH)2 */
  R(0,E_fe)   -= volume[0]*(N_Fe(0) - N_Fen(0)) + dt*surf*W_Fe ;
  R(1,E_fe)   -= volume[1]*(N_Fe(1) - N_Fen(1)) - dt*surf*W_Fe ;
  /* CONSERVATION DE LA CHARGE */
  R(0,E_q)    -= + dt*surf*W_q ;
  R(1,E_q)    -= - dt*surf*W_q ;
  /* CONSERVATION DE L'OXYGENE O2 */
  R(0,E_o2)   -= volume[0]*(N_O2(0) - N_O2n(0)) + dt*surf*W_O2 ;
  R(1,E_o2)   -= volume[1]*(N_O2(1) - N_O2n(1)) - dt*surf*W_O2 ;
  /* CONSERVATION DES CATIONS K+  */
  R(0,E_cat)  -= volume[0]*(N_cat(0) - N_catn(0)) + dt*surf*W_cat ;
  R(1,E_cat)  -= volume[1]*(N_cat(1) - N_catn(1)) - dt*surf*W_cat ;
  /* CONSERVATION DES ANIONS A-   */
  R(0,E_ani)  -= volume[0]*(N_ani(0) - N_anin(0)) + dt*surf*W_ani ;
  R(1,E_ani)  -= volume[1]*(N_ani(1) - N_anin(1)) - dt*surf*W_ani ;
  /* EQUILIBRE CHIMIQUE DE LA DISSOCIATION DE Fe(OH)2_l EN Fe2+ et 2OH- */
  for(i=0;i<2;i++) {
    c_fe      = C_Fe(i) ;
    c_cat     = C_cat(i) ;
    c_ani     = C_ani(i) ;
    c_feoh2   = C_FeOH2(i) ;

    c_oh      = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
    c_feoh2_l = concentration_feoh2_l(c_feoh2,s_feoh2) ;
    R(i,E_feoh2_l) -= volume[i]*(c_feoh2_l - c_fe*c_oh*c_oh/k_feoh2_l) ;
  }

#undef R
}

int so22(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double p_l,p_c,s_l,psi ;
  double c_oh,c_h,c_fe,c_feoh2,c_feoh2_l,c_cat,c_ani,c_o2_l ;
  double r_elec ;
  double grd_p_l,grd_h,grd_oh,grd_fe,grd_feoh2_l,grd_o2,grd_cat,grd_ani,grd_psi,grd_feoh2 ;
  double w_h2o,w_h,w_oh,w_fe,w_feoh2_l,w_o2_l,w_o2_g,w_o2,w_cat,w_ani,i_ionic,i_ohmic ;
  double dx ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;

/* Donnees */
  phi     = el.mat->pr[pm("phi")] ;

  /* initialisation */
  nso = 21 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pression */
  p_l       = param(u,h_s,el.nn,I_p_l) ;
  p_c       = p_g - p_l ;
  /* saturation */
  s_l       = SATURATION(p_c,el.mat->cb[0]) ;
  /* concentrations */
  c_fe      = param(u,h_s,el.nn,I_c_fe) ;
  c_o2_l    = param(u,h_s,el.nn,I_c_o2) ;
  c_cat     = param(u,h_s,el.nn,I_c_cat) ;
  c_ani     = param(u,h_s,el.nn,I_c_ani) ;
  c_feoh2   = param(u,h_s,el.nn,I_c_feoh2) ;

  c_oh      = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
  c_h       = k_eau/c_oh ;
  c_feoh2_l = concentration_feoh2_l(c_feoh2,s_feoh2) ;
  /* resistivite electrique */
  r_elec    = 1./(FARADAY*(z_h*KE_H + z_fe*KE_Fe + z_cat*KE_cat + z_oh*KE_OH + z_ani*KE_ani)) ;
  /* potentiel electrique */
  psi       = param(u,h_s,el.nn,I_psi) ;

  /* gradients */
  dx        = x[1][0] - x[0][0] ;
  
  grd_p_l   = (P_l(1) - P_l(0))/dx ;
  grd_fe    = (C_Fe(1) - C_Fe(0))/dx ;
  grd_o2    = (C_O2(1) - C_O2(0))/dx ;
  grd_cat   = (C_cat(1) - C_cat(0))/dx ;
  grd_ani   = (C_ani(1) - C_ani(0))/dx ;
  grd_feoh2 = (C_FeOH2(1) - C_FeOH2(0))/dx ;
  grd_psi   = (PSI(1) - PSI(0))/dx ;

  {
    double c_fe,c_feoh2,c_cat,c_ani ;
    double c_oh[2],c_h[2],c_feoh2_l[2] ;
    for(i=0;i<2;i++) {
      c_fe      = C_Fe(i) ;
      c_cat     = C_cat(i) ;
      c_ani     = C_ani(i) ;
      c_feoh2   = C_FeOH2(i) ;
      
      c_oh[i]   = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
      c_h[i]    = k_eau/c_oh[i] ;
      c_feoh2_l[i] = concentration_feoh2_l(c_feoh2,s_feoh2) ;
    }
    grd_oh      = (c_oh[1] - c_oh[0])/dx ;
    grd_h       = (c_h[1] - c_h[0])/dx ;
    grd_feoh2_l = (c_feoh2_l[1] - c_feoh2_l[0])/dx ;
  }

  /* flux */
  w_h       = - KD_H*grd_p_l - KF_H*grd_h - KE_H*grd_psi ;
  w_oh      = - KD_OH*grd_p_l - KF_OH*grd_oh - KE_OH*grd_psi;
  w_fe      = - KD_Fe*grd_p_l - KF_Fe*grd_fe - KE_Fe*grd_psi ;
  w_feoh2_l = - KD_FeOH2_l*grd_p_l - KF_FeOH2_l*grd_feoh2_l ;
  w_o2_l    = - KD_O2*grd_p_l - KF_O2_l*grd_o2 ;
  w_o2_g    =                 - KF_O2_g/k_hen*grd_o2 ;
  w_o2      = w_o2_l + w_o2_g ;
  w_cat     = - KD_cat*grd_p_l - KF_cat*grd_cat - KE_cat*grd_psi ;
  w_ani     = - KD_ani*grd_p_l - KF_ani*grd_ani - KE_ani*grd_psi ;
  i_ionic   = FARADAY*(z_fe*w_fe + z_h*w_h + z_cat*w_cat + z_oh*w_oh + z_ani*w_ani) ;
  i_ohmic   = - grd_psi/r_elec ;
  w_h2o     = - KD_m*grd_p_l - (M_o2*w_o2_l + M_h*w_h + M_oh*w_oh + M_fe*w_fe + M_feoh2*w_feoh2_l + M_cat*w_cat + M_ani*w_ani) ;
  w_h2o    /= M_h2o ;

  /* quantites exploitees */
  strcpy(r[0].text,"pression_liquide") ; r[0].n = 1 ;
  r[0].v[0] = p_l ;
  strcpy(r[1].text,"saturation") ; r[1].n = 1 ;
  r[1].v[0] = s_l ;
  strcpy(r[2].text,"molarite_O2") ; r[2].n = 1 ;
  r[2].v[0] = c_o2_l ;
  strcpy(r[3].text,"molarite_H+") ; r[3].n = 1 ;
  r[3].v[0] = c_h ;
  strcpy(r[4].text,"molarite_OH-") ; r[4].n = 1 ;
  r[4].v[0] = c_oh ;
  strcpy(r[5].text,"molarite_Fe2+") ; r[5].n = 1 ;
  r[5].v[0] = c_fe ;
  strcpy(r[6].text,"molarite_Fe(OH)2") ; r[6].n = 1 ;
  r[6].v[0] = c_feoh2 ;
  strcpy(r[7].text,"molarite_cations") ; r[7].n = 1 ;
  r[7].v[0] = c_cat ;
  strcpy(r[8].text,"molarite_anions") ; r[8].n = 1 ;
  r[8].v[0] = c_ani ;
  strcpy(r[9].text,"potentiel_electrique") ; r[9].n = 1 ;
  r[9].v[0] = psi ;
  strcpy(r[10].text,"flux_H2O") ; r[10].n = 3 ;
  r[10].v[0] = w_h2o ;
  strcpy(r[11].text,"flux_H+") ; r[11].n = 3 ;
  r[11].v[0] = w_h ;
  strcpy(r[12].text,"flux_OH-") ; r[12].n = 3 ;
  r[12].v[0] = w_oh ;
  strcpy(r[13].text,"flux_Fe2+") ; r[13].n = 3 ;
  r[13].v[0] = w_fe ;
  strcpy(r[14].text,"flux_Fe(OH)2") ; r[14].n = 3 ;
  r[14].v[0] = w_feoh2_l ;
  strcpy(r[15].text,"flux_O2") ; r[15].n = 3 ;
  r[15].v[0] = w_o2 ;
  strcpy(r[16].text,"flux_cations") ; r[16].n = 3 ;
  r[16].v[0] = w_cat ;
  strcpy(r[17].text,"flux_anions") ; r[17].n = 3 ;
  r[17].v[0] = w_ani ;
  strcpy(r[18].text,"courant_ionique") ; r[18].n = 3 ;
  r[18].v[0] = i_ionic ;
  strcpy(r[19].text,"courant_ohmique") ; r[19].n = 3 ;
  r[19].v[0] = i_ohmic ;
  strcpy(r[20].text,"resistivite_electrique") ; r[20].n = 1 ;
  r[20].v[0] = r_elec ;
  return (nso) ;
}



void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Calcul des flux (f) */ 
{
  double c_oh[2],c_h[2],c_fe,c_feoh2,c_feoh2_l[2],c_cat,c_ani ;
  double grd_p_l,grd_h,grd_oh,grd_fe,grd_feoh2_l,grd_o2,grd_cat,grd_ani,grd_psi,grd_feoh2 ;
  double w_fe,w_feoh2,w_o2_l,w_o2_g ;
  double dx ;
  int    i ;

  /* gradients */
  dx      = x[1][0] - x[0][0] ;
  
  grd_p_l   = (P_l(1) - P_l(0))/dx ;
  grd_fe    = (C_Fe(1) - C_Fe(0))/dx ;
  grd_o2    = (C_O2(1) - C_O2(0))/dx ;
  grd_cat   = (C_cat(1) - C_cat(0))/dx ;
  grd_ani   = (C_ani(1) - C_ani(0))/dx ;
  grd_feoh2 = (C_FeOH2(1) - C_FeOH2(0))/dx ;
  grd_psi   = (PSI(1) - PSI(0))/dx ;

  for(i=0;i<2;i++) {
    c_fe      = C_Fe(i) ;
    c_cat     = C_cat(i) ;
    c_ani     = C_ani(i) ;
    c_feoh2   = C_FeOH2(i) ;

    c_oh[i]   = concentration_oh(z_fe*c_fe + z_cat*c_cat + z_ani*c_ani,k_eau) ;
    c_h[i]    = k_eau/c_oh[i] ;
    c_feoh2_l[i] = concentration_feoh2_l(c_feoh2,s_feoh2) ;
  }
  grd_oh      = (c_oh[1] - c_oh[0])/dx ;
  grd_h       = (c_h[1] - c_h[0])/dx ;
  grd_feoh2_l = (c_feoh2_l[1] - c_feoh2_l[0])/dx ;


  /* flux */
  w_fe    = - KD_Fe*grd_p_l      - KF_Fe*grd_fe           - KE_Fe*grd_psi ;
  w_feoh2 = - KD_FeOH2_l*grd_p_l - KF_FeOH2_l*grd_feoh2_l ;
  w_o2_l  = - KD_O2*grd_p_l      - KF_O2_l*grd_o2 ;
  w_o2_g  =                      - KF_O2_g/k_hen*grd_o2 ;

  W_O2    = w_o2_l + w_o2_g ;
  W_Fe    = w_fe + w_feoh2 ;
  W_cat   = - KD_cat*grd_p_l - KF_cat*grd_cat - KE_cat*grd_psi ;
  W_ani   = - KD_ani*grd_p_l - KF_ani*grd_ani - KE_ani*grd_psi ;
  W_q     = - z_h*KF_H*grd_h - z_oh*KF_OH*grd_oh - z_fe*KF_Fe*grd_fe - z_cat*KF_cat*grd_cat - z_ani*KF_ani*grd_ani
         -   (z_h*KE_H + z_oh*KE_OH + z_fe*KE_Fe + z_cat*KE_cat + z_ani*KE_ani)*grd_psi ;
  W_m     = - KD_m*grd_p_l + M_o2*w_o2_g ;
} 


#undef P_l
#undef C_Fe
#undef PSI
#undef C_O2
#undef C_cat
#undef C_ani

#undef P_ln
#undef C_Fen
#undef PSIn
#undef C_O2n
#undef C_catn
#undef C_anin

#undef M
#undef N_Fe
#undef N_O2
#undef N_cat
#undef N_ani
#undef W_m
#undef W_Fe
#undef W_O2
#undef W_cat
#undef W_ani

#undef M_n
#undef N_Fen
#undef N_O2n
#undef N_catn
#undef N_anin

#undef KD_m
#undef KD_H
#undef KD_OH
#undef KD_Fe
#undef KD_FeOH2_l
#undef KD_O2
#undef KD_cat
#undef KD_ani

#undef KF_OH
#undef KF_H
#undef KF_Fe
#undef KF_FeOH2_l
#undef KF_O2_l
#undef KF_O2_g
#undef KF_cat
#undef KF_ani

#undef KE_H
#undef KE_OH
#undef KE_Fe
#undef KE_cat
#undef KE_ani



#undef z_fe
#undef z_oh
#undef z_h
#undef z_cat
#undef z_ani

#undef M_fe
#undef M_o
#undef M_oh
#undef M_h
#undef M_h2o
#undef M_o2
#undef M_cat
#undef M_ani
#undef M_feoh2

#undef d_oh
#undef d_h
#undef d_o2_l
#undef d_o2_g
#undef d_fe
#undef d_feoh2_l
#undef d_cat
#undef d_ani

#undef k_eau
#undef k_feoh2_l
#undef s_feoh2

#undef FARADAY
#undef RT 

#undef mu_l
#undef c_h2o



double tortuosite_l(double p_c,mate_t *mat)
{
  double phi = mat->pr[pm("phi")] ;
  double tau_l_sat = 0.296e-3*exp(9.95*phi)/phi ;
  return(tau_l_sat*courbe(p_c,mat->cb[2])) ;
}

double tortuosite_g(double p_c,mate_t *mat)
{
  double phi = mat->pr[pm("phi")] ;
  double tau_g_sat = pow(phi,1.74) ;
  return(tau_g_sat*courbe(p_c,mat->cb[3])) ;
}

double concentration_oh(double q,double k_eau)
{
  return (0.5*(q + sqrt(q*q + 4*k_eau))) ;
}

double concentration_feoh2_l(double c_feoh2,double s_feoh2)
{
  return((c_feoh2 < s_feoh2) ? c_feoh2 : s_feoh2) ;
}

double saturation(double pc,crbe_t cb)
/* Degre de saturation prolongee */
{
  double sl ;
  int    ni = cb.n - 1 ;
  double pc1 = cb.a[0],pc2 = cb.a[1] ;
  double sl1 = cb.f[0],sl2 = cb.f[ni] ;
  
  if(pc < pc1) {
    if(sl1 >= 1.) sl = 1. ;
    else {
      double b  = - dcourbe(pc1,cb)/(1. - sl1) ;
      sl = 1. - (1. - sl1) * exp(b*(pc - pc1)) ;
    }
  } else if(pc > pc2) {
    if(sl2 <= 0.) sl = 0. ;
    else {
      double b  = - dcourbe(pc2,cb)/sl2 ;
      sl =  sl2 * exp(-b*(pc - pc2)) ;
    }
  } else {
    sl = courbe(pc,cb) ;
  }
  return(sl) ;
}

double dsaturation(double pc,crbe_t cb)
{
  double dsl ;
  int    ni = cb.n - 1 ;
  double pc1 = cb.a[0],pc2 = cb.a[1] ;
  double sl1 = cb.f[0],sl2 = cb.f[ni] ;
  
  if(pc < pc1) {
    if(sl1 >= 1.) dsl = 0. ;
    else {
      double b  = - dcourbe(pc1,cb)/(1. - sl1) ;
      dsl = - b * (1. - sl1) * exp(b*(pc - pc1)) ;
    }
  } else if(pc > pc2) {
    if(sl2 <= 0.) dsl = 0. ;
    else {
      double b  = - dcourbe(pc2,cb)/sl2 ;
      dsl =  - b * sl2 * exp(-b*(pc - pc2)) ;
    }
  } else {
    dsl = dcourbe(pc,cb) ;
  }
  return(dsl) ;
}

#undef SATURATION
#undef DSATURATION
