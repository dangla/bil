#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  70
#define TITLE   "Sulfuric acid attack of concrete (06/2011)"
#define AUTHORS "Yuan"

#include "OldMethods.h"

/* Macros */
#define NEQ     (3)
#define NVE     (16)
#define NVI     (13)
#define NVE_TR  (12)

#define E_S     (0)
#define E_q     (1)
#define E_Ca    (2)

#define I_H2SO4 (0)
#define I_psi   (1)
#define I_Ca_S  (2)

#define RHO      1
#define LOG_RHO  2
#define ZRHO     3
#define LOG_ZRHO 4
#define Ln10    2.302585093
#define U_H2SO4 LOG_RHO

#if (U_H2SO4 == LOG_RHO)
  #define ZC_H2SO4(n)  (exp(Ln10*u[(n)][I_H2SO4])/c_h2so4_eq)
  #define ZC_H2SO4n(n) (exp(Ln10*u_n[(n)][I_H2SO4])/c_h2so4_eq)
#elif (U_H2SO4 == LOG_ZRHO)
  #define ZC_H2SO4(n)  (exp(Ln10*u[(n)][I_H2SO4]))
  #define ZC_H2SO4n(n) (exp(Ln10*u_n[(n)][I_H2SO4]))
#elif (U_H2SO4 == ZRHO)
  #define ZC_H2SO4(n)  (u[(n)][I_H2SO4])
  #define ZC_H2SO4n(n) (u_n[(n)][I_H2SO4])
#else
  #define ZC_H2SO4(n)  (u[(n)][I_H2SO4]/c_h2so4_eq)
  #define ZC_H2SO4n(n) (u_n[(n)][I_H2SO4]/c_h2so4_eq)
#endif

#define ZN_Ca_S(n)  (u[(n)][I_Ca_S])
#define PSI(n)      (u[(n)][I_psi])

#define ZN_Ca_Sn(n) (u_n[(n)][I_Ca_S])
#define PSIn(n)     (u_n[(n)][I_psi])

#define N_S(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define W_S        (f[6])
#define W_q        (f[7])
#define W_Ca       (f[8])
#define N_CH(n)    (f[(9+n)])
#define N_CSH2(n)  (f[(11+n)])

#define N_Sn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define W_Sn       (f_n[6])
#define W_qn       (f_n[7])
#define W_Can      (f_n[8])
#define N_CHn(n)   (f_n[(9+n)])
#define N_CSH2n(n) (f_n[(11+n)])


#define KF_OH       (va[(0)])
#define KF_H        (va[(1)])
#define KF_SO4      (va[(2)])
#define KF_HSO4     (va[(3)])
#define KF_Ca       (va[(4)])
#define KF_H2SO4    (va[(5)])

#define Kpsi_OH     (va[(6)])
#define Kpsi_H      (va[(7)])
#define Kpsi_SO4    (va[(8)])
#define Kpsi_HSO4   (va[(9)])
#define Kpsi_Ca     (va[(10)])
#define Kpsi_q      (va[(11)])

#define N_CH0(n)    (va[(12+n)])
#define N_CSH20(n)  (va[(14+n)])

/*
  Aqueous Solution
*/

/* valence numbers */
#define z_ca       (2.)
#define z_h        (1.)
#define z_oh       (-1.)
#define z_so4      (-2.)
#define z_hso4     (-1.)

/* partial molar volumes of ions (dm3/mole) */
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_so4      (22.9e-3)
#define v_hso4     (0)
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */

/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s  */
#define d_oh       (1.22e-7)    /* 1.22e-7 (radius = 1.75e-10 m) */
#define d_h        (9.310e-7)    /* 4.76e-8 (radius = 4.5e-10 m) */
#define d_so4      (1.065e-7)     /*  */
#define d_hso4     (1.385e-7)
#define d_h2so4    (1.5e-7)
#define d_ca       (7.92e-8)

/* Porosity */
#define PHII(n)     (phi0 + V_CH*(N_CH0(n) - n_ch) + V_CSH2*(N_CSH20(n) - n_csh2))
#define PHI(n)      (phi0) /* (MAX(PHII(n),0.01)) */
#define DPHISDN_CH(n)     (0.) /* ((PHII(n) > 0.01) ? -V_CH : 0.) */
#define DPHISDN_CSH2(n)   (0.) /* ((PHII(n) > 0.01) ? -V_CSH2 : 0.) */
/* tortuosity */
#define  EXP1(x)           ((x > 1.) ? exp(x) : x*exp(1.))
#define  TORTUOSITY(phi)   (2.9e-4*EXP1(9.95*phi))

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_h2o      (1.e-14)          /* autoprotolyse de l'eau */

#define K_h2so4    (1.e3)        /* H2SO4[0] = HSO4[-] + H[+] */
#define K_hso4     (1.e-2)       /* HSO4[-]  = SO4[2-] + H[+] */

/*
  Solids
  CH  = Portlandite
  CSH2 = Gypsum
*/

/* CH */
/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* CH = Ca[2+] + 2OH[-] */
/* volumes molaires solides (dm3/mole) */
#define V_CH       (33.e-3)     /* (33.e-3) */
/* CSH2 */
/* Equilibrium Constant (ref = 1 mol/L) */
#define K_CSH2     (2.5e-5)          /* CSH2 = Ca[2+] + SO4[2-] + 2H2O */
/* Molar Volume  (dm3/mole) */
#define V_CSH2     (75.e-3)

#define C_H2SO4_eq   (K_h2o*K_h2o*K_CSH2/(K_h2so4*K_hso4*K_CH))

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */

/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,double) ;
static double poly4(double,double,double,double,double) ;
static double poly4_sansc(double,double,double,double) ;
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   concentrations(double,double,double*) ;

#define MIN(a,b)   ((a < b) ? a : b)
#define MAX(a,b)   ((a > b) ? a : b)
#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)
#define DNEGEXP(y) ((y < 0.) ? exp(y) : 0.)

/* Ion Activity Products */
#define IAP_CSH2(zc_h2so4,zn_ca_s)         (K_CSH2*NEGEXP(zn_ca_s)*MIN(zc_h2so4,1.))
#define IAP_CH(zc_h2so4,zn_ca_s)           (K_CH*NEGEXP(zn_ca_s)/MAX(zc_h2so4,1.))

/* Parametres */
static double phi0,c_h2so4_eq,t_ch,t_csh2 ;
static double n_ca_ref ;
#define  GET(a)     (el.mat->pr[pm(a)])

static double var[6] ;  
#define  c_oh       (var[0])
#define  c_h        (var[1])
#define  c_hso4     (var[2])
#define  c_so4      (var[3])
#define  c_ca       (var[4])
#define  c_h2so4    (var[5])

#define  CONCENTRATIONS(zc_h2so4,zn_ca_s)     concentrations(zc_h2so4,zn_ca_s,var)


int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"C_H2SO4_eq") == 0) return (2) ;
  else if(strcmp(s,"T_CH") == 0)     return (3) ;
  else if(strcmp(s,"T_CSH2") == 0)     return (4) ;
  else return(-1) ;
}


int dm70(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 5 ;
  
  mat->neq      = NEQ ;

  strcpy(mat->eqn[E_S],   "sulfur") ;
  strcpy(mat->eqn[E_Ca],  "calcium") ;
  strcpy(mat->eqn[E_q],   "charge") ;
   

#if (U_H2SO4 == LOG_RHO)
  strcpy(mat->inc[I_H2SO4],"logc_h2so4") ;
#elif (U_H2SO4 == LOG_ZRHO)
  strcpy(mat->inc[I_H2SO4],"logz_h2so4") ;
#elif (U_H2SO4 == ZRHO)
  strcpy(mat->inc[I_H2SO4],"z_h2so4") ;
#else
  strcpy(mat->inc[I_H2SO4],"c_h2so4") ;
#endif
  strcpy(mat->inc[I_Ca_S], "z_ca") ;
  strcpy(mat->inc[I_psi],  "psi") ;

  {
    /* initialisation automatique */
    double t_ch        = 600. ;
    double t_csh2      = 0. ;
    double n_ca_ref    = 1. ;

    mat->pr[pm("N_CH")]   = n_ca_ref ;
    mat->pr[pm("T_CH")]   = t_ch ;
    mat->pr[pm("T_CSH2")] = t_csh2 ;

    dmat(mat,ficd,pm,n_donnees) ;

    t_ch      = mat->pr[pm("T_CH")] ;
    t_csh2    = mat->pr[pm("T_CSH2")] ;

    if(t_csh2  == 0.) mat->pr[pm("T_CSH2")]  = t_ch ;

    mat->pr[pm("C_H2SO4_eq")] = C_H2SO4_eq ;
  }
  
  return(mat->n) ;
}


int qm70(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;
  printf("The system is composed of 3 equations:\n") ;
#if (U_H2SO4 == LOG_RHO)
  printf("\t- Mass balance of S       (logc_h2so4)\n") ;
#elif (U_H2SO4 == LOG_ZRHO)
  printf("\t- Mass balance of S       (logz_h2so4)\n") ;
#elif (U_H2SO4 == ZRHO)
  printf("\t- Mass balance of S       (z_h2so4)\n") ;
#else
  printf("\t- Mass balance of S       (c_h2so4)\n") ;
#endif
  printf("\t- Mass balance of charge  (psi)\n") ;
  printf("\t- Mass balance of Ca      (zn_ca_s)\n") ;

  printf("\n\
Take care of unities : \n\
\t lenght : dm !\n\
\t time   : s !\n") ;

  printf("Example of material data\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1       # Contenu en CH (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5      # Cinetique de dissolution de CH (s)\n") ;

  return(NEQ) ;
}

void tb70(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch70(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in70(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;
  
  /*
    Donnees
  */
  phi0        = GET("porosite") ;
  n_ca_ref    = GET("N_CH") ;
  c_h2so4_eq  = GET("C_H2SO4_eq") ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarities */
    double zc_h2so4   = ZC_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
  
    CONCENTRATIONS(zc_h2so4,zn_ca_s) ;
    {
    /* equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_csh2_eq  = n_ca_ref*MAX(zn_ca_s,0.) ;
    
    /* solid contents : CH, CSH2 */
    double n_ch       = (zc_h2so4 <= 1) ? n_ch_eq  : 0 ;
    double n_csh2     = (zc_h2so4 >  1) ? n_csh2_eq  : 0 ;
    double n_ca_s     = n_ch + n_csh2 ;
    double n_s_s      = n_csh2 ;

    /* porosity */
    double phi = phi0 ;

    /* molar contents */
    double n_so4      = phi*c_so4 ;
    double n_hso4     = phi*c_hso4 ;
    double n_h2so4    = phi*c_h2so4 ;
    double n_ca       = phi*c_ca ;
    
    N_S(i)  = n_so4 + n_hso4 + n_h2so4 + n_s_s ;
    N_Ca(i) = n_ca + n_ca_s ;

    /* charge density */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hso4*c_hso4 + z_so4*c_so4 + z_ca*c_ca ;

    N_CH(i)    = n_ch ;
    N_CSH2(i)  = n_csh2 ;

    N_CH0(i)   = n_ch ;
    N_CSH20(i) = n_csh2 ;
    }
  }
  
  if(el.dim < dim) return ;

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex70(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  
  if(el.dim < dim) return(0) ;
  
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct70(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;
  
  /*
    Donnees
  */
  phi0        = GET("porosite") ;
  n_ca_ref    = GET("N_CH") ;
  c_h2so4_eq  = GET("C_H2SO4_eq") ;
  t_ch        = GET("T_CH") ;
  t_csh2      = GET("T_CSH2") ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarities */
    double zc_h2so4   = ZC_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
  
    CONCENTRATIONS(zc_h2so4,zn_ca_s) ;
    {
    /* kinetics */
    double n_chn      = N_CHn(i) ;
    double n_csh2n    = N_CSH2n(i) ;
    double n_ch_ci    = n_chn*pow(zc_h2so4,-dt/t_ch) ;    /* if zc_h2so4 > 1 */
    double n_csh2_ci  = n_csh2n*pow(zc_h2so4,dt/t_csh2) ; /* if zc_h2so4 < 1 */
    
    /* equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_csh2_eq  = n_ca_ref*MAX(zn_ca_s,0.) ;
    
    /* solid contents : CH, CSH2 */
    double n_ch       = (zc_h2so4 <= 1) ? n_ch_eq   : n_ch_ci ;
    double n_csh2     = (zc_h2so4 >  1) ? n_csh2_eq : n_csh2_ci ;
    double n_ca_s     = n_ch + n_csh2 ;
    double n_s_s      = n_csh2 ;

    /* porosity */
    double phi        = PHI(i) ;

    /* molar contents */
    double n_so4      = phi*c_so4 ;
    double n_hso4     = phi*c_hso4 ;
    double n_h2so4    = phi*c_h2so4 ;
    double n_ca       = phi*c_ca ;
    
    N_S(i)  = n_so4 + n_hso4 + n_h2so4 + n_s_s ;
    N_Ca(i) = n_ca + n_ca_s ;

    /* charge density */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hso4*c_hso4 + z_so4*c_so4 + z_ca*c_ca ;

    N_CH(i)    = n_ch ;
    N_CSH2(i)  = n_csh2 ;

    if(c_h2so4 < 0. || n_ca_s < 0. || n_s_s < 0. || c_oh < 0.) {
      printf("x         = %e\n",x[i][0]) ;
      printf("c_h2so4   = %e\n",c_h2so4) ;
      printf("n_csh2    = %e\n",n_csh2) ;
      printf("n_ca_s    = %e\n",n_ca_s) ;
      printf("n_s_s     = %e\n",n_s_s) ;
      printf("zn_ca_s   = %e\n",zn_ca_s) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
    if(phi < 0.) {
      printf("phi = %e\n",phi) ;
      return(-1) ;
    }
    }
  }
  
  if(el.dim < dim) return(0) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
  
  return(0) ;
}

int mx70(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double c[2] ;

  double Dc_hSDzc_h2so4[2]        ;
  double Dc_ohSDzc_h2so4[2]       ;
  double Dc_h2so4SDzc_h2so4[2]    ;
  double Dc_hso4SDzc_h2so4[2]     ;
  double Dc_so4SDzc_h2so4[2]      ;
  double Dc_caSDzc_h2so4[2]       ;

  double Dc_hSDzn_ca_s[2]      ;
  double Dc_ohSDzn_ca_s[2]     ;
  double Dc_h2so4SDzn_ca_s[2]  ;
  double Dc_hso4SDzn_ca_s[2]   ;
  double Dc_so4SDzn_ca_s[2]    ;
  double Dc_caSDzn_ca_s[2]     ;
  
  /*
    Donnees
  */
  phi0        = GET("porosite") ;
  n_ca_ref    = GET("N_CH") ;
  c_h2so4_eq  = GET("C_H2SO4_eq") ;
  t_ch        = GET("T_CH") ;
  t_csh2      = GET("T_CSH2") ;

  /*
    Initialisation 
  */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    CALCUL DE volume ET DE surf 
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<el.nn;i++) {
    volume[i] = fabs(dx)*0.5 ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = 2*M_PI*xm ; else surf = 1 ;
  /*
    termes d'accumulation
  */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double zc_h2so4   = ZC_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double P_CSH2     = IAP_CSH2(zc_h2so4,zn_ca_s) ;
  
    CONCENTRATIONS(zc_h2so4,zn_ca_s) ;
    
    {
    /* solid contents : CH, CSH2 */
    double n_ch       = N_CH(i) ;
    double n_csh2     = N_CSH2(i) ;
    
    /* kinetics */
    double n_chn      = N_CHn(i) ;
    double n_csh2n    = N_CSH2n(i) ;
    double n_ch_ci    = n_chn*pow(zc_h2so4,-dt/t_ch) ;    /* if zc_h2so4 > 1 */
    double n_csh2_ci  = n_csh2n*pow(zc_h2so4,dt/t_csh2) ; /* if zc_h2so4 < 1 */

    /* porosity */
    double phi        = PHI(i) ;
    double dphisdn_ch   = DPHISDN_CH(i) ;
    double dphisdn_csh2 = DPHISDN_CSH2(i) ;

    /* derivatives ... */
    /* ... with respect to zc_h2so4 */
    double dzc_h2so4          = 1.e-4*zc_h2so4*((zc_h2so4 > ZC_H2SO4n(i)) ? 1 : -1) ;
    double zc_h2so42          = zc_h2so4 + dzc_h2so4 ;
    double dc_h2so4sdzc_h2so4 = c_h2so4_eq ;
    double c_oh2              = concentration_oh(zc_h2so42,zn_ca_s) ;   
    double dc_ohsdzc_h2so4    = (c_oh2 - c_oh)/dzc_h2so4 ;
    double dc_hsdzc_h2so4     = - c_h*dc_ohsdzc_h2so4/c_oh ;
    double dc_hso4sdzc_h2so4  = c_hso4*(dc_h2so4sdzc_h2so4/c_h2so4 \
                              + dc_ohsdzc_h2so4/c_oh) ;
    double dc_so4sdzc_h2so4   = c_so4*(dc_h2so4sdzc_h2so4/c_h2so4 \
                              + 2*dc_ohsdzc_h2so4/c_oh) ;
			                     
    double dP_CSH2sdzc_h2so4  = (zc_h2so4 < 1.) ? P_CSH2/zc_h2so4 : 0. ;
    
    double dc_casdzc_h2so4    = c_ca*(dP_CSH2sdzc_h2so4/P_CSH2 - dc_so4sdzc_h2so4/c_so4) ;
		                    
    double dn_ch_cisdzc_h2so4    = n_ch_ci*(-dt/t_ch)/zc_h2so4 ;
    double dn_csh2_cisdzc_h2so4  = n_csh2_ci*(dt/t_csh2)/zc_h2so4 ;

    double dn_ch_eqsdzc_h2so4    = 0 ;
    double dn_csh2_eqsdzc_h2so4  = 0 ;

    double dn_chsdzc_h2so4     = (zc_h2so4 <= 1) ? dn_ch_eqsdzc_h2so4   : dn_ch_cisdzc_h2so4 ;
    double dn_csh2sdzc_h2so4   = (zc_h2so4 >  1) ? dn_csh2_eqsdzc_h2so4 : dn_csh2_cisdzc_h2so4 ;
    
    double dn_ca_ssdzc_h2so4   = dn_chsdzc_h2so4 + dn_csh2sdzc_h2so4 ;
    double dn_s_ssdzc_h2so4    = dn_csh2sdzc_h2so4 ;

    double dphisdzc_h2so4      = dphisdn_ch*dn_chsdzc_h2so4 \
                               + dphisdn_csh2*dn_csh2sdzc_h2so4 ;

    double dn_h2so4sdzc_h2so4  = phi*dc_h2so4sdzc_h2so4 + dphisdzc_h2so4*c_h2so4 ;
    double dn_hso4sdzc_h2so4   = phi*dc_hso4sdzc_h2so4 + dphisdzc_h2so4*c_hso4 ;
    double dn_so4sdzc_h2so4    = phi*dc_so4sdzc_h2so4 + dphisdzc_h2so4*c_so4 ;
    double dn_casdzc_h2so4     = phi*dc_casdzc_h2so4 + dphisdzc_h2so4*c_ca ;
    
    /* with respect to zn_ca_s */
    /* double dzn_ca_s             = ((zn_ca_s > 0.) ? 1 : -1)*1.e-2 ; */
    double dzn_ca_s            = 1.e-6*((zn_ca_s > ZN_Ca_Sn(i)) ? 1 : -1) ;
    double zn_ca_s2            = zn_ca_s + dzn_ca_s ;
    double c_oh3               = concentration_oh(zc_h2so4,zn_ca_s2) ;
    double dc_ohsdzn_ca_s      = (c_oh3 - c_oh)/dzn_ca_s ;
    double dc_hsdzn_ca_s       = - c_h*dc_ohsdzn_ca_s/c_oh ;
    double dc_hso4sdzn_ca_s    = c_hso4*(dc_ohsdzn_ca_s/c_oh) ;
    double dc_so4sdzn_ca_s     = c_so4*2*dc_ohsdzn_ca_s/c_oh ;

    double dP_CSH2sdzn_ca_s    = (zn_ca_s < 0) ? P_CSH2 : 0 ;
    
    double dc_casdzn_ca_s      = c_ca*(dP_CSH2sdzn_ca_s/P_CSH2 - dc_so4sdzn_ca_s/c_so4) ;

    double dn_ch_cisdzn_ca_s   = 0 ;
    double dn_csh2_cisdzn_ca_s = 0 ;


    double dn_ch_eqsdzn_ca_s   = (zn_ca_s > 0) ? n_ca_ref : 0 ;
    double dn_csh2_eqsdzn_ca_s = (zn_ca_s > 0) ? n_ca_ref : 0 ;
     
    double dn_chsdzn_ca_s      = (zc_h2so4 <= 1) ? dn_ch_eqsdzn_ca_s   : dn_ch_cisdzn_ca_s ;
    double dn_csh2sdzn_ca_s    = (zc_h2so4 >  1) ? dn_csh2_eqsdzn_ca_s : dn_csh2_cisdzn_ca_s ;
    
    double dn_ca_ssdzn_ca_s    = dn_chsdzn_ca_s + dn_csh2sdzn_ca_s ;
    double dn_s_ssdzn_ca_s     = dn_csh2sdzn_ca_s ;

    double dphisdzn_ca_s       = dphisdn_ch*dn_chsdzn_ca_s \
                               + dphisdn_csh2*dn_csh2sdzn_ca_s ;

    double dn_h2so4sdzn_ca_s   = dphisdzn_ca_s*c_h2so4 ;
    double dn_hso4sdzn_ca_s    = phi*dc_hso4sdzn_ca_s + dphisdzn_ca_s*c_hso4 ;
    double dn_so4sdzn_ca_s     = phi*dc_so4sdzn_ca_s + dphisdzn_ca_s*c_so4 ;
    double dn_casdzn_ca_s      = phi*dc_casdzn_ca_s + dphisdzn_ca_s*c_ca ;

    j = i*NEQ ;
    /*
      Conservation de S (sulfure) : (n_S1 - n_Sn) + dt * div(w_S) = 0
    */
    K(E_S+j,I_H2SO4+j)  += volume[i]*(dn_h2so4sdzc_h2so4 + dn_hso4sdzc_h2so4 \
                        + dn_so4sdzc_h2so4 + dn_s_ssdzc_h2so4) ;
    K(E_S+j,I_Ca_S+j)   += volume[i]*(dn_h2so4sdzn_ca_s + dn_hso4sdzn_ca_s \
                        + dn_so4sdzn_ca_s + dn_s_ssdzn_ca_s) ;
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_H2SO4+j) += volume[i]*(dn_casdzc_h2so4 + dn_ca_ssdzc_h2so4) ;
    K(E_Ca+j,I_Ca_S+j)  += volume[i]*(dn_casdzn_ca_s + dn_ca_ssdzn_ca_s) ;
    /*
      Conservation de la charge  : div(w_q) = 0
    */

    /* sauvegardes pour les termes de transport */
    Dc_hSDzc_h2so4[i]        = dc_hsdzc_h2so4 ;
    Dc_ohSDzc_h2so4[i]       = dc_ohsdzc_h2so4 ;
    Dc_h2so4SDzc_h2so4[i]    = dc_h2so4sdzc_h2so4 ;
    Dc_hso4SDzc_h2so4[i]     = dc_hso4sdzc_h2so4 ;
    Dc_so4SDzc_h2so4[i]      = dc_so4sdzc_h2so4 ;
    Dc_caSDzc_h2so4[i]       = dc_casdzc_h2so4 ;

    Dc_hSDzn_ca_s[i]        = dc_hsdzn_ca_s ;
    Dc_ohSDzn_ca_s[i]       = dc_ohsdzn_ca_s ;
    Dc_h2so4SDzn_ca_s[i]    = 0. ;
    Dc_hso4SDzn_ca_s[i]     = dc_hso4sdzn_ca_s ;
    Dc_so4SDzn_ca_s[i]      = dc_so4sdzn_ca_s ;
    Dc_caSDzn_ca_s[i]       = dc_casdzn_ca_s ; 
    }
  }

  /* termes d'ecoulement */
  {
  double tr        = dt*surf/dx ;

  double trf_oh       = tr*KF_OH ;
  double trf_h        = tr*KF_H ;
  double trf_h2so4    = tr*KF_H2SO4 ;
  double trf_hso4     = tr*KF_HSO4 ;
  double trf_so4      = tr*KF_SO4 ;
  double trf_ca       = tr*KF_Ca ;
  
  double tre_hso4     = tr*Kpsi_HSO4 ;
  double tre_so4      = tr*Kpsi_SO4 ;
  double tre_ca       = tr*Kpsi_Ca ; 

  double tre_q        = tr*Kpsi_q ;
  
  /*
    Conservation de S (sulfure) : (n_S1 - n_Sn) + dt * div(w_S) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h2so4*Dc_h2so4SDzc_h2so4[i] + trf_hso4*Dc_hso4SDzc_h2so4[i] + trf_so4*Dc_so4SDzc_h2so4[i] ;
  }
  K(E_S,I_H2SO4)          += + c[0] ;
  K(E_S,I_H2SO4+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_H2SO4)      += - c[0] ;
  K(E_S+NEQ,I_H2SO4+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hso4*Dc_hso4SDzn_ca_s[i] + trf_so4*Dc_so4SDzn_ca_s[i] ;
  }
  K(E_S,I_Ca_S)           += + c[0] ;
  K(E_S,I_Ca_S+NEQ)       += - c[1] ;
  K(E_S+NEQ,I_Ca_S)       += - c[0] ;
  K(E_S+NEQ,I_Ca_S+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hso4 + tre_so4 ;
  }
  K(E_S,I_psi)          += + c[0] ;
  K(E_S,I_psi+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_psi)      += - c[0] ;
  K(E_S+NEQ,I_psi+NEQ)  += + c[1] ;

  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzc_h2so4[i] ;
  }
  K(E_Ca,I_H2SO4)         += + c[0] ;
  K(E_Ca,I_H2SO4+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_H2SO4)     += - c[0] ;
  K(E_Ca+NEQ,I_H2SO4+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_ca_s[i] ;
  }
  K(E_Ca,I_Ca_S)         += + c[0] ;
  K(E_Ca,I_Ca_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Ca_S)     += - c[0] ;
  K(E_Ca+NEQ,I_Ca_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_ca ;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzc_h2so4[i] + z_oh*trf_oh*Dc_ohSDzc_h2so4[i] + z_hso4*trf_hso4*Dc_hso4SDzc_h2so4[i] + z_so4*trf_so4*Dc_so4SDzc_h2so4[i] \
           + z_ca*trf_ca*Dc_caSDzc_h2so4[i] ;
  }
  K(E_q,I_H2SO4)           += + c[0] ;
  K(E_q,I_H2SO4+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_H2SO4)       += - c[0] ;
  K(E_q+NEQ,I_H2SO4+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_ca_s[i] + z_oh*trf_oh*Dc_ohSDzn_ca_s[i] + z_hso4*trf_hso4*Dc_hso4SDzn_ca_s[i] + z_so4*trf_so4*Dc_so4SDzn_ca_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_ca_s[i] ;
  }
  K(E_q,I_Ca_S)           += + c[0] ;
  K(E_q,I_Ca_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Ca_S)       += - c[0] ;
  K(E_q+NEQ,I_Ca_S+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_psi)          += + c[0] ;
  K(E_q,I_psi+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_psi)      += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)  += + c[1] ;
  
 }
  

#if (U_H2SO4 == LOG_RHO || U_H2SO4 == LOG_ZRHO)
  for(i=0;i<el.nn*NEQ;i++){
    K(i,I_H2SO4)     *= Ln10*ZC_H2SO4(0) ;
    K(i,I_H2SO4+NEQ) *= Ln10*ZC_H2SO4(1) ;
  }
#elif (U_H2SO4 == RHO)
  for(i=0;i<el.nn*NEQ;i++){
    K(i,I_H2SO4)     /= c_h2so4_eq ;
    K(i,I_H2SO4+NEQ) /= c_h2so4_eq ;
  }
#endif

  return(0) ;

#undef K
}


void rs70(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
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
  for(i=0;i<el.nn;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  /*
    Conservation de S (sulfure) : (n_S1 - n_Sn) + dt * div(w_S) = 0
  */
  R(0,E_S) -= volume[0]*(N_S(0) - N_Sn(0)) + dt*surf*W_S ;
  R(1,E_S) -= volume[1]*(N_S(1) - N_Sn(1)) - dt*surf*W_S ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;

#undef R
}


int so70(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  /* if(el.dim < dim) return(0) ; */

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
  
  /*
    Donnees
  */
  phi0     = GET("porosite") ;
  c_h2so4_eq  = GET("C_H2SO4_eq") ;

  /* initialisation */
  nso = 15 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* output quantities */
  {
    /* molarities */
#if (U_H2SO4 == LOG_RHO)
    double zc_h2so4   = exp(Ln10*param(u,h_s,el.nn,I_H2SO4))/c_h2so4_eq ;
#elif (U_H2SO4 == LOG_ZRHO)
    double zc_h2so4    = exp(Ln10*param(u,h_s,el.nn,I_H2SO4)) ;
#elif (U_H2SO4 == ZRHO)
    double zc_h2so4   = param(u,h_s,el.nn,I_H2SO4) ;
#else
    double zc_h2so4   = param(u,h_s,el.nn,I_H2SO4)/c_h2so4_eq ;
#endif
    double zn_ca_s    = param(u,h_s,el.nn,I_Ca_S) ;
    
    int j = (el.dim < dim) ? 0 : ((s[0] < (x[0][0] + x[1][0])*0.5) ? 0 : 1) ;

    /* charge density */
    double c_q = N_q(j) ;
    /* solid contents */
    double n_ch       = N_CH(j) ;
    double n_csh2     = N_CSH2(j) ;
    /* porosity */
    double phi        = PHI(j) ;

    double psi        = param(u,h_s,el.nn,I_psi) ;
      
    CONCENTRATIONS(zc_h2so4,zn_ca_s) ;

    i = 0 ;
    strcpy(r[i].text,"c_h2so4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2so4 ;
    strcpy(r[i].text,"ph") ; r[i].n = 1 ;
    r[i++].v[0] = 14 + log(c_oh)/log(10.) ;
    strcpy(r[i].text,"c_hso4") ; r[i].n = 1 ;
    r[i++].v[0] = c_hso4 ;
    strcpy(r[i].text,"c_so4") ; r[i].n = 1 ;
    r[i++].v[0] = c_so4 ;
    strcpy(r[i].text,"c_ca") ; r[i].n = 1 ;
    r[i++].v[0] = c_ca ;
    strcpy(r[i].text,"n_ch") ; r[i].n = 1 ;
    r[i++].v[0] = n_ch ;
    strcpy(r[i].text,"n_csh2") ; r[i].n = 1 ;
    r[i++].v[0] = n_csh2 ;
    strcpy(r[i].text,"porosite") ; r[i].n = 1 ;
    r[i++].v[0] = phi ;
    strcpy(r[i].text,"c_oh") ; r[i].n = 1 ;
    r[i++].v[0] = c_oh ;
    strcpy(r[i].text,"potentiel_electrique") ; r[i].n = 1 ;
    r[i++].v[0] = psi ;
    strcpy(r[i].text,"charge") ; r[i].n = 1 ;
    r[i++].v[0] = c_q ;
    strcpy(r[i].text,"zn_ca_s") ; r[i].n = 1 ;
    r[i++].v[0] = zn_ca_s ;
    strcpy(r[i].text,"zc_h2so4") ; r[i].n = 1 ;
    r[i++].v[0] = zc_h2so4 ;
    strcpy(r[i].text,"flux de S") ; r[i].n = 1 ;
    r[i++].v[0] = (el.dim < dim) ? 0 : W_S ;
    strcpy(r[i].text,"flux de Ca") ; r[i].n = 1 ;
    r[i++].v[0] = (el.dim < dim) ? 0 : W_Ca ;
  }
  
  if(i != nso) arret("so70") ;

  return(nso) ;
}


void transfert(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom) 
/* Termes explicites (va)  */
{
  int    i ;
  /*
    Donnees
  */
  phi0     = GET("porosite") ;

  /* initialisation */
  for(i=0;i<NVE_TR;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<el.nn;i++) {
    /* molarities */
    double zc_h2so4   = ZC_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
  
    CONCENTRATIONS(zc_h2so4,zn_ca_s) ;
    
    {
    /* solid contents */
    double n_csh2     = N_CSH2(i) ;
    double n_ch       = N_CH(i) ;

    /* porosity */
    double phi        = PHI(i) ;

    /* tortuosite liquide */
    double iff    = TORTUOSITY(phi) ;
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_SO4        += d_so4*iff ;
    KF_H2SO4      += d_h2so4*iff ;
    KF_HSO4       += d_hso4*iff ;

    KF_Ca         += d_ca*iff ;


    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HSO4     += FsRT*KF_HSO4*z_hso4*c_hso4 ;
    Kpsi_SO4      += FsRT*KF_SO4*z_so4*c_so4 ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;

    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hso4*Kpsi_HSO4 + z_so4*Kpsi_SO4 + z_ca*Kpsi_Ca  ;
    }
  }
  
  /* moyenne */
  for(i=0;i<NVE_TR;i++) va[i] *= 0.5 ;
}


void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les flux (f) */
{
  double r_h[2],r_oh[2] ;
  double r_h2so4[2],r_hso4[2],r_so4[2] ;
  double r_ca[2] ;

  int    i ;

  for(i=0;i<el.nn;i++) {
    /* molarities */
    double zc_h2so4   = ZC_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
  
    CONCENTRATIONS(zc_h2so4,zn_ca_s) ;
    
    r_oh[i]       = c_oh ;
    r_h[i]        = c_h ;
    r_h2so4[i]    = c_h2so4 ;
    r_hso4[i]     = c_hso4 ;
    r_so4[i]      = c_so4 ;
    r_ca[i]       = c_ca ;
  }

  /* Gradients */
  {
    double dx = x[1][0] - x[0][0] ;
    double grd_h        = (r_h[1]        - r_h[0]       )/dx ;
    double grd_oh       = (r_oh[1]       - r_oh[0]      )/dx ;
    double grd_h2so4    = (r_h2so4[1]    - r_h2so4[0]   )/dx ;
    double grd_hso4     = (r_hso4[1]     - r_hso4[0]    )/dx ;
    double grd_so4      = (r_so4[1]      - r_so4[0]     )/dx ;
    double grd_ca       = (r_ca[1]       - r_ca[0]      )/dx ;
    
    double grd_psi      = (PSI(1)        - PSI(0)       )/dx ;
    
    /* Flux */
    double w_h2so4    = - KF_H2SO4*grd_h2so4   ;
    double w_hso4     = - KF_HSO4*grd_hso4          - Kpsi_HSO4*grd_psi ;
    double w_so4      = - KF_SO4*grd_so4            - Kpsi_SO4*grd_psi  ;
    double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
    
    double w_q        = - z_h*KF_H*grd_h		        \
                        - z_oh*KF_OH*grd_oh		      \
                        - z_hso4*KF_HSO4*grd_hso4   \
                        - z_so4*KF_SO4*grd_so4		  \
                        - z_ca*KF_Ca*grd_ca		      \
                        - Kpsi_q*grd_psi ;

    W_S     = w_h2so4 + w_hso4 + w_so4 ;
    W_Ca    = w_ca ;
    W_q     = w_q ;
  }
}

void   concentrations(double zc_h2so4,double zn_ca_s,double *var)
{
  double P_CSH2     = IAP_CSH2(zc_h2so4,zn_ca_s) ;
  /* molarities */
  c_h2so4    = zc_h2so4*C_H2SO4_eq ;
  c_oh       = concentration_oh(zc_h2so4,zn_ca_s) ;
  c_h        = K_h2o/c_oh ;
  c_hso4     = K_h2so4*c_h2so4/c_h ;
  c_so4      = K_hso4*c_hso4/c_h ;
  c_ca       = P_CSH2/c_so4 ;
}


double concentration_oh(double zc_h2so4,double zn_ca_s)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* zc_h2so4 et zn_ca_s sont fixes */
  double P_CSH2     = IAP_CSH2(zc_h2so4,zn_ca_s) ;

  /*
  rappel des expressions c_i = A_i*(c_h)**n   : n
     c_oh       = K_h2o/c_h                      : -1
     c_hso4     = K_h2so4*c_h2so4/c_h            : -1
     c_so4      = K_hso4*c_hso4/c_h              : -2
     c_ca       = P_CSH2/c_so4                   : +2   
  */
  double A_h2so4    = zc_h2so4*C_H2SO4_eq ;
  double A_hso4     = K_h2so4*A_h2so4 ;
  double A_so4      = K_hso4*A_hso4 ;
  double A_ca       = P_CSH2/A_so4 ;

  double a = z_ca*A_ca ;
  double b = z_h ;
  double c = 0. ;
  double d = z_oh*K_h2o + z_hso4*A_hso4 ;
  double e = z_so4*A_so4 ;

  /* a > 0 ; b > 0 ; c > 0 ; d < 0 ; e < 0 */
  return(K_h2o/poly4_sansc(a,b,d,e)) ;
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

double poly4_sansc(double a,double b,double d,double e)
/* on resout ax^4 + bx^3 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double y0 = -e/b ;
  double y1 = -(d+e)/(a+b) ;
  double x0 = pow(y0,1./3) ;
  double x1 = pow(y1,1./3) ;
  double x  = x1 ;
  int    i = 0 ;
  
  if(y1 > 1.) {
    printf("a = %e\n",a) ;
    printf("b = %e\n",b) ;
    printf("d = %e\n",d) ;
    printf("e = %e\n",e) ;
    printf("x0 = %e\n",x0) ;
    printf("x1  = %e\n",x1) ;
    printf("y0  = %e\n",y0) ;
    printf("y1  = %e\n",y1) ;
    return(-1.) ;
    arret("poly4_sansc : valeur non realiste") ;
  }

  do {
    double x2 = x*x,x3 = x2*x ;
    double f  = x3 + (d*x + e)/(a*x + b) ;
    double df = 3*x2 + (d*b - e*a)/((a*x + b)*(a*x + b)) ;
    double dx = -f/df ;
    err = fabs(dx/x) ;
    x += dx ;
    if(i++ > 20) {
      printf("a = %e\n",a) ;
      printf("b = %e\n",b) ;
      printf("d = %e\n",d) ;
      printf("e = %e\n",e) ;
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      return(-1.) ;
      arret("poly4_sansc : non convergence") ;
    }
  } while(err > tol) ;
  return(x) ;
}
