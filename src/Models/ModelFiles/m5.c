#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  5
#define TITLE   "Sechage isotherme"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ      (2)

#define NVI      (9)
#define NVE      (5)

#define E_eau    (0)
#define E_air    (1)

#define I_p_l    (0)
#define I_p_a    (1)

#define P_l(n)   (u[(n)][I_p_l])
#define P_a(n)   (u[(n)][I_p_a])

#define M_l(n)   (f[(n)])
#define M_v(n)   (f[(n+2)])
#define M_a(n)   (f[(n+4)])
#define W_l      (f[(6)])
#define W_v      (f[(7)])
#define W_a      (f[(8)])

#define M_ln(n)  (f_n[(n)])
#define M_vn(n)  (f_n[(n+2)])
#define M_an(n)  (f_n[(n+4)])

#define K_l      (va[(0)])
#define KD_v     (va[(1)])
#define KF_v     (va[(2)])
#define KD_a     (va[(3)])
#define KF_a     (va[(4)])

/* Masses molaires (kg) */
#define M_h2o    (18.e-3)
#define M_air    (28.8e-3)

/* Coefficient de diffusion moleculaire (m2/s) */
#define D_av0    (2.48e-5)

/* Constantes physiques */
#define RT       (2436.)     /* produit de R et T (J/mol) */
                             /* R = 8.3143 J/K/mol cste des gaz parfait */
                             /* T = 293 K */

/* Viscosites (Pa.s) */
#define mu_l     (1.e-3)
#define mu_g     (1.8e-5)

/* masse volumique (kg/m3) */
#define rho_l    (1000.)


/* Fonctions */
static int    pm(const char *s) ;
static double saturation(double,double,crbe_t) ;
static double dsaturation(double,double,crbe_t) ;
/* Parametres */
static double gravite,phi,k_int,p_l0,p_v0,p_a0,p_c3 ;

/* Macros pour les fonctions saturation */
#define SATURATION(x)             saturation(x,p_c3,el.mat->cb[0])
#define DSATURATION(x)            dsaturation(x,p_c3,el.mat->cb[0])
#define PERMEABILITYTOLIQUID(x)   courbe(x,el.mat->cb[1])
#define PERMEABILITYTOGAS(x)      courbe(x,el.mat->cb[2])
#define TORTUOSITYTOGAS(f,sg)     (pow(f,aa)*pow(sg,bb))
#define aa                        (1./3)  /* 1/3 Millington, 2.74 */
#define bb                        (7./3)  /* 7/3 Millington,  4.2 */



int pm(const char *s)
{
  if(!strcmp(s,"gravite")) return (0) ;
  else if(!strcmp(s,"phi")) return (1) ;
  else if(!strcmp(s,"k_int")) return (2) ;
  else if(!strcmp(s,"p_l0")) return (3) ;
  else if(!strcmp(s,"p_v0")) return (4) ;
  else if(!strcmp(s,"p_a0")) return (5) ;
  else if(!strcmp(s,"p_c3")) return (6) ;
  else if(!strcmp(s,"courbes")) return (7) ;
  else
    { printf("donnee \"%s\" non connue (pm5)\n",s) ; exit(0) ; }
}

int dm5(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 7,n_courbes = 3 ;

  mat->neq    = NEQ ;
  strcpy(mat->eqn[E_eau],"eau") ;
  strcpy(mat->eqn[E_air],"air") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_p_a],"p_a") ;

  lit_mate(mat,ficd,pm,n_donnees,n_courbes) ;
  return(mat->n) ;
}

int qm5(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est formee :\n\
\t 1. Conservation de la masse d\'eau (p_l = pression de liquide)\n\
\t 2. Conservation de la masse d\'air (p_a = pression d\'air)\n\
") ;

  printf("\n\
Exemple de donnees (en unites MKS)\n\n") ;

  fprintf(ficd,"gravite = 0     # Gravite\n") ;
  fprintf(ficd,"phi = 0.3       # Porosite\n") ;
  fprintf(ficd,"k_int = 1.e-20  # Permeabilite intrinseque\n") ;
  fprintf(ficd,"p_l0 = 100000   # Pression de liquide de reference\n") ;
  fprintf(ficd,"p_v0 = 2460     # Pression de vapeur de reference\n") ;
  fprintf(ficd,"p_a0 = 97540    # Pression d\'air de reference\n") ;
  fprintf(ficd,"p_c3 = 1.e6     # Pression capillaire limite\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}

void tb5(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}

void ch5(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int i ;
  
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<el.nn*NEQ;i++) r[i] = -r[i] ;
}

void in5(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double p_l,p_a,p_c,p_v,p_g,pg[2],s_l,s_g ;
  double rho_v,rho_a,rho_g ;
  double cv[2],ca[2] ;
  double dx ;
  int    i ;
  double un = 1. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  
  /* masses d'eau et d'air */
  for(i = 0 ; i < 2 ; i++) {
    p_l     = P_l(i) ;
    p_v     = p_v0*exp(M_h2o/RT*(p_l - p_l0)/rho_l) ;
    p_a     = (P_a(i) > 0) ? P_a(i) : (-P_a(i)) - p_v ;
    p_g     = p_v + p_a ;
    p_c     = p_g - p_l ;

    pg[i]   = p_g ;

    s_l     = SATURATION(p_c) ;
    s_g     = un - s_l ;

    rho_v   = p_v*M_h2o/RT ;
    rho_a   = p_a*M_air/RT ;
    rho_g   = rho_v + rho_a ;

    cv[i]   = rho_v/rho_g ;
    ca[i]   = un - cv[i] ;

    M_l(i)  = rho_l*phi*s_l ;
    M_v(i)  = rho_v*phi*s_g ;
    M_a(i)  = rho_a*phi*s_g ;
    
    P_a(i)  = p_a ;
  }

  /* coefficient de transfert */
  {
    /* ex_t ex5 ; */
    ex5(x,u,f,va,el,dim,geom,0.) ;
  }

  /* flux */
  dx    = x[1][0] - x[0][0] ;
  W_l   = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
  W_v   = - KD_v*(pg[1] - pg[0])/dx - KF_v*(cv[1] - cv[0])/dx ;
  W_a   = - KD_a*(pg[1] - pg[0])/dx - KF_a*(ca[1] - ca[0])/dx ;
}

int ex5(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double p_l,p_a,p_c,s_l,s_g,p_v,p_g,c_v,c_a ;
  double rho_v,rho_a,rho_g ;
  double D_av,D_eff,tau,kh_l,kh_g ;
  double un = 1.,deux = 2. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /*
    COEFFICIENTS DE TRANSFERT
  */
  p_l    = (P_l(0) + P_l(1))/deux ;
  p_a    = (P_a(0) + P_a(1))/deux ;
  p_v    = p_v0*exp(M_h2o/RT*(p_l - p_l0)/rho_l) ;
  p_g    = p_v + p_a ;
  p_c    = p_g - p_l ;

  s_l    = SATURATION(p_c) ;
  s_g    = un - s_l ;

  rho_v  = M_h2o*p_v/RT ;
  rho_a  = M_air*p_a/RT ;
  rho_g  = rho_v + rho_a ;

  c_v    = rho_v/rho_g ;
  c_a    = un - c_v ;
  
  kh_l  = k_int/mu_l*PERMEABILITYTOLIQUID(p_c) ;  /* permeabilite liquide */
  K_l   = rho_l*kh_l ;                            /* Darcy liquide */

  kh_g  = k_int/mu_g*PERMEABILITYTOGAS(p_c) ;  /* permeabilite gaz */
  KD_v  = rho_v*kh_g ;                            /* Darcy vapeur */
  KD_a  = rho_a*kh_g ;                            /* Darcy air */

  tau   = TORTUOSITYTOGAS(phi,s_g) ;           /* tortuosite */
  D_av  = D_av0*(p_v0+p_a0)/p_g ;                 /* diffusion moleculaire */
  D_eff = phi*s_g*tau*D_av ;                      /* diffusion effective */
  KF_v  = rho_g*D_eff ;                           /* Fick vapeur */
  KF_a  = KF_v ;                                  /* Fick air */
  KD_v += D_eff*c_v*c_a*(M_air - M_h2o)/RT ; /* diffusion barometrique vap */
  KD_a += D_eff*c_v*c_a*(M_h2o - M_air)/RT ; /* diffusion barometrique air */
  return(0) ;
}

int ct5(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double rho_v,rho_a,rho_g ;
  double p_l,p_a,p_c,s_l,s_g,p_v,p_g,pg[2] ;
  double cv[2],ca[2] ;
  double dx ;
  int    i ;
  double zero = 0.,un = 1. ;
  
  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /* masses d'eau et d'air */
  for(i=0;i<2;i++) {
    if(P_a(i) <= zero)  {
      printf("ct5 : pression d\'air = %e\n",P_a(i)) ;
      return(1) ;
    }
    p_l     = P_l(i) ;
    p_a     = P_a(i) ;
    p_v     = p_v0*exp(M_h2o/RT*(p_l - p_l0)/rho_l) ;
    p_g     = p_v + p_a ;
    p_c     = p_g - p_l ;

    pg[i]   = p_g ;

    s_l     = SATURATION(p_c) ;
    s_g     = un - s_l ;

    rho_v   = p_v*M_h2o/RT ;
    rho_a   = p_a*M_air/RT ;
    rho_g   = rho_v + rho_a ;

    cv[i]   = rho_v/rho_g ;
    ca[i]   = un - cv[i] ;

    M_l(i) = rho_l*phi*s_l ;
    M_v(i) = rho_v*phi*s_g ;
    M_a(i) = rho_a*phi*s_g ;
  }
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  W_l  = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
  W_v  = - KD_v*(pg[1] - pg[0])/dx - KF_v*(cv[1] - cv[0])/dx ;
  W_a  = - KD_a*(pg[1] - pg[0])/dx - KF_a*(ca[1] - ca[0])/dx ;
  return(0) ;
}

int mx5(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
  /* Numeros d'ordre des equations et des inconnues locales */
#define NEQ1     (NEQ+1)
#define E_vap    NEQ
#define I_p_v    NEQ
#define E_liq    E_eau

#define KE(i,j)  (ke[(i)*2*NEQ1+(j)])
#define K(i,j)   (k[(i)*2*NEQ+(j)])
  double p_c,s_l,s_g,p_v,p_l,p_a,p_g,dslsdpc,pv[2],pa[2],pg[2],dpvsdpl[2] ;
  double rho_v,rho_a,rho_g ;
  double cv[2],ca[2] ;
  double tr_l,trd_v,trf_v,trd_a,trf_a ;
  double dx ,xm ;
  double volume[2],surf ;
  int    i,j,n ;
  double zero = 0.,un = 1.,deux = 2. ;
  double ke[4*NEQ1*NEQ1] ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /*
    INITIALISATION DE LA MATRICE
  */
  for(i=0;i<4*NEQ1*NEQ1;i++) ke[i] = zero ;
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
  
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    p_l   = P_l(i) ;
    p_a   = P_a(i) ;
    p_v   = p_v0*exp(M_h2o/RT*(p_l - p_l0)/rho_l) ;
    p_g   = p_v + p_a ;
    p_c   = p_g - p_l ;

    pv[i] = p_v ;
    pa[i] = p_a ;
    pg[i] = p_g ;

    s_l   = SATURATION(p_c) ;
    s_g   = un - s_l ;

    rho_v  = p_v*M_h2o/RT ;
    rho_a  = p_a*M_air/RT ;
    rho_g  = rho_v + rho_a ;

    cv[i]  = rho_v/rho_g ;
    ca[i]  = un - cv[i] ;

    dslsdpc = DSATURATION(p_c) ;
    dpvsdpl[i] = rho_v/rho_l ;
    /*
      CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l) = +A
    */
    KE(E_liq+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*phi*rho_l*(-dslsdpc) ;
    KE(E_liq+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*rho_l*dslsdpc ;
    KE(E_liq+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*phi*rho_l*dslsdpc ;
    /*
      CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v) = -A
    */
    KE(E_vap+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*phi*rho_v*dslsdpc ;
    KE(E_vap+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*(M_h2o/RT*s_g - rho_v*dslsdpc) ;
    KE(E_vap+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*phi*rho_v*(-dslsdpc) ;
    /*
      CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a) = 0
    */
    KE(E_air+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*phi*rho_a*dslsdpc ;
    KE(E_air+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*rho_a*(-dslsdpc) ;
    KE(E_air+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*phi*(M_air/RT*s_g - rho_a*dslsdpc) ;
  }
  /*
    termes d'ecoulement
  */
  tr_l  = dt*surf/dx*K_l ;
  trd_v = dt*surf/dx*KD_v ;
  trd_a = dt*surf/dx*KD_a ;
  trf_v = dt*surf/dx*KF_v ;
  trf_a = dt*surf/dx*KF_a ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l) = +A
  */
  KE(E_liq,I_p_l)           += + tr_l ;
  KE(E_liq,I_p_l+NEQ1)      += - tr_l ;
  KE(E_liq+NEQ1,I_p_l)      += - tr_l ;
  KE(E_liq+NEQ1,I_p_l+NEQ1) += + tr_l ;
  /*
    CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v) = -A
  */
  KE(E_vap,I_p_v)      += + trd_v + trf_v*ca[0]*cv[0]/pv[0] ;
  KE(E_vap,I_p_v+NEQ1) += - trd_v - trf_v*ca[1]*cv[1]/pv[1] ;
  KE(E_vap,I_p_a)      += + trd_v + trf_v*(-ca[0]*cv[0]/pa[0]) ;
  KE(E_vap,I_p_a+NEQ1) += - trd_v - trf_v*(-ca[1]*cv[1]/pa[1]) ;
  
  KE(E_vap+NEQ1,I_p_v+NEQ1) += + trd_v + trf_v*ca[1]*cv[1]/pv[1] ;
  KE(E_vap+NEQ1,I_p_v)      += - trd_v - trf_v*ca[0]*cv[0]/pv[0] ;
  KE(E_vap+NEQ1,I_p_a+NEQ1) += + trd_v + trf_v*(-ca[1]*cv[1]/pa[1]) ;
  KE(E_vap+NEQ1,I_p_a)      += - trd_v - trf_v*(-ca[0]*cv[0]/pa[0]) ;
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a) = 0
  */
  KE(E_air,I_p_a)      += + trd_a + trf_a*ca[0]*cv[0]/pa[0] ;
  KE(E_air,I_p_a+NEQ1) += - trd_a - trf_a*ca[1]*cv[1]/pa[1] ;
  KE(E_air,I_p_v)      += + trd_a + trf_a*(-ca[0]*cv[0]/pv[0]) ;
  KE(E_air,I_p_v+NEQ1) += - trd_a - trf_a*(-ca[1]*cv[1]/pv[1]) ;
  
  KE(E_air+NEQ1,I_p_a+NEQ1) += + trd_a + trf_a*ca[1]*cv[1]/pa[1] ;
  KE(E_air+NEQ1,I_p_a)      += - trd_a - trf_a*ca[0]*cv[0]/pa[0] ;
  KE(E_air+NEQ1,I_p_v+NEQ1) += + trd_a + trf_a*(-ca[1]*cv[1]/pv[1]) ;
  KE(E_air+NEQ1,I_p_v)      += - trd_a - trf_a*(-ca[0]*cv[0]/pv[0]) ;

  /*
    ELIMINATION DE LA COLONNE I_p_v :
    COLONNE I_p_l += COLONNE I_p_v * dpvsdpl
  */
  for(j=0;j<NEQ1;j++) {
    for(i=0;i<2;i++) {
      KE(j+i*NEQ1,I_p_l+i*NEQ1) += KE(j+i*NEQ1,I_p_v+i*NEQ1) * dpvsdpl[i] ;
    }
    KE(j,I_p_l+NEQ1) += KE(j,I_p_v+NEQ1) * dpvsdpl[1] ;
    KE(j+NEQ1,I_p_l) += KE(j+NEQ1,I_p_v) * dpvsdpl[0] ;
  }
  /*
    ELIMINATION DE LA LIGNE E_vap :
    LIGNE E_liq += LIGNE E_vap
  */
  for(j=0;j<NEQ;j++) {
    for(i=0;i<2;i++) KE(E_liq+i*NEQ1,j+i*NEQ1) += KE(E_vap+i*NEQ1,j+i*NEQ1) ;
    KE(E_liq,j+NEQ1) += KE(E_vap,j+NEQ1) ;
    KE(E_liq+NEQ1,j) += KE(E_vap+NEQ1,j) ;
  }

  /*
    ASSEMBLAGE DANS K
  */
  for(i=0;i<NEQ;i++) {
    for(j=0;j<NEQ;j++) {
      for(n=0;n<2;n++)  K(i+n*NEQ,j+n*NEQ) += KE(i+n*NEQ1,j+n*NEQ1) ;
      K(i,j+NEQ) += KE(i,j+NEQ1) ;
      K(i+NEQ,j) += KE(i+NEQ1,j) ;
    }
  }

  return(0) ;

#undef NEQ1
#undef E_vap
#undef I_p_v
#undef E_liq
#undef K
#undef KE
}

void rs5(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  /* Numeros d'ordre des equations et des inconnues locales */
#define NEQ1     (NEQ+1)
#define E_vap    NEQ
#define E_liq    E_eau
#define R(n,i)    (r[(n)*NEQ+(i)])
#define RE(n,i)   (re[(n)*NEQ1+(i)])
  double dx ,xm ;
  double volume[2],surf ;
  int    i,j ;
  double zero = 0.,un = 1.,deux = 2. ;
  double re[NEQ1*2] ;
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;

  /*
    INITIALISATION DU RESIDU
  */
  for(i=0;i<NEQ1*2;i++) re[i] = zero ;
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
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l) = +A
  */
  RE(0,E_liq) -= volume[0]*(M_l(0) - M_ln(0)) + dt*surf*W_l ;
  RE(1,E_liq) -= volume[1]*(M_l(1) - M_ln(1)) - dt*surf*W_l ;
  /*
    CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v) = -A
  */
  RE(0,E_vap) -= volume[0]*(M_v(0) - M_vn(0)) + dt*surf*W_v ;
  RE(1,E_vap) -= volume[1]*(M_v(1) - M_vn(1)) - dt*surf*W_v ;
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a) = 0
  */
  RE(0,E_air) -= volume[0]*(M_a(0) - M_an(0)) + dt*surf*W_a ;
  RE(1,E_air) -= volume[1]*(M_a(1) - M_an(1)) - dt*surf*W_a ;
  /*
    ELIMINATION DE LA LIGNE E_vap : LIGNE E_liq += LIGNE E_vap
  */
  for(i=0;i<2;i++) RE(i,E_liq) += RE(i,E_vap) ;
  /*
    ASSEMBLAGE DANS R
  */
  for(j=0;j<NEQ;j++) for(i=0;i<2;i++) R(i,j) += RE(i,j) ;

#undef NEQ1
#undef E_vap
#undef E_liq
#undef R
#undef RE
}

int so5(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double p_l,p_v,p_a,p_g,p_c,s_l,w_g ;
  int    i,j,nso ;                                                           
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /* initialisation */
  nso = 8 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pressions */
  p_l =  param(u,h_s,el.nn,I_p_l) ;
  p_a =  param(u,h_s,el.nn,I_p_a) ;
  /* autres pressions */
  p_v = p_v0*exp(M_h2o/RT*(p_l - p_l0)/rho_l) ;
  p_g = p_v + p_a ;
  p_c = p_g - p_l ;
  /* saturation */
  s_l = SATURATION(p_c) ;
  /* flux */
  w_g  = W_v + W_a ;

  /* quantites exploitees */
  strcpy(r[0].text,"p_l") ; r[0].n = 1 ;
  r[0].v[0] = p_l ;
  strcpy(r[1].text,"p_a") ; r[1].n = 1 ;
  r[1].v[0] = p_a ;
  strcpy(r[2].text,"p_v") ; r[2].n = 1 ;
  r[2].v[0] = p_v ;
  strcpy(r[3].text,"w_l") ; r[3].n = 1 ;
  r[3].v[0] = W_l ;
  strcpy(r[4].text,"w_a") ; r[4].n = 1 ;
  r[4].v[0] = W_a ;
  strcpy(r[5].text,"w_v") ; r[5].n = 1 ;
  r[5].v[0] = W_v ;
  strcpy(r[6].text,"w_g") ; r[6].n = 1 ;
  r[6].v[0] = w_g ;
  strcpy(r[7].text,"saturation") ; r[7].n = 1 ;
  r[7].v[0] = s_l ;
  return (nso) ;
}

double saturation(double pc,double pc3,crbe_t cb)
/* Degre de saturation regularise autour de 1 */
{
  double sl,sl3 ;
  double pc1 = cb.a[0],sl1 = cb.f[0] ;
  
  if(pc >= pc3 || pc1 >= pc3) sl = courbe(pc,cb) ;
  else {
    sl3 = courbe(pc3,cb) ;
    sl  = sl1 - (sl1 - sl3) * exp(-(pc - pc3)/(pc1 - pc3)) ;
  }
  return(sl) ;
}

double dsaturation(double pc,double pc3,crbe_t cb)
{
  int    n_i = cb.n - 1 ;
  double dpc = (cb.a[1] - cb.a[0])/n_i ;

  return((saturation(pc + dpc,pc3,cb) - saturation(pc,pc3,cb))/dpc) ;
}
