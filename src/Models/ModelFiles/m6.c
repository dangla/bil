#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  6
#define TITLE   "Sechage anisotherme"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ      (3)

#define NVI      (18)
#define NVE      (6)

#define E_eau    (0)
#define E_air    (1)
#define E_the    (2)

#define I_p_l    (0)
#define I_p_a    (1)
#define I_tem    (2)

#define P_l(n)   (u[(n)][I_p_l])
#define P_a(n)   (u[(n)][I_p_a])
#define TEM(n)   (u[(n)][I_tem])

#define M_l(n)   (f[(n)])
#define M_v(n)   (f[(n+2)])
#define M_a(n)   (f[(n+4)])
#define S(n)     (f[(n+6)])
#define S_l(n)   (f[(n+8)])
#define S_v(n)   (f[(n+10)])
#define S_a(n)   (f[(n+12)])
#define W_l      (f[(14)])
#define W_v      (f[(15)])
#define W_a      (f[(16)])
#define Q        (f[(17)])

#define P_ln(n)  (u_n[(n)][I_p_l])
#define P_an(n)  (u_n[(n)][I_p_a])
#define TEM_n(n) (u_n[(n)][I_tem])

#define M_ln(n)  (f_n[(n)])
#define M_vn(n)  (f_n[(n+2)])
#define M_an(n)  (f_n[(n+4)])
#define S_n(n)   (f_n[(n+6)])
#define S_ln(n)  (f_n[(n+8)])
#define S_vn(n)  (f_n[(n+10)])
#define S_an(n)  (f_n[(n+12)])
#define W_ln     (f_n[(14)])
#define W_vn     (f_n[(15)])
#define W_an     (f_n[(16)])
#define Q_n      (f_n[(17)])

#define K_l      (va[(0)])
#define KD_v     (va[(1)])
#define KF_v     (va[(2)])
#define KD_a     (va[(3)])
#define KF_a     (va[(4)])
#define KTH      (va[(5)])

/* Fonctions */
static int    pm(const char *s) ;
static double saturation(double,double,crbe_t) ;
static double dsaturation(double,double,crbe_t) ;
static double dSsdT(double,double,double,double,crbe_t) ;
/* Parametres */
static double rho_l,mu_l,mu_g,M_vsR,M_asR,D_av0,T_0 ; /* Pourrait etre en dur */
static double lam_l,lam_g,C_pl,C_pv,C_pa,L_0 ;        /* Ca aussi ! */
static double lam_s,C_s ;
static double gravite,phi,k_int,p_l0,p_v0,p_a0,alpha,p_c3 ;

/* Macros pour les fonctions saturation */
#define SATURATION(x)   saturation(x,p_c3,el.mat->cb[0])
#define DSATURATION(x)  dsaturation(x,p_c3,el.mat->cb[0])

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"mu_g") == 0) return (5) ;
  else if(strcmp(s,"p_l0") == 0) return (6) ;
  else if(strcmp(s,"p_v0") == 0) return (7) ;
  else if(strcmp(s,"p_a0") == 0) return (8) ;
  else if(strcmp(s,"T_0") == 0) return (9) ;
  else if(strcmp(s,"M_vsR") == 0) return (10) ;
  else if(strcmp(s,"M_asR") == 0) return (11) ;
  else if(strcmp(s,"D_av") == 0) return (12) ;
  else if(strcmp(s,"lam_s") == 0) return (13) ;
  else if(strcmp(s,"lam_l") == 0) return (14) ;
  else if(strcmp(s,"lam_g") == 0) return (15) ;
  else if(strcmp(s,"C_s") == 0) return (16) ;
  else if(strcmp(s,"C_pl") == 0) return (17) ;
  else if(strcmp(s,"C_pv") == 0) return (18) ;
  else if(strcmp(s,"C_pa") == 0) return (19) ;
  else if(strcmp(s,"L_0") == 0) return (20) ;
  else if(strcmp(s,"a") == 0) return (21) ;
  else if(strcmp(s,"p_c3") == 0) return (22) ;
  else if(strcmp(s,"courbes") == 0) return (23) ;
  else
    { printf("donnee \"%s\" non connue (pm6)\n",s) ; exit(0) ; }
}

int dm6(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 24 ;

  mat->neq    = NEQ ;

  strcpy(mat->eqn[E_eau],"eau") ;
  strcpy(mat->eqn[E_air],"air") ;
  strcpy(mat->eqn[E_the],"the") ;

  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_p_a],"p_a") ;
  strcpy(mat->inc[I_tem],"tem") ;
  
  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm6(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est formee :\n\
\t 1. l\'equation de la chaleur       (tem)\n\
\t 2. Conservation de la masse d\'eau (p_l)\n\
\t 3. Conservation de la masse d\'air (p_a)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 0     # Gravite\n") ;
  fprintf(ficd,"phi = 0.3       # Porosite\n") ;
  fprintf(ficd,"rho_l = 1000    # Masse volumique du liquide\n") ;
  fprintf(ficd,"M_asR = 0.00346 # Masse molaire de l\'air sur R\n") ;
  fprintf(ficd,"M_vsR = 0.00216 # Masse molaire de la vapeur sur R\n") ;
  fprintf(ficd,"k_int = 1.e-20  # Permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001    # Viscosite du liquide\n") ;
  fprintf(ficd,"mu_g = 1.8e-05  # Viscosite du gaz\n") ;
  fprintf(ficd,"p_l0 = 100000   # Pression de liquide de reference\n") ;
  fprintf(ficd,"p_v0 = 2460     # Pression de vapeur de reference\n") ;
  fprintf(ficd,"p_a0 = 97540    # Pression d\'air de reference\n") ;
  fprintf(ficd,"T_0 = 293       # Temperature\n") ;
  fprintf(ficd,"D_av = 0.248e-4 # Coefficient de diffusion air-vapeur\n") ;
  fprintf(ficd,"lam_s = 1.12    # Conductivite thermique du solide\n") ;
  fprintf(ficd,"lam_l = 0.6     # Conductivite thermique du liquide\n") ;
  fprintf(ficd,"lam_g = 0.026   # Conductivite thermique du gaz\n") ;
  fprintf(ficd,"C_s = 2.3e+06   # Chaleur volumique du solide\n") ;
  fprintf(ficd,"C_pl = 4180     # Chaleur specifique du liquide\n") ;
  fprintf(ficd,"C_pv = 1800     # Chaleur specifique de la vapeur\n") ;
  fprintf(ficd,"C_pa = 1000     # Chaleur specifique de l\'air\n") ;
  fprintf(ficd,"L_0 = 2.45e+06  # Chaleur latente de changement de phase\n") ;
  fprintf(ficd,"a = 0.003       # S_l(p_c,T) = S_l(p_c/(1-a(T-T_0))\n") ;
  fprintf(ficd,"p_c3 = 1.e6     # Pression capillaire limite\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}

void tb6(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}

void ch6(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int    i ;
  
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;

  for(i=0;i<el.nn*NEQ;i++) r[i] = -r[i] ;

  if(isdigit(cg.eqn[0]) && (atoi(cg.eqn) - 1) == E_the) {
    R(0,E_the) /= TEM_n(0) ;
  } else if(!strcmp(cg.eqn,el.mat->eqn[E_the])) {
    R(0,E_the) /= TEM_n(0) ;
  }
#undef R
}

void in6(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double c_v[2],c_a[2] ;
  double p_g[2] ;
  double dx ;
  int    i,j ;
  int    np = 3 ;
  double a[3] = {0.93246951420,0.66120938646,0.23861918608} ;
  double w[3] = {0.17132449237,0.36076157304,0.46791393457} ;
  double zero = 0.,un = 1.,deux = 2. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_vsR   = el.mat->pr[pm("M_vsR")] ;
  M_asR   = el.mat->pr[pm("M_asR")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_av0   = el.mat->pr[pm("D_av")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  lam_l   = el.mat->pr[pm("lam_l")] ;
  lam_g   = el.mat->pr[pm("lam_g")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  C_pl    = el.mat->pr[pm("C_pl")] ;
  C_pv    = el.mat->pr[pm("C_pv")] ;
  C_pa    = el.mat->pr[pm("C_pa")] ;
  L_0     = el.mat->pr[pm("L_0")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  alpha   = el.mat->pr[pm("a")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  
  dx = x[1][0] - x[0][0] ;
  /* masses d'eau et d'air, entropies */
  for(i=0;i<2;i++) {
    double theta  = TEM(i) - T_0 ;
    double pl = P_l(i) ;
    double pa = P_a(i) ;
    double pv = p_v0*exp(M_vsR/TEM(i)*((pl - p_l0)/rho_l + L_0*theta/T_0
	   + (C_pl - C_pv)*(theta - TEM(i)*log(TEM(i)/T_0)))) ;
    double pg = pv + pa ;
    double pc = pg - pl ;
    double sl = SATURATION(pc/(un-alpha*theta)) ;
    double sg = un - sl ;
    double rho_v  = pv*M_vsR/TEM(i) ;
    double rho_a  = pa*M_asR/TEM(i) ;
    double rho_g  = rho_v + rho_a ;
    double cv = rho_v/rho_g ;
    double ca = un - cv ;
    double dUsdT  = zero ;
    for(j=0;j<np;j++) {
      dUsdT += w[j]*(
	       dSsdT(pc/deux*(un+a[j]),p_c3,theta,alpha,el.mat->cb[0])
	     + dSsdT(pc/deux*(un-a[j]),p_c3,theta,alpha,el.mat->cb[0])) ;
    }
    dUsdT *= pc/deux ;
    M_l(i) = rho_l*phi*sl ;
    M_v(i) = rho_v*phi*sg ;
    M_a(i) = rho_a*phi*sg ;
    S_l(i) = C_pl*log(TEM(i)/T_0) ;
    S_v(i) = C_pv*log(TEM(i)/T_0) - log(pv/p_v0)/M_vsR + L_0/T_0 ;
    S_a(i) = C_pa*log(TEM(i)/T_0) - log(P_a(i)/p_a0)/M_asR  ;
    S(i)   = C_s*log(TEM(i)/T_0) - phi*dUsdT
           + M_l(i)*S_l(i) + M_v(i)*S_v(i) + M_a(i)*S_a(i) ;
           
    
    p_g[i] = pg ;
    c_v[i] = cv ;
    c_a[i] = ca ;
  }
  /* coefficient de transfert */
  {
    ex_t ex6 ;
    ex6(x,u,f,va,el,dim,geom,0.) ;
  }
  /* flux */
  {
    double Tij   = (TEM(0) + TEM(1))*0.5 ;
    double S_lij = (S_l(0) + S_l(1))*0.5 ;
    double S_vij = (S_v(0) + S_v(1))*0.5 ;
    double S_aij = (S_a(0) + S_a(1))*0.5 ;
    W_l   = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
    W_v   = - KD_v*(p_g[1] - p_g[0])/dx - KF_v*(c_v[1] - c_v[0])/dx ;
    W_a   = - KD_a*(p_g[1] - p_g[0])/dx - KF_a*(c_a[1] - c_a[0])/dx ;
    Q     = - KTH/Tij*(TEM(1) - TEM(0))/dx 
          +   S_lij*W_l + S_vij*W_v + S_aij*W_a ;
  }
}

int ex6(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double p_lij,p_aij,Tij ;
  double pc,sl,sg,pv,pg,theta,cv,ca,pc0 ;
  double rho_v,rho_a,rho_g ;
  double D_av,D_eff,tau,kh_l,kh_g ;
  double un = 1.,deux = 2. ;

  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_vsR   = el.mat->pr[pm("M_vsR")] ;
  M_asR   = el.mat->pr[pm("M_asR")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_av0   = el.mat->pr[pm("D_av")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  lam_l   = el.mat->pr[pm("lam_l")] ;
  lam_g   = el.mat->pr[pm("lam_g")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  C_pl    = el.mat->pr[pm("C_pl")] ;
  C_pv    = el.mat->pr[pm("C_pv")] ;
  C_pa    = el.mat->pr[pm("C_pa")] ;
  L_0     = el.mat->pr[pm("L_0")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  alpha   = el.mat->pr[pm("a")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /*
    COEFFICIENTS DE TRANSFERT
  */
  p_lij = (P_l(0) + P_l(1))/deux ;
  p_aij = (P_a(0) + P_a(1))/deux ; 
  Tij   = (TEM(0) + TEM(1))/deux ;  
  
  theta = Tij - T_0 ;
  pv    = p_v0*exp(M_vsR/Tij*((p_lij - p_l0)/rho_l + L_0*theta/T_0
	+ (C_pl - C_pv)*(theta - Tij*log(Tij/T_0)))) ;
  pg    = pv + p_aij ;
  pc    = pg - p_lij ;
  pc0   = pc/(un-alpha*theta) ;
  sl    = SATURATION(pc0) ;
  sg    = un - sl ;
  rho_v = M_vsR*pv/Tij ;
  rho_a = M_asR*p_aij/Tij ;
  rho_g = rho_v + rho_a ;
  cv    = rho_v/rho_g ;
  ca    = un - cv ;

  kh_l  = k_int/mu_l*courbe(pc0,el.mat->cb[1]) ;    /* permeabilite liquide */
  K_l   = rho_l*kh_l ;                              /* Darcy liquide */

  kh_g  = k_int/mu_g*courbe(pc0,el.mat->cb[2]) ;    /* permeabilite gaz */
  KD_v  = rho_v*kh_g ;                              /* Darcy vapeur */
  KD_a  = rho_a*kh_g ;                              /* Darcy air */

  tau   = pow(phi,1./3.)*pow(sg,7./3.) ;            /* tortuosite */
  D_av  = D_av0*(p_v0+p_a0)/pg*pow(Tij/T_0,1.88) ;  /* diffusion moleculaire */
  D_eff = phi*sg*tau*D_av ;                         /* diffusion effective */
  KF_v  = rho_g*D_eff ;                             /* Fick vapeur */
  KF_a  = KF_v ;                                    /* Fick air */
  KD_v += D_eff*cv*ca*(M_asR - M_vsR)/Tij ; /* diffusion barometrique vap */
  KD_a += D_eff*cv*ca*(M_vsR - M_asR)/Tij ; /* diffusion barometrique air */

  KTH   = pow(lam_s,un-phi)*pow(lam_l,phi*sl)*pow(lam_g,phi*sg) ; /* Fourier */
  return(0) ;
}

int ct6(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double c_v[2],c_a[2] ;
  double Tij,S_lij,S_vij,S_aij ;
  double pc,pc0,sl,sg,pv,dUsdT,theta,p_g[2] ;
  double rho_v,rho_a,rho_g ;
  double dx ;
  int    i,j ;
  int    np = 3 ;
  double a[3] = {0.93246951420,0.66120938646,0.23861918608} ;
  double w[3] = {0.17132449237,0.36076157304,0.46791393457} ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_vsR   = el.mat->pr[pm("M_vsR")] ;
  M_asR   = el.mat->pr[pm("M_asR")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_av0   = el.mat->pr[pm("D_av")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  lam_l   = el.mat->pr[pm("lam_l")] ;
  lam_g   = el.mat->pr[pm("lam_g")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  C_pl    = el.mat->pr[pm("C_pl")] ;
  C_pv    = el.mat->pr[pm("C_pv")] ;
  C_pa    = el.mat->pr[pm("C_pa")] ;
  L_0     = el.mat->pr[pm("L_0")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  alpha   = el.mat->pr[pm("a")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /* masses d'eau et d'air, entropie */
  for(i=0;i<2;i++) {
    if(P_a(i) <= zero)  {
      printf("pression d\'air = %e\n",P_a(i)) ;
      return(1) ;
    }
    theta  = TEM(i) - T_0 ;
    pv     = p_v0*exp(M_vsR/TEM(i)*((P_l(i) - p_l0)/rho_l + L_0*theta/T_0
	   + (C_pl - C_pv)*(theta - TEM(i)*log(TEM(i)/T_0)))) ;
    p_g[i] = pv + P_a(i) ;
    pc     = p_g[i] - P_l(i) ;
    pc0    = pc/(un-alpha*theta) ;
    sl     = SATURATION(pc0) ;
    sg     = un - sl ;
    rho_v  = pv*M_vsR/TEM(i) ;
    rho_a  = P_a(i)*M_asR/TEM(i) ;
    rho_g  = rho_v + rho_a ;
    c_v[i] = rho_v/rho_g ;
    c_a[i] = un - c_v[i] ;
    dUsdT  = zero ;
    for(j=0;j<np;j++) {
      dUsdT += w[j]*(
	       dSsdT(pc/deux*(un+a[j]),p_c3,theta,alpha,el.mat->cb[0])
	     + dSsdT(pc/deux*(un-a[j]),p_c3,theta,alpha,el.mat->cb[0])) ;
    }
    dUsdT  *= pc/deux ;
    M_l(i) = rho_l*phi*sl ;
    M_v(i) = rho_v*phi*sg ;
    M_a(i) = rho_a*phi*sg ;
    S_l(i) = C_pl*log(TEM(i)/T_0) ;
    S_v(i) = C_pv*log(TEM(i)/T_0) - log(pv/p_v0)/M_vsR + L_0/T_0 ;
    S_a(i) = C_pa*log(TEM(i)/T_0) - log(P_a(i)/p_a0)/M_asR  ;
    S(i)  = C_s*log(TEM(i)/T_0) - phi*dUsdT
            + M_l(i)*S_l(i) + M_v(i)*S_v(i) + M_a(i)*S_a(i) ;
  }
  /* flux */
  dx    = x[1][0] - x[0][0] ;
  Tij   = (TEM_n(0) + TEM_n(1))/deux ;
  S_lij = (S_ln(0) + S_ln(1))/deux ;
  S_vij = (S_vn(0) + S_vn(1))/deux ;
  S_aij = (S_an(0) + S_an(1))/deux ;
  
  W_l  = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
  W_v  = - KD_v*(p_g[1] - p_g[0])/dx - KF_v*(c_v[1] - c_v[0])/dx ;
  W_a  = - KD_a*(p_g[1] - p_g[0])/dx - KF_a*(c_a[1] - c_a[0])/dx ;
  Q   = - KTH/Tij*(TEM(1) - TEM(0))/dx 
        +   S_lij*W_l + S_vij*W_v + S_aij*W_a ;
  return(0) ;
}

int mx6(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define NEQ1      (4)
#define KE(i,j)   (ke[(i)*2*NEQ1+(j)])
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double p_v[2],p_a[2],p_g[2],dpvsdpl[2],dpvsdT[2] ;
  double c_v[2],c_a[2] ;
  double Tij,S_lij,S_vij,S_aij ;
  double pc,pc0,sl,sg,theta,dslsdpc,dslsdT,lat ;
  double rho_v,rho_a,rho_g ;
  double tr_l,trd_v,trf_v,trd_a,trf_a,trth ;
  double dx ,xm ;
  double volume[2],surf ;
  double c ;
  int    i,j,n ;
  double zero = 0.,un = 1.,deux = 2. ;
  /* Numeros d'ordre des equations et des inconnues locales */
  int E_liq = E_eau,E_vap = 3,I_p_v = 3 ;
  double ke[4*NEQ1*NEQ1] ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_vsR   = el.mat->pr[pm("M_vsR")] ;
  M_asR   = el.mat->pr[pm("M_asR")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_av0   = el.mat->pr[pm("D_av")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  lam_l   = el.mat->pr[pm("lam_l")] ;
  lam_g   = el.mat->pr[pm("lam_g")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  C_pl    = el.mat->pr[pm("C_pl")] ;
  C_pv    = el.mat->pr[pm("C_pv")] ;
  C_pa    = el.mat->pr[pm("C_pa")] ;
  L_0     = el.mat->pr[pm("L_0")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  alpha   = el.mat->pr[pm("a")] ;
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
    theta  = TEM(i) - T_0 ;
    p_v[i] = p_v0*exp(M_vsR/TEM(i)*((P_l(i) - p_l0)/rho_l + L_0*theta/T_0
	   + (C_pl - C_pv)*(theta - TEM(i)*log(TEM(i)/T_0)))) ;
    p_a[i] = P_a(i) ;
    p_g[i] = p_v[i] + P_a(i) ;
    pc     = p_g[i] - P_l(i) ;
    pc0    = pc/(un-alpha*theta) ;
    sl     = SATURATION(pc0) ;
    sg     = un - sl ;
    rho_v  = p_v[i]*M_vsR/TEM(i) ;
    rho_a  = P_a(i)*M_asR/TEM(i) ;
    rho_g  = rho_v + rho_a ;
    c_v[i] = rho_v/rho_g ;
    c_a[i] = un - c_v[i] ;
    dslsdpc = DSATURATION(pc0)/(un-alpha*theta) ;
    dslsdT  = dSsdT(pc,p_c3,theta,alpha,el.mat->cb[0]) ;
    dpvsdpl[i] = rho_v/rho_l ;
    lat        = L_0 - (C_pl - C_pv)*theta - (P_l(i) - p_l0)/rho_l ;
    dpvsdT[i]  = rho_v/TEM(i)*lat ;
    /*
      CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = +A
    */
    KE(E_liq+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*phi*rho_l*(-dslsdpc) ;
    KE(E_liq+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*rho_l*dslsdpc ;
    KE(E_liq+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*phi*rho_l*dslsdpc ;
    KE(E_liq+i*NEQ1,I_tem+i*NEQ1) += volume[i]*phi*rho_l*dslsdT ;
    /*
      CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v1) = -A
    */
    KE(E_vap+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*phi*rho_v*dslsdpc ;
    KE(E_vap+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*(M_vsR/TEM(i)*sg - rho_v*dslsdpc) ;
    KE(E_vap+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*phi*rho_v*(-dslsdpc) ;
    KE(E_vap+i*NEQ1,I_tem+i*NEQ1) += volume[i]*phi*rho_v*(-dslsdT - sg/TEM(i)) ;
    /*
      CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a1) = 0
    */
    KE(E_air+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*phi*rho_a*dslsdpc ;
    KE(E_air+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*rho_a*(-dslsdpc) ;
    KE(E_air+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*phi*(M_asR/TEM(i)*sg - rho_a*dslsdpc) ;
    KE(E_air+i*NEQ1,I_tem+i*NEQ1) += volume[i]*phi*rho_a*(-dslsdT - sg/TEM(i)) ;
    /*
      BILAN D'ENTROPIE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
    */
    KE(E_the+i*NEQ1,I_p_l+i*NEQ1) += volume[i]*(
				   - phi*alpha*pc/(un-alpha*theta)*(-dslsdpc)
				   + S_l(i)*phi*rho_l*(-dslsdpc)
				   + S_v(i)*phi*rho_v*dslsdpc
				   + S_a(i)*phi*rho_a*dslsdpc) ;
    KE(E_the+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*(
                                   - phi*alpha*pc/(un-alpha*theta)*dslsdpc
				   + S_l(i)*phi*rho_l*dslsdpc
				   + S_v(i)*phi*(M_vsR/TEM(i)*sg - rho_v*dslsdpc)
				   - M_v(i)/M_vsR/p_v[i]
				   + S_a(i)*phi*rho_a*(-dslsdpc)) ;
  KE(E_the+i*NEQ1,I_p_a+i*NEQ1) += volume[i]*(
				 - phi*alpha*pc/(un-alpha*theta)*dslsdpc
				 + S_l(i)*phi*rho_l*dslsdpc
				 + S_v(i)*phi*rho_v*(-dslsdpc)
				 + S_a(i)*phi*(M_asR/TEM(i)*sg - rho_a*dslsdpc)
				 - M_a(i)/M_asR/P_a(i)) ;
  KE(E_the+i*NEQ1,I_tem+i*NEQ1) += volume[i]*(
                                 + C_s/TEM(i)
				 - phi*alpha*pc/(un-alpha*theta)*dslsdT
				 + S_l(i)*phi*rho_l*dslsdT
				 + C_pl/TEM(i)*M_l(i)
				 + S_v(i)*phi*rho_v*(-dslsdT - sg/TEM(i))
				 + C_pv/TEM(i)*M_v(i)
				 + S_a(i)*phi*rho_a*(-dslsdT - sg/TEM(i))
				 + C_pa/TEM(i)*M_a(i)) ;
  }
  /*
    termes d'ecoulement
  */
  tr_l  = dt*surf/dx*K_l ;
  trd_v = dt*surf/dx*KD_v ;
  trd_a = dt*surf/dx*KD_a ;
  trf_v = dt*surf/dx*KF_v ;
  trf_a = dt*surf/dx*KF_a ;
  trth  = dt*surf/dx*KTH ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = +A
  */
  KE(E_liq,I_p_l)           += + tr_l ;
  KE(E_liq,I_p_l+NEQ1)      += - tr_l ;
  KE(E_liq+NEQ1,I_p_l)      += - tr_l ;
  KE(E_liq+NEQ1,I_p_l+NEQ1) += + tr_l ;
  /*
    CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v1) = -A
  */
  c = trd_v + trf_v*c_a[0]*c_v[0]/p_v[0] ;
  KE(E_vap,I_p_v)           += + c ;
  KE(E_vap+NEQ1,I_p_v)      += - c ;

  c = trd_v + trf_v*c_a[1]*c_v[1]/p_v[1] ;
  KE(E_vap,I_p_v+NEQ1)      += - c ;  
  KE(E_vap+NEQ1,I_p_v+NEQ1) += + c ;


  c = trd_v + trf_v*(-c_a[0]*c_v[0]/p_a[0]) ;
  KE(E_vap,I_p_a)           += + c ;
  KE(E_vap+NEQ1,I_p_a)      += - c ;

  c = trd_v + trf_v*(-c_a[1]*c_v[1]/p_a[1]) ;
  KE(E_vap,I_p_a+NEQ1)      += - c ;
  KE(E_vap+NEQ1,I_p_a+NEQ1) += + c ;
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a1) = 0
  */
  c = trd_a + trf_a*c_a[0]*c_v[0]/p_a[0] ;
  KE(E_air,I_p_a)           += + c ;
  KE(E_air+NEQ1,I_p_a)      += - c ;

  c = trd_a + trf_a*c_a[1]*c_v[1]/p_a[1] ;
  KE(E_air,I_p_a+NEQ1)      += - c ;  
  KE(E_air+NEQ1,I_p_a+NEQ1) += + c ;


  c = trd_a + trf_a*(-c_a[0]*c_v[0]/p_v[0]) ;
  KE(E_air,I_p_v)           += + c ;
  KE(E_air+NEQ1,I_p_v)      += - c ;

  c = trd_a + trf_a*(-c_a[1]*c_v[1]/p_v[1]) ;
  KE(E_air,I_p_v+NEQ1)      += - c ;
  KE(E_air+NEQ1,I_p_v+NEQ1) += + c ;
  /*
    BILAN D'ENTROPIE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
  */
  S_lij = (S_ln(0) + S_ln(1))/deux ;
  S_vij = (S_vn(0) + S_vn(1))/deux ;
  S_aij = (S_an(0) + S_an(1))/deux ;
  Tij   = (TEM_n(0) + TEM_n(1))/deux ;
  
  c = S_lij*tr_l ;
  KE(E_the,I_p_l)           += + c ;
  KE(E_the,I_p_l+NEQ1)      += - c ;
  KE(E_the+NEQ1,I_p_l+NEQ1) += + c ;
  KE(E_the+NEQ1,I_p_l)      += - c ;


  c = S_vij*(trd_v + trf_v*c_a[0]*c_v[0]/p_v[0])
    + S_aij*(trd_a + trf_a*(-c_a[0]*c_v[0]/p_v[0])) ;
  KE(E_the,I_p_v)           += + c ;
  KE(E_the+NEQ1,I_p_v)      += - c ;

  c = S_vij*(trd_v + trf_v*c_a[1]*c_v[1]/p_v[1])
    + S_aij*(trd_a + trf_a*(-c_a[1]*c_v[1]/p_v[1])) ;
  KE(E_the,I_p_v+NEQ1)      += - c ;
  KE(E_the+NEQ1,I_p_v+NEQ1) += + c ;


  c = S_vij*(trd_v + trf_v*(-c_a[0]*c_v[0]/p_a[0]))
    + S_aij*(trd_a + trf_a*c_a[0]*c_v[0]/p_a[0]) ;
  KE(E_the,I_p_a)           += + c ;
  KE(E_the+NEQ1,I_p_a)      += - c ;

  c = S_vij*(trd_v + trf_v*(-c_a[1]*c_v[1]/p_a[1]))
    + S_aij*(trd_a + trf_a*c_a[1]*c_v[1]/p_a[1]) ;
  KE(E_the,I_p_a+NEQ1)      += - c ;
  KE(E_the+NEQ1,I_p_a+NEQ1) += + c ;


  KE(E_the,I_tem)           += + trth/Tij ;
  KE(E_the,I_tem+NEQ1)      += - trth/Tij ;
  KE(E_the+NEQ1,I_tem+NEQ1) += + trth/Tij ;
  KE(E_the+NEQ1,I_tem)      += - trth/Tij ;

  /*
    ELIMINATION DE LA COLONNE I_p_v :
    COLONNE I_p_l += COLONNE I_p_v * dpvsdpl
    COLONNE I_tem += COLONNE I_p_v * dpvsdT
  */
  for(j=0;j<NEQ1;j++) {
    for(i=0;i<2;i++) {
      KE(j+i*NEQ1,I_p_l+i*NEQ1) += KE(j+i*NEQ1,I_p_v+i*NEQ1) * dpvsdpl[i] ;
      KE(j+i*NEQ1,I_tem+i*NEQ1) += KE(j+i*NEQ1,I_p_v+i*NEQ1) * dpvsdT[i] ;
    }
    KE(j,I_p_l+NEQ1) += KE(j,I_p_v+NEQ1) * dpvsdpl[1] ;
    KE(j+NEQ1,I_p_l) += KE(j+NEQ1,I_p_v) * dpvsdpl[0] ;
    KE(j,I_tem+NEQ1) += KE(j,I_p_v+NEQ1) * dpvsdT[1] ;
    KE(j+NEQ1,I_tem) += KE(j+NEQ1,I_p_v) * dpvsdT[0] ;
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
#undef K
#undef KE
}

void rs6(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
#define RE(n,i)   (re[(n)*NEQ1+(i)])
  double dx ,xm ;
  double volume[2],surf ;
  int    i,j ;
  double zero = 0.,un = 1.,deux = 2. ;
  /* Numeros d'ordre des equations et des inconnues locales */
  int E_liq = E_eau,E_vap = 3 ;
  enum {NEQ1 = 4} ;
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
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = +A
  */
  RE(0,E_liq) -= volume[0]*(M_l(0) - M_ln(0)) + dt*surf*W_l ;
  RE(1,E_liq) -= volume[1]*(M_l(1) - M_ln(1)) - dt*surf*W_l ;
  /*
    CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v1) = -A
  */
  RE(0,E_vap) -= volume[0]*(M_v(0) - M_vn(0)) + dt*surf*W_v ;
  RE(1,E_vap) -= volume[1]*(M_v(1) - M_vn(1)) - dt*surf*W_v ;
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a1) = 0
  */
  RE(0,E_air) -= volume[0]*(M_a(0) - M_an(0)) + dt*surf*W_a ;
  RE(1,E_air) -= volume[1]*(M_a(1) - M_an(1)) - dt*surf*W_a ;
  /*
    BILAN D'ENTROPIE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
  */
  RE(0,E_the) -= volume[0]*(S(0) - S_n(0)) + dt*surf*Q ;
  RE(1,E_the) -= volume[1]*(S(1) - S_n(1)) - dt*surf*Q ;
  /*
    ELIMINATION DE LA LIGNE E_vap : LIGNE E_liq += LIGNE E_vap
  */
  for(i=0;i<2;i++) RE(i,E_liq) += RE(i,E_vap) ;
  /*
    ASSEMBLAGE DANS R
  */
  for(j=0;j<NEQ;j++) for(i=0;i<2;i++) R(i,j) += RE(i,j) ;

#undef R
#undef RE
}

int so6(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double p_v[2],p_g[2] ;
  double c_v[2],c_a[2] ;
  double pl,pv,pa,pg,pc,theta,sl,tem ;
  double rho_v,rho_a,rho_g ;
  double w_l,w_v,w_a,w_g,q_T ;
  double dx ;
  int    i,j,nso ;                                                           
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0.,un = 1. ;

  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_vsR   = el.mat->pr[pm("M_vsR")] ;
  M_asR   = el.mat->pr[pm("M_asR")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_v0    = el.mat->pr[pm("p_v0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_av0   = el.mat->pr[pm("D_av")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  lam_l   = el.mat->pr[pm("lam_l")] ;
  lam_g   = el.mat->pr[pm("lam_g")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  C_pl    = el.mat->pr[pm("C_pl")] ;
  C_pv    = el.mat->pr[pm("C_pv")] ;
  C_pa    = el.mat->pr[pm("C_pa")] ;
  L_0     = el.mat->pr[pm("L_0")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  alpha   = el.mat->pr[pm("a")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /* initialisation */
  nso = 10 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pressions */
  pl =  param(u,h_s,el.nn,I_p_l) ;
  pa =  param(u,h_s,el.nn,I_p_a) ;
  /* temperature */
  tem   =  param(u,h_s,el.nn,I_tem) ;
  theta = tem - T_0 ;
  /* autres pressions */
  pv = p_v0*exp(M_vsR/tem*((pl - p_l0)/rho_l + L_0*theta/T_0
			   + (C_pl - C_pv)*(theta - tem*log(tem/T_0)))) ;
  pg = pv + pa ;
  pc = pg - pl ;
  /* saturation */
  sl = SATURATION(pc/(un-alpha*theta)) ;

  /* Calcul de certaines quantites aux noeuds */
  for(i=0;i<2;i++) {
    theta  = TEM(i) - T_0 ;
    p_v[i] = p_v0*exp(M_vsR/TEM(i)*((P_l(i) - p_l0)/rho_l + L_0*theta/T_0
           + (C_pl - C_pv)*(theta - TEM(i)*log(TEM(i)/T_0)))) ;
    p_g[i] = p_v[i] + P_a(i) ;
    rho_v  = p_v[i]*M_vsR/TEM(i) ;
    rho_a  = P_a(i)*M_asR/TEM(i) ;
    rho_g  = rho_v + rho_a ;
    c_v[i] = rho_v/rho_g ;
    c_a[i] = un - c_v[i] ;
  }
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  w_l  = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
  w_g  = - (KD_v+KD_a)*(p_g[1] - p_g[0])/dx ;
  w_v  = - KD_v*(p_g[1] - p_g[0])/dx - KF_v*(c_v[1] - c_v[0])/dx ;
  w_a  = - KD_a*(p_g[1] - p_g[0])/dx - KF_a*(c_a[1] - c_a[0])/dx ;
  q_T  = - KTH*(TEM(1) - TEM(0))/dx ;

  /* quantites exploitees */
  strcpy(r[0].text,"p_l") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"p_a") ; r[1].n = 1 ;
  r[1].v[0] = pa ;
  strcpy(r[2].text,"T") ; r[2].n = 1 ;
  r[2].v[0] = tem ;
  strcpy(r[3].text,"w_l") ; r[3].n = 1 ;
  r[3].v[0] = w_l ;
  strcpy(r[4].text,"w_g") ; r[4].n = 1 ;
  r[4].v[0] = w_g ;
  strcpy(r[5].text,"w_v") ; r[5].n = 1 ;
  r[5].v[0] = w_v ;
  strcpy(r[6].text,"w_a") ; r[6].n = 1 ;
  r[6].v[0] = w_a ;
  strcpy(r[7].text,"q_T") ; r[7].n = 1 ;
  r[7].v[0] = q_T ;
  strcpy(r[8].text,"saturation") ; r[8].n = 1 ;
  r[8].v[0] = sl ;
  strcpy(r[9].text,"p_v") ; r[9].n = 1 ;
  r[9].v[0] = pv ;
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

double dSsdT(double pc,double pc3,double theta,double alpha,crbe_t cb)
{
  double dsl,pc0 ;
  double un = 1. ;
  
  pc0  = pc/(un-alpha*theta) ;
  dsl  = dsaturation(pc0,pc3,cb)*pc0*alpha/(un-alpha*theta) ;
  return(dsl) ;
}

#undef NEQ
#undef NVI
#undef NVE

#undef P_l
#undef P_a
#undef TEM

#undef M_l
#undef M_v
#undef M_a
#undef S
#undef S_l
#undef S_v
#undef S_a
#undef W_l
#undef W_v
#undef W_a
#undef Q

#undef P_ln
#undef P_an
#undef TEM_n

#undef M_ln
#undef M_vn
#undef M_an
#undef S_n
#undef S_ln
#undef S_vn
#undef S_an
#undef W_ln
#undef W_vn
#undef W_an
#undef Q_n

#undef K_l
#undef KD_v
#undef KF_v
#undef KD_a
#undef KF_a
#undef KTH

#undef SATURATION
#undef DSATURATION
